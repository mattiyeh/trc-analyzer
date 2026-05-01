#!/usr/bin/env python3
"""End-to-end SigProfiler workflow for the TRC pan-cancer study.

Setup (1-3):
   1.  Convert each ICGC TSV (SBS + indel) -> minimal MAF.
   2.  Run SigProfilerMatrixGenerator on every (tumor x group) -> SBS96
       and ID83 matrices.
   3.  Aggregate the 17 tumor-type "promoter_high" SBS96 + ID83 matrices
       into the pan-cancer PPP-HTG cohort.

Pan-cancer Extractors -- compartment specificity (4-7):
   4.  Pan-cancer PPP-HTG Extractor (canonical): de novo NMF on the
       pan-cancer PPP-HTG SBS96 + ID83 aggregates.
   5.  Pan-cancer PPP-LTG Extractor (PPP-HTG specificity control).
   6.  Pan-cancer non-promoter SBS96 Extractor (PPP specificity control).
   7.  Pan-cancer CFS Extractor (unified mechanism, all 17 tumors).

Pan-cancer Extractors -- compositional bias (8-10):
   8.  Liver-excluded pan-cancer Extractor.
   9.  Pancreas-excluded pan-cancer Extractor.
  10.  Equal-weight pan-cancer Extractor, cap 50 donors per tumor seed=42.

Assignment (11-12):
  11.  Constrained SigProfilerAssignment (SBS + ID83) on every
       (tumor x group) using the step-4 Extractor-derived COSMIC signatures.
  12.  Unconstrained Assignment (full COSMIC database) on the pan-cancer
       PPP-HTG aggregates -- the constrained-vs-unconstrained sensitivity
       comparison; required by step 14.

Pan-cancer Extractors -- hypermutator subsets (13-15):
  13.  No-hyper count-based: drop samples > mean+2SD of PPP-HTG SBS count
       and rerun Extractor (both kinds).
  14.  No-MMR/POLE mechanism-based: drop donors with MMR >=30% OR POLE
       >=30% attribution in the unconstrained SBS96 Assignment (step 12)
       and rerun Extractor (PCAWG-consistent alternative to step 13).
  15.  No-WGS-hyper: drop donors whose whole-genome unfiltered SBS96 count
       > mean+2SD over the PPP-HTG cohort. Single donor exclusion set
       applied to both SBS96 and ID83.

KO validation (22):
  22.  Per-subclone constrained Assignment of Zou 2021 (Nik-Zainal) KO
       iPSC subclones (RECQL5, SETX, ATP2B4 + comparator panel) against
       COSMIC v3.5 + SBS96-TRC + SBS96B. Tests whether knocking out a
       TRC-resolution helicase recovers the SBS96-TRC signature in vitro
       -- the mechanistic-causation half of the validation logic.
       Standalone runner: validation/run_validation.py.

Per-tumor + subset analyses (16-17):
  16.  Per-tumor SBS96 Extractor on each tumor's promoter_high and
       promoter_low matrices (PPP-HTG vs PPP-LTG validation).
  17.  Ovary HRD exclusion: remove ovary donors with SBS3 attribution > 50%
       in the per-tumor ovary Extractor's COSMIC-decomposed Activities
       (produced by step 16); rerun MatrixGenerator + Extractor on the
       filtered ovary cohort.

Resumable: each step skips work whose outputs already exist.
"""

import csv
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


# --- Configuration ----------------------------------------------------------

REPO_ROOT     = Path("/data/research/projects/trc_signatures")
ICGC_OUTPUTS  = REPO_ROOT / "outputs/java/2026.04.27"
TIMESTAMP_TAG = "2026.04.27_23.34.10"
WORKDIR       = REPO_ROOT / "outputs/sigprofiler"
REFERENCE     = "GRCh37"
COSMIC_VER    = 3.5

TUMORS = [
    "bladder", "bone", "breast", "cervix", "colorectal", "esophagus",
    "gallbladder", "kidney", "liver", "lung", "melanoma", "ovary",
    "pancreas", "pnet", "prostate", "stomach", "uterus",
]

# (input subdir under each DCO__* dir, file suffix on disk, group label)
GROUPS = [
    ("promoter/sbs", "promoter_sbs_mutations_high", "promoter_high"),
    ("promoter/sbs", "promoter_sbs_mutations_low",  "promoter_low"),
    ("promoter/sbs", "promoter_sbs_mutations_mid",  "promoter_mid"),
    ("promoter/sbs", "promoter_sbs_mutations_zero", "promoter_zero"),
    ("promoter/sbs", "promoter_sbs_mutations",      "promoter_all"),
    ("promoter/sbs", "non_promoter_sbs_mutations",  "non_promoter"),
    ("cfs/sbs",      "cfs_sbs_mutations",           "cfs"),
    ("cfs/sbs",      "non_cfs_sbs_mutations",       "non_cfs"),
]

INDEL_GROUPS = [
    ("promoter/indel", "promoter_indel_mutations_high", "promoter_high"),
    ("promoter/indel", "promoter_indel_mutations_low",  "promoter_low"),
    ("promoter/indel", "promoter_indel_mutations_mid",  "promoter_mid"),
    ("promoter/indel", "promoter_indel_mutations_zero", "promoter_zero"),
    ("promoter/indel", "promoter_indel_mutations",      "promoter_all"),
    ("promoter/indel", "non_promoter_indel_mutations",  "non_promoter"),
    ("cfs/indel",      "cfs_indel_mutations",           "cfs"),
    ("cfs/indel",      "non_cfs_indel_mutations",       "non_cfs"),
]

PRIMARY_GROUP         = "promoter_high"
PANCANCER_LABEL       = "pancancer_promoter_high"
PANCANCER_INDEL_LABEL = "pancancer_promoter_high_indel"

# Per-tumor extraction is run on these groups (SBS96 only).
PER_TUMOR_GROUPS = ["promoter_high", "promoter_low"]

# Ovary HRD exclusion threshold: donors with SBS3 fraction above this are
# considered HRD-dominant and removed for the ovary sensitivity analysis.
HRD_SBS3_THRESHOLD = 0.50

# Equal-weight (step 10) compositional rebalancing parameters.
EQUAL_WEIGHT_CAP  = 50
EQUAL_WEIGHT_SEED = 42

# Tumors to exclude in single-tumor sensitivity steps.
LIVER_EXCLUDED_TUMORS    = ["liver"]
PANCREAS_EXCLUDED_TUMORS = ["pancreas"]

# Which de novo column in the canonical pan-cancer Extractor output is
# the TRC signature. SigProfiler's NMF output assigns "A"/"B" labels by
# proportional mass, so the mapping can swap if the canonical run is
# rerun (it did from 2023 to 2026 when SBS96-TRC mass crossed 50%).
# These constants reflect the current canonical run -- inspect
# docs/EXTRACTOR_RESULTS_2026.md and docs/ID83_EXTRACTOR_RESULTS_2026.md
# if the canonical pan-cancer extraction is ever rerun.
PANCANCER_TRC_COLUMN = {
    "SBS96": "SBS96A",   # SBS96-TRC; T[C>G]T-dominant. See EXTRACTOR_RESULTS_2026.md
    "ID83":  "ID83B",    # ID83-TRC; ≥5 bp microhomology deletion-dominant.
}

# Minimal MAF columns understood by SigProfilerMatrixGenerator.
MAF_COLUMNS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
    "Chromosome", "Start_position", "End_position", "Strand",
    "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
]

# Per-context dispatch: file paths, parameter names, and output conventions
# differ between SBS96 and ID83.
KIND = {
    "SBS96": {
        "matgen_subdir":       "SBS",
        "matrix_ext":          "SBS96.all",
        "extractor_subdir":    "SBS96",
        "extractor_kdir_pre":  "SBS96",
        "cosmic_subdir":       "COSMIC_SBS96_Decomposed_Solution",
        "cosmic_filename":     "COSMIC_SBS96_Signatures.txt",
        "context_type":        "SBS96",
        "assign_context_type": "96",
        "collapse_to_SBS96":   True,
        "n_channels":          96,
    },
    "ID83": {
        "matgen_subdir":       "ID",
        "matrix_ext":          "ID83.all",
        "extractor_subdir":    "ID83",
        "extractor_kdir_pre":  "ID83",
        "cosmic_subdir":       "COSMIC_ID83_Decomposed_Solution",
        "cosmic_filename":     "COSMIC_ID83_Signatures.txt",
        "context_type":        "ID83",
        "assign_context_type": "ID",
        "collapse_to_SBS96":   False,   # ID Assignment crashes if True
        "n_channels":          83,
    },
}

# Extractor parameters per the best-practices guide. KL divergence is the
# default objective in SigProfilerExtractor's underlying NMF and is not
# exposed as a top-level parameter. context_type is set per-call.
EXTRACTOR_KW_BASE = dict(
    reference_genome=REFERENCE,
    opportunity_genome=REFERENCE,
    minimum_signatures=1,
    maximum_signatures=10,
    nmf_replicates=100,
    nmf_init="random",
    matrix_normalization="gmm",
    precision="single",
    resample=True,
    cosmic_version=COSMIC_VER,
    cpu=-1,
)


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


# --- Sanity gates -----------------------------------------------------------
# After-Step-3 and after-Step-11 byte-identity checks for the post-Java-rerun
# (2026.04.27) workflow. The pipeline auto-aborts if a freshly-rebuilt file
# diverges from a saved backup, so the long unattended run is safe.
#
# To enable a gate, drop a backup file at the matching path in
# WORKDIR/_sanity_backup/ before launching the pipeline. With no backup
# present, the gate logs "no backup, skipping" and continues -- so this
# infrastructure is invisible to normal future runs.
#
# Expected backup files (post-Java-rerun workflow):
#   _sanity_backup/pancancer_promoter_high.SBS96.all
#   _sanity_backup/pancancer_promoter_high_indel.ID83.all
#   _sanity_backup/bone_promoter_high_Assignment_Solution_Activities.txt

SANITY_BACKUP_DIR = None  # set after WORKDIR is known

def _sanity_backup_dir():
    return WORKDIR / "_sanity_backup"

def _sanity_check_pair(label, current_path, backup_path):
    if not Path(backup_path).exists():
        log(f"  sanity [{label}]: no backup at {backup_path.name}, skipping")
        return
    if not Path(current_path).exists():
        sys.exit(f"ABORT sanity [{label}]: rebuilt file missing at {current_path}")
    import filecmp
    if not filecmp.cmp(str(current_path), str(backup_path), shallow=False):
        sys.exit(
            f"ABORT sanity [{label}]: rebuilt file does NOT match backup.\n"
            f"  current: {current_path}\n"
            f"  backup:  {backup_path}\n"
            f"Investigate before re-running. Downstream steps not executed."
        )
    log(f"  sanity [{label}]: PASS (byte-identical to backup)")


def _sanity_check_after_step3():
    log("Sanity gate after Step 3: pan-cancer aggregates vs backup")
    bdir = _sanity_backup_dir()
    _sanity_check_pair(
        "aggregate-SBS96",
        WORKDIR / "pancancer" / "pancancer_promoter_high"
            / "pancancer_promoter_high.SBS96.all",
        bdir / "pancancer_promoter_high.SBS96.all",
    )
    _sanity_check_pair(
        "aggregate-ID83",
        WORKDIR / "pancancer" / "pancancer_promoter_high_indel"
            / "pancancer_promoter_high_indel.ID83.all",
        bdir / "pancancer_promoter_high_indel.ID83.all",
    )


def _sanity_check_after_step11():
    log("Sanity gate after Step 11: bone PPP-HTG Assignment vs backup")
    bdir = _sanity_backup_dir()
    _sanity_check_pair(
        "bone-PPP-HTG-Assignment",
        WORKDIR / "assignment" / "bone_promoter_high"
            / "Assignment_Solution" / "Activities"
            / "Assignment_Solution_Activities.txt",
        bdir / "bone_promoter_high_Assignment_Solution_Activities.txt",
    )


# --- Step 1: TSV -> MAF -----------------------------------------------------

def icgc_tsv_path(tumor, subdir, suffix):
    return (
        ICGC_OUTPUTS
        / f"DCO__{TIMESTAMP_TAG}_{tumor}"
        / subdir
        / f"{tumor}_{suffix}.tsv"
    )


def sbs_maf_path(project):
    return WORKDIR / "maf" / project / f"{project}.maf"


def indel_maf_path(project):
    return WORKDIR / "maf" / project / f"{project}_indel.maf"


def convert_sbs_tsv_to_maf(tsv_path, maf_path):
    """Convert one ICGC SBS TSV to a minimal MAF. Returns row count."""
    with open(tsv_path) as fin, open(maf_path, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        ix = {c: i for i, c in enumerate(header)}
        writer.writerow(MAF_COLUMNS)
        n = 0
        for row in reader:
            if not row:
                continue
            chrom = row[ix["chromosome"]]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            ref    = row[ix["reference_genome_allele"]]
            alt    = row[ix["mutated_to_allele"]]
            sample = row[ix["icgc_donor_id"]]
            start  = row[ix["chromosome_start"]]
            end    = row[ix["chromosome_end"]]
            writer.writerow([
                "Unknown", "0", "ICGC", REFERENCE,
                chrom, start, end, "+",
                "Unknown", "SNP",
                ref, ref, alt, "", "",
                sample,
            ])
            n += 1
    return n


def convert_indel_tsv_to_maf(tsv_path, maf_path):
    """Convert one ICGC indel TSV to a standard MAF.

    ICGC indel coordinates are already MAF-compatible:
      - Deletion:  ref = deleted bases, alt = "-",
                   chromosome_start = first deleted base (1-based).
      - Insertion: ref = "-", alt = inserted bases,
                   chromosome_start = base before insertion site.
    SigProfilerMatrixGenerator's MAF parser handles indels internally
    (it decrements Start_position by 1 to convert to its pre-base
    simple-file convention), so we pass coordinates and alleles through
    verbatim. Adding a leading '-' to ref/alt or pre-decrementing the
    start is what the simple-format converter wants, not the MAF parser
    -- doing that here causes matgen to silently report 0 indels.

    Returns row count.
    """
    with open(tsv_path) as fin, open(maf_path, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        ix = {c: i for i, c in enumerate(header)}
        writer.writerow(MAF_COLUMNS)
        n = 0
        for row in reader:
            if not row:
                continue
            chrom = row[ix["chromosome"]]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            ref    = row[ix["reference_genome_allele"]]
            alt    = row[ix["mutated_to_allele"]]
            sample = row[ix["icgc_donor_id"]]
            start  = row[ix["chromosome_start"]]
            end    = row[ix["chromosome_end"]]
            if alt == "-" and ref != "-":
                vt = "DEL"
            elif ref == "-" and alt != "-":
                vt = "INS"
            else:
                # Not an indel (or malformed row); skip.
                continue
            writer.writerow([
                "Unknown", "0", "ICGC", REFERENCE,
                chrom, start, end, "+",
                "Unknown", vt,
                ref, ref, alt, "", "",
                sample,
            ])
            n += 1
    return n


def step1_convert_all():
    log("Step 1: TSV -> MAF (SBS + indel)")
    for tumor in TUMORS:
        for subdir, suffix, label in GROUPS:
            tsv = icgc_tsv_path(tumor, subdir, suffix)
            project = f"{tumor}_{label}"
            maf = sbs_maf_path(project)
            if maf.exists():
                continue
            if not tsv.exists():
                log(f"  MISSING SBS TSV: {tsv}")
                continue
            maf.parent.mkdir(parents=True, exist_ok=True)
            n = convert_sbs_tsv_to_maf(tsv, maf)
            log(f"  SBS  {project}: {n} rows -> {maf.name}")
        for subdir, suffix, label in INDEL_GROUPS:
            tsv = icgc_tsv_path(tumor, subdir, suffix)
            project = f"{tumor}_{label}"
            maf = indel_maf_path(project)
            if maf.exists():
                continue
            if not tsv.exists():
                log(f"  MISSING indel TSV: {tsv}")
                continue
            maf.parent.mkdir(parents=True, exist_ok=True)
            n = convert_indel_tsv_to_maf(tsv, maf)
            log(f"  ID   {project}: {n} rows -> {maf.name}")


# --- Step 2: SigProfilerMatrixGenerator -------------------------------------

def matgen_matrix_path(project, kind):
    info = KIND[kind]
    return (
        WORKDIR / "maf" / project
        / "output" / info["matgen_subdir"]
        / f"{project}.{info['matrix_ext']}"
    )


def step2_matrix_generator():
    log("Step 2: SigProfilerMatrixGenerator (SBS96 + ID83)")
    import shutil
    from SigProfilerMatrixGenerator.scripts import (
        SigProfilerMatrixGeneratorFunc as matGen,
    )
    # Build the union of all (tumor, group) pairs from SBS + indel groups.
    all_labels = sorted({l for _, _, l in GROUPS} | {l for _, _, l in INDEL_GROUPS})
    for tumor in TUMORS:
        for label in all_labels:
            project = f"{tumor}_{label}"
            project_dir = WORKDIR / "maf" / project
            sbs_maf = sbs_maf_path(project)
            id_maf  = indel_maf_path(project)
            sbs96   = matgen_matrix_path(project, "SBS96")
            id83    = matgen_matrix_path(project, "ID83")

            need_sbs = sbs_maf.exists() and not sbs96.exists()
            need_id  = id_maf.exists()  and not id83.exists()
            if not (need_sbs or need_id):
                continue
            if not (sbs_maf.exists() or id_maf.exists()):
                continue
            # matgen copies inputs to <project_dir>/input/ on the first run
            # and reuses that cache on subsequent runs -- so a MAF added
            # later is ignored unless we clear the cache first.
            cache = project_dir / "input"
            if cache.exists():
                shutil.rmtree(cache)
            log(f"  matgen {project} (need_sbs={need_sbs}, need_id={need_id})")
            try:
                matGen.SigProfilerMatrixGeneratorFunc(
                    project, REFERENCE, str(project_dir),
                    plot=False, exome=False, bed_file=None,
                    chrom_based=False, tsb_stat=False, seqInfo=False,
                )
            except Exception as exc:
                log(f"    SKIP {project}: {exc}")


# --- Step 3: aggregate pan-cancer high (SBS96 + ID83) -----------------------

def pancancer_matrix_path(kind):
    label = PANCANCER_LABEL if kind == "SBS96" else PANCANCER_INDEL_LABEL
    return WORKDIR / "pancancer" / label / f"{label}.{KIND[kind]['matrix_ext']}"


def aggregate_high(kind):
    label = PANCANCER_LABEL if kind == "SBS96" else PANCANCER_INDEL_LABEL
    out_path = pancancer_matrix_path(kind)
    if out_path.exists():
        log(f"  [{kind}] exists: {out_path}")
        return out_path

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for tumor in TUMORS:
        project = f"{tumor}_{PRIMARY_GROUP}"
        m = matgen_matrix_path(project, kind)
        if not m.exists():
            log(f"  [{kind}] missing matrix for {project} -- skipping")
            continue
        df = pd.read_csv(m, sep="\t", index_col=0)
        if df.shape[1] == 0:
            log(f"  [{kind}] empty matrix for {project}")
            continue
        df.columns = [f"{tumor}__{c}" for c in df.columns]
        frames.append((tumor, df))
        log(f"  [{kind}] {project}: {df.shape[1]} samples")

    if not frames:
        log(f"  [{kind}] no high matrices found -- skipping aggregation")
        return None

    # All matrices must share the same channel index in the same order;
    # otherwise pd.concat would silently misalign rows.
    ref_tumor, ref_df = frames[0]
    for tumor, df in frames[1:]:
        if not df.index.equals(ref_df.index):
            sys.exit(
                f"{kind} channel mismatch: '{tumor}' row index differs from "
                f"'{ref_tumor}'. Refusing to aggregate -- inputs were not "
                f"produced by a consistent MatrixGenerator run."
            )

    agg = pd.concat([df for _, df in frames], axis=1).fillna(0).astype(int)
    agg.index.name = "MutationType"
    agg.to_csv(out_path, sep="\t")
    log(f"  [{kind}] wrote {out_path} "
        f"({agg.shape[1]} samples, {agg.shape[0]} channels)")
    return out_path


def step3_aggregate_high():
    log(f"Step 3: aggregate '{PRIMARY_GROUP}' across {len(TUMORS)} tumor types")
    sbs_path = aggregate_high("SBS96")
    id_path  = aggregate_high("ID83")
    return sbs_path, id_path


# --- Generalized aggregator for sensitivity-step pan-cancer matrices -------

def _aggregate_group(group_label, kind, out_path, tumors=None):
    """Aggregate per-tumor matrices for one (group, kind) combination.

    Mirrors aggregate_high() but parameterized over the per-tumor
    matrix subdir, so it can serve promoter_low / non_promoter / cfs.
    `tumors` defaults to the full TUMORS list. Skips tumors whose
    per-tumor matrix is missing or empty (header-only).
    """
    if out_path.exists():
        log(f"  [{kind}] exists: {out_path}")
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    iter_tumors = TUMORS if tumors is None else tumors
    frames = []
    for tumor in iter_tumors:
        project = f"{tumor}_{group_label}"
        m = matgen_matrix_path(project, kind)
        if not m.exists():
            log(f"  [{kind}] missing matrix for {project} -- skipping")
            continue
        df = pd.read_csv(m, sep="\t", index_col=0)
        if df.shape[1] == 0:
            log(f"  [{kind}] empty matrix for {project}")
            continue
        df.columns = [f"{tumor}__{c}" for c in df.columns]
        frames.append((tumor, df))
        log(f"  [{kind}] {project}: {df.shape[1]} samples")

    if not frames:
        log(f"  [{kind}] no '{group_label}' matrices found -- skipping aggregation")
        return None

    ref_tumor, ref_df = frames[0]
    for tumor, df in frames[1:]:
        if not df.index.equals(ref_df.index):
            sys.exit(
                f"{kind} channel mismatch: '{tumor}' row index differs from "
                f"'{ref_tumor}' for group '{group_label}'. Refusing to aggregate."
            )

    agg = pd.concat([df for _, df in frames], axis=1).fillna(0).astype(int)
    agg.index.name = "MutationType"
    agg.to_csv(out_path, sep="\t")
    log(f"  [{kind}] wrote {out_path} "
        f"({agg.shape[1]} samples, {agg.shape[0]} channels)")
    return out_path


def aggregate_low(kind):
    out_path = (WORKDIR / "sensitivity" / "pancancer_low" / "maf"
                / f"pancancer_promoter_low.{KIND[kind]['matrix_ext']}")
    return _aggregate_group("promoter_low", kind, out_path)


def aggregate_nonpromoter(kind):
    out_path = (WORKDIR / "sensitivity" / "pancancer_nonpromoter" / "maf"
                / f"pancancer_non_promoter.{KIND[kind]['matrix_ext']}")
    return _aggregate_group("non_promoter", kind, out_path)


def aggregate_cfs(kind):
    """All 17 tumors -- CFS regions don't depend on expression data."""
    out_path = (WORKDIR / "sensitivity" / "pancancer_cfs" / "maf"
                / f"pancancer_cfs.{KIND[kind]['matrix_ext']}")
    return _aggregate_group("cfs", kind, out_path)


def make_excluded_aggregate(kind, tumors_to_exclude, out_dir, out_label):
    """Drop columns whose names start with any of `<tumor>__` from the
    canonical pan-cancer aggregate matrix and write the result.
    Returns the new matrix path or None if the source matrix is missing.
    """
    src = pancancer_matrix_path(kind)
    out_path = out_dir / "maf" / f"{out_label}.{KIND[kind]['matrix_ext']}"
    if out_path.exists():
        log(f"  [{kind}] exists: {out_path}")
        return out_path
    if not src.exists():
        log(f"  [{kind}] ABORT: source pan-cancer matrix missing: {src}")
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(src, sep="\t", index_col=0)
    n_in = df.shape[1]
    prefixes = tuple(f"{t}__" for t in tumors_to_exclude)
    keep = [c for c in df.columns if not c.startswith(prefixes)]
    drop = sorted(set(df.columns) - set(keep))
    log(f"  [{kind}] excluding tumors {tumors_to_exclude}: "
        f"dropping {len(drop)}/{n_in} samples")
    df[keep].to_csv(out_path, sep="\t")
    log(f"  [{kind}] wrote {out_path} ({len(keep)} samples)")
    return out_path


def make_capped_aggregate(kind, cap, seed, out_dir, out_label):
    """Subsample columns of the canonical pan-cancer aggregate per tumor:
    if a tumor has > `cap` donors, randomly sample `cap` with the given
    `seed`. Tumors with <= `cap` donors are kept whole.
    """
    import numpy as np
    src = pancancer_matrix_path(kind)
    out_path = out_dir / "maf" / f"{out_label}.{KIND[kind]['matrix_ext']}"
    if out_path.exists():
        log(f"  [{kind}] exists: {out_path}")
        return out_path
    if not src.exists():
        log(f"  [{kind}] ABORT: source pan-cancer matrix missing: {src}")
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(src, sep="\t", index_col=0)

    # Group columns by tumor prefix.
    by_tumor = {}
    for c in df.columns:
        tumor = c.split("__", 1)[0]
        by_tumor.setdefault(tumor, []).append(c)

    rng = np.random.default_rng(seed)
    keep = []
    for tumor in sorted(by_tumor):
        cols = by_tumor[tumor]
        if len(cols) > cap:
            # Sort first for deterministic ordering, then sample.
            cols_sorted = sorted(cols)
            picked = list(rng.choice(cols_sorted, size=cap, replace=False))
            log(f"  [{kind}] {tumor}: capping {len(cols)} -> {cap}")
        else:
            picked = sorted(cols)
            log(f"  [{kind}] {tumor}: keeping all {len(cols)}")
        keep.extend(sorted(picked))

    df[keep].to_csv(out_path, sep="\t")
    log(f"  [{kind}] wrote {out_path} ({len(keep)} samples, "
        f"cap={cap}, seed={seed})")
    return out_path


def make_donor_excluded_aggregate(kind, donors_to_exclude, out_dir, out_label):
    """Drop specific sample columns (full pan-cancer column names like
    'ovary__DO35442') from the canonical pan-cancer aggregate matrix.
    Returns the new matrix path or None if the source matrix is missing.
    """
    src = pancancer_matrix_path(kind)
    out_path = out_dir / "maf" / f"{out_label}.{KIND[kind]['matrix_ext']}"
    if out_path.exists():
        log(f"  [{kind}] exists: {out_path}")
        return out_path
    if not src.exists():
        log(f"  [{kind}] ABORT: source pan-cancer matrix missing: {src}")
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(src, sep="\t", index_col=0)
    n_in = df.shape[1]
    drop_set = set(donors_to_exclude)
    keep = [c for c in df.columns if c not in drop_set]
    matched = sorted(set(df.columns) & drop_set)
    requested = len(drop_set)
    log(f"  [{kind}] excluding {len(matched)}/{requested} requested donors: "
        f"{len(keep)}/{n_in} samples remain")
    df[keep].to_csv(out_path, sep="\t")
    log(f"  [{kind}] wrote {out_path}")
    return out_path


# --- TRC component identification & cosine reporting ------------------------

def _pancancer_de_novo_signatures_path(kind):
    """Canonical pan-cancer Extractor de novo signatures file.
    Used as the cosine reference for sensitivity-step extracted sigs.
    """
    label = PANCANCER_LABEL if kind == "SBS96" else PANCANCER_INDEL_LABEL
    sub = KIND[kind]["extractor_subdir"]
    return (WORKDIR / "pancancer" / label / "extractor" / sub
            / "Suggested_Solution"
            / f"{sub}_De-Novo_Solution" / "Signatures"
            / f"{sub}_De-Novo_Signatures.txt")


def _sensitivity_de_novo_signatures_path(extractor_out, kind):
    sub = KIND[kind]["extractor_subdir"]
    return (Path(extractor_out) / sub / "Suggested_Solution"
            / f"{sub}_De-Novo_Solution" / "Signatures"
            / f"{sub}_De-Novo_Signatures.txt")


def _cosine(a, b):
    import numpy as np
    a = a.astype(float).to_numpy()
    b = b.astype(float).to_numpy()
    na, nb = (a * a).sum() ** 0.5, (b * b).sum() ** 0.5
    if na == 0 or nb == 0:
        return float("nan")
    return float((a * b).sum() / (na * nb))


def report_trc_cosine(extractor_out, label, kind):
    """Log a pairwise cosine matrix between the sensitivity-step de novo
    signatures and the canonical pan-cancer de novo signatures, and
    identify the TRC component (sensitivity sig with highest cosine to
    PANCANCER_TRC_COLUMN[kind]).

    Defensive: silently skips if either file is missing (e.g. canonical
    pan-cancer Extractor hasn't run yet).
    """
    sens_path = _sensitivity_de_novo_signatures_path(extractor_out, kind)
    ref_path  = _pancancer_de_novo_signatures_path(kind)
    if not sens_path.exists():
        log(f"  [{kind}] TRC cosine: sensitivity sig file missing ({sens_path})")
        return
    if not ref_path.exists():
        log(f"  [{kind}] TRC cosine: canonical pan-cancer ref missing ({ref_path})")
        return

    sens = pd.read_csv(sens_path, sep="\t", index_col=0)
    ref  = pd.read_csv(ref_path,  sep="\t", index_col=0)
    if not sens.index.equals(ref.index):
        # Defensive: same channel index expected.
        log(f"  [{kind}] TRC cosine: channel index mismatch -- skipping")
        return

    trc_col = PANCANCER_TRC_COLUMN[kind]
    if trc_col not in ref.columns:
        log(f"  [{kind}] TRC cosine: '{trc_col}' not in canonical ref columns "
            f"{list(ref.columns)}; PANCANCER_TRC_COLUMN may need updating")
        return

    log(f"  [{kind}] TRC cosines for {label} (vs canonical pan-cancer):")
    log(f"    {'sens_sig':<14} | " + " | ".join(f"{c:^8}" for c in ref.columns))
    log(f"    {'-'*14}-+-" + "-+-".join("-" * 8 for _ in ref.columns))
    best_sens, best_cos = None, -2.0
    for sens_col in sens.columns:
        cosines = {ref_col: _cosine(sens[sens_col], ref[ref_col])
                   for ref_col in ref.columns}
        log(f"    {sens_col:<14} | " + " | ".join(
            f"{cosines[c]:8.4f}" for c in ref.columns))
        if cosines[trc_col] > best_cos:
            best_cos = cosines[trc_col]
            best_sens = sens_col
    log(f"  [{kind}] TRC component for {label}: {best_sens} "
        f"(cosine to canonical {trc_col} = {best_cos:.4f})")


# --- Step 4 helpers: Extractor de novo (used by step 4 + sensitivity steps) -

def cosmic_signatures_path(extractor_output, kind):
    info = KIND[kind]
    return (
        Path(extractor_output)
        / info["extractor_subdir"] / "Suggested_Solution"
        / info["cosmic_subdir"] / "Signatures"
        / info["cosmic_filename"]
    )


def run_extractor(matrix_path, output_dir, label, kind):
    log(f"Extractor [{kind}]: de novo on {label}")
    if matrix_path is None or not Path(matrix_path).exists():
        log(f"  ABORT: matrix missing for {label}")
        return None
    cosmic_db = cosmic_signatures_path(output_dir, kind)
    if cosmic_db.exists():
        log(f"  exists: {cosmic_db}")
        return cosmic_db
    from SigProfilerExtractor import sigpro as sig
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    kw = dict(EXTRACTOR_KW_BASE, context_type=KIND[kind]["context_type"])
    log(f"  params: {kw}")
    try:
        sig.sigProfilerExtractor(
            "matrix", str(output_dir), str(matrix_path),
            **kw,
        )
    except Exception as exc:
        log(f"  Extractor FAILED for {label}: {exc}")
        return None
    return cosmic_db


# --- Real per-signature minimum silhouette ---------------------------------
# The "Minimum Stability" column in All_solutions_stat.csv is mislabeled:
# it is actually the average across signatures. The true per-signature minimum
# must be recovered from the per-k signature-level files.

def _k_from_dirname(path):
    try:
        return int(path.name.split("_")[1])
    except (ValueError, IndexError):
        return None


def _scan_min_silhouette(k_dir, k):
    best = None
    for p in k_dir.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in (".csv", ".txt", ".tsv"):
            continue
        try:
            df = pd.read_csv(p, sep=None, engine="python")
        except Exception:
            continue
        for col in df.columns:
            cl = str(col).strip().lower()
            if "silhouette" in cl or cl == "stability":
                vals = pd.to_numeric(df[col], errors="coerce").dropna()
                if 0 < len(vals) <= k and vals.between(-1.0, 1.0).all():
                    m = float(vals.min())
                    if best is None or m < best:
                        best = m
    return best


def report_real_min_silhouette(extractor_output, label, kind):
    log(f"Real per-signature minimum silhouette ({label}, {kind}):")
    sub_dir = Path(extractor_output) / KIND[kind]["extractor_subdir"] / "All_Solutions"
    if not sub_dir.exists():
        log("  (no All_Solutions/ directory)")
        return
    prefix = KIND[kind]["extractor_kdir_pre"]
    rows = []
    for k_dir in sorted(sub_dir.glob(f"{prefix}_*_Signatures"),
                        key=lambda d: _k_from_dirname(d) or 0):
        k = _k_from_dirname(k_dir)
        if k is None:
            continue
        rows.append((k, _scan_min_silhouette(k_dir, k)))
    if not rows:
        log("  (no per-k subdirectories found)")
        return
    log("  k   | min per-signature silhouette")
    log("  ----+-------------------------------")
    for k, m in rows:
        log(f"  {k:<3} | {m:.4f}" if m is not None else
            f"  {k:<3} | n/a (file not located)")


# --- Step 4: pan-cancer Extractor de novo ----------------------------------

def step4_pancancer_extractor(sbs_panc, id_panc):
    """Run de novo NMF on the pan-cancer PPP-HTG SBS96 + ID83 aggregates.

    Returns (sbs_cosmic_db, id_cosmic_db) -- paths to the COSMIC-decomposed
    signature files used by step 11 (constrained Assignment).
    """
    sbs_extractor_out = WORKDIR / "pancancer" / PANCANCER_LABEL / "extractor"
    id_extractor_out  = WORKDIR / "pancancer" / PANCANCER_INDEL_LABEL / "extractor"
    sbs_cosmic_db = run_extractor(sbs_panc, sbs_extractor_out, PANCANCER_LABEL, "SBS96")
    report_real_min_silhouette(sbs_extractor_out, PANCANCER_LABEL, "SBS96")
    id_cosmic_db  = run_extractor(id_panc,  id_extractor_out,  PANCANCER_INDEL_LABEL, "ID83")
    report_real_min_silhouette(id_extractor_out, PANCANCER_INDEL_LABEL, "ID83")
    return sbs_cosmic_db, id_cosmic_db


# --- Steps 5-10: pan-cancer Extractors on alternative aggregates ----------
# Compartment specificity (5-7):
#   step  5: pan-cancer PPP-LTG aggregate Extractor (PPP-HTG specificity control)
#   step  6: pan-cancer non-promoter SBS96 Extractor (PPP specificity control)
#   step  7: pan-cancer CFS aggregate Extractor (unified mechanism, all 17)
# Compositional bias (8-10):
#   step  8: liver-excluded pan-cancer Extractor
#   step  9: pancreas-excluded pan-cancer Extractor
#   step 10: equal-weight pan-cancer Extractor (cap=50, seed=42)
# All are resumable -- skip work whose Extractor cosmic-DB output exists.
# All use the KIND dispatch table for SBS96/ID83 parameterization.

def _run_sensitivity_extractor(matrix_path, out_dir, label, kind):
    """Wrapper that runs Extractor + reports min silhouette + reports
    TRC cosine. Skips silently if matrix_path is None (upstream
    aggregator returned no data)."""
    if matrix_path is None:
        return
    extractor_out = out_dir / f"extractor_{kind}"
    if cosmic_signatures_path(extractor_out, kind).exists():
        log(f"  [{kind}] exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, label, kind)
        report_trc_cosine(extractor_out, label, kind)
        return
    run_extractor(matrix_path, extractor_out, label, kind)
    report_real_min_silhouette(extractor_out, label, kind)
    report_trc_cosine(extractor_out, label, kind)


def step5_pancancer_low_extractor():
    log("Step 5: pan-cancer PPP-LTG Extractor (PPP-HTG specificity control, SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_low"
    sbs_matrix = aggregate_low("SBS96")
    id_matrix  = aggregate_low("ID83")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_promoter_low", "SBS96")
    _run_sensitivity_extractor(id_matrix,  out_dir, "pancancer_promoter_low", "ID83")


def step6_pancancer_nonpromoter_extractor():
    log("Step 6: pan-cancer non-promoter SBS96 Extractor (PPP specificity)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_nonpromoter"
    sbs_matrix = aggregate_nonpromoter("SBS96")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_non_promoter", "SBS96")


def step7_pancancer_cfs_extractor():
    log(f"Step 7: pan-cancer CFS Extractor (all {len(TUMORS)} tumors, SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_cfs"
    sbs_matrix = aggregate_cfs("SBS96")
    id_matrix  = aggregate_cfs("ID83")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_cfs", "SBS96")
    _run_sensitivity_extractor(id_matrix,  out_dir, "pancancer_cfs", "ID83")


def step8_liver_excluded():
    log("Step 8: liver-excluded pan-cancer Extractor (SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_no_liver"
    sbs_matrix = make_excluded_aggregate(
        "SBS96", LIVER_EXCLUDED_TUMORS, out_dir,
        "pancancer_promoter_high_no_liver",
    )
    id_matrix = make_excluded_aggregate(
        "ID83", LIVER_EXCLUDED_TUMORS, out_dir,
        "pancancer_promoter_high_indel_no_liver",
    )
    _run_sensitivity_extractor(sbs_matrix, out_dir,
        "pancancer_promoter_high_no_liver", "SBS96")
    _run_sensitivity_extractor(id_matrix, out_dir,
        "pancancer_promoter_high_indel_no_liver", "ID83")


def step9_pancreas_excluded():
    log("Step 9: pancreas-excluded pan-cancer Extractor (SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_no_pancreas"
    sbs_matrix = make_excluded_aggregate(
        "SBS96", PANCREAS_EXCLUDED_TUMORS, out_dir,
        "pancancer_promoter_high_no_pancreas",
    )
    id_matrix = make_excluded_aggregate(
        "ID83", PANCREAS_EXCLUDED_TUMORS, out_dir,
        "pancancer_promoter_high_indel_no_pancreas",
    )
    _run_sensitivity_extractor(sbs_matrix, out_dir,
        "pancancer_promoter_high_no_pancreas", "SBS96")
    _run_sensitivity_extractor(id_matrix, out_dir,
        "pancancer_promoter_high_indel_no_pancreas", "ID83")


def step10_equal_weight():
    log(f"Step 10: equal-weight pan-cancer Extractor "
        f"(cap={EQUAL_WEIGHT_CAP}, seed={EQUAL_WEIGHT_SEED}, SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_capped50"
    sbs_matrix = make_capped_aggregate(
        "SBS96", EQUAL_WEIGHT_CAP, EQUAL_WEIGHT_SEED, out_dir,
        f"pancancer_promoter_high_capped{EQUAL_WEIGHT_CAP}_seed{EQUAL_WEIGHT_SEED}",
    )
    id_matrix = make_capped_aggregate(
        "ID83", EQUAL_WEIGHT_CAP, EQUAL_WEIGHT_SEED, out_dir,
        f"pancancer_promoter_high_indel_capped{EQUAL_WEIGHT_CAP}_seed{EQUAL_WEIGHT_SEED}",
    )
    _run_sensitivity_extractor(sbs_matrix, out_dir,
        f"pancancer_capped{EQUAL_WEIGHT_CAP}", "SBS96")
    _run_sensitivity_extractor(id_matrix, out_dir,
        f"pancancer_capped{EQUAL_WEIGHT_CAP}_indel", "ID83")


# --- Step 11: constrained Assignment ----------------------------------------

def _assign(label, matrix_path, out_dir, cosmic_db, kind):
    if matrix_path is None or not Path(matrix_path).exists():
        return
    if (out_dir / "Assignment_Solution").exists():
        return
    if cosmic_db is None or not Path(cosmic_db).exists():
        log(f"  SKIP {label} [{kind}]: signature DB missing")
        return
    info = KIND[kind]
    log(f"  assign {label} [{kind}]")
    out_dir.mkdir(parents=True, exist_ok=True)
    from SigProfilerAssignment import Analyzer as Analyze
    try:
        Analyze.cosmic_fit(
            samples=str(matrix_path),
            output=str(out_dir),
            input_type="matrix",
            context_type=info["assign_context_type"],
            genome_build=REFERENCE,
            cosmic_version=COSMIC_VER,
            signature_database=str(cosmic_db),
            collapse_to_SBS96=info["collapse_to_SBS96"],
            nnls_add_penalty=0.05,
            nnls_remove_penalty=0.01,
        )
    except Exception as exc:
        log(f"    FAILED {label} [{kind}]: {exc}")


def assignment_out_dir(project, kind):
    suffix = "" if kind == "SBS96" else "_id"
    return WORKDIR / "assignment" / f"{project}{suffix}"


def step11_constrained_assignment(sbs_cosmic_db, id_cosmic_db):
    log("Step 11: constrained Assignment (SBS + ID83)")

    # Pan-cancer aggregates (constrained refit on the same input the
    # Extractor ran on -- methodological transparency).
    _assign(
        PANCANCER_LABEL, pancancer_matrix_path("SBS96"),
        assignment_out_dir(PANCANCER_LABEL, "SBS96"),
        sbs_cosmic_db, "SBS96",
    )
    _assign(
        PANCANCER_INDEL_LABEL, pancancer_matrix_path("ID83"),
        assignment_out_dir(PANCANCER_INDEL_LABEL, "ID83"),
        id_cosmic_db, "ID83",
    )

    # Per-(tumor, group) cohorts: SBS for every group, ID for every group.
    sbs_labels = sorted({l for _, _, l in GROUPS})
    id_labels  = sorted({l for _, _, l in INDEL_GROUPS})
    for tumor in TUMORS:
        for label in sbs_labels:
            project = f"{tumor}_{label}"
            _assign(
                project, matgen_matrix_path(project, "SBS96"),
                assignment_out_dir(project, "SBS96"),
                sbs_cosmic_db, "SBS96",
            )
        for label in id_labels:
            project = f"{tumor}_{label}"
            _assign(
                project, matgen_matrix_path(project, "ID83"),
                assignment_out_dir(project, "ID83"),
                id_cosmic_db, "ID83",
            )


# --- Step 12: unconstrained Assignment (sensitivity comparison) -------------

def _unconstrained_one(label, matrix_path, out_dir, kind):
    if matrix_path is None or not Path(matrix_path).exists():
        log(f"  ABORT [{kind}]: matrix missing for {label}")
        return
    if (out_dir / "Assignment_Solution").exists():
        log(f"  exists: {out_dir}")
        return
    info = KIND[kind]
    log(f"  assign {label} [{kind}] (full COSMIC database, signature_database=None)")
    out_dir.mkdir(parents=True, exist_ok=True)
    from SigProfilerAssignment import Analyzer as Analyze
    Analyze.cosmic_fit(
        samples=str(matrix_path),
        output=str(out_dir),
        input_type="matrix",
        context_type=info["assign_context_type"],
        genome_build=REFERENCE,
        cosmic_version=COSMIC_VER,
        signature_database=None,
        collapse_to_SBS96=info["collapse_to_SBS96"],
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
    )


def step12_unconstrained_assignment():
    log("Step 12: unconstrained Assignment on pan-cancer SBS96 + ID83")
    _unconstrained_one(
        PANCANCER_LABEL, pancancer_matrix_path("SBS96"),
        WORKDIR / "assignment" / f"{PANCANCER_LABEL}_unconstrained",
        "SBS96",
    )
    _unconstrained_one(
        PANCANCER_INDEL_LABEL, pancancer_matrix_path("ID83"),
        WORKDIR / "assignment" / f"{PANCANCER_INDEL_LABEL}_unconstrained",
        "ID83",
    )


# --- Step 13: no-hypermutator sensitivity (count-based) -------------------

def _no_hyper_one(panc_matrix, kind, label):
    if panc_matrix is None or not Path(panc_matrix).exists():
        log(f"  [{kind}] ABORT: pan-cancer matrix missing")
        return
    no_hyper_dir = WORKDIR / "pancancer" / f"{label}_no_hyper"
    no_hyper_matrix = no_hyper_dir / f"{label}_no_hyper.{KIND[kind]['matrix_ext']}"
    extractor_out = no_hyper_dir / "extractor"

    if cosmic_signatures_path(extractor_out, kind).exists():
        log(f"  [{kind}] exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, f"{label}_no_hyper", kind)
        report_trc_cosine(extractor_out, f"{label}_no_hyper", kind)
        return

    no_hyper_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(panc_matrix, sep="\t", index_col=0)
    counts = df.sum(axis=0)
    mu, sd = float(counts.mean()), float(counts.std())
    cutoff = mu + 2 * sd
    keep = counts[counts <= cutoff].index
    drop = sorted(set(counts.index) - set(keep))
    log(f"  [{kind}] cohort: n={len(counts)}, mean={mu:.0f}, sd={sd:.0f}, "
        f"cutoff={cutoff:.0f}")
    log(f"  [{kind}] dropping {len(drop)} hypermutator(s): {drop}")

    df[keep].to_csv(no_hyper_matrix, sep="\t")
    log(f"  [{kind}] filtered matrix: {no_hyper_matrix} ({len(keep)} samples)")

    run_extractor(no_hyper_matrix, extractor_out, f"{label}_no_hyper", kind)
    report_real_min_silhouette(extractor_out, f"{label}_no_hyper", kind)
    report_trc_cosine(extractor_out, f"{label}_no_hyper", kind)


def step13_no_hyper(panc_sbs_matrix, panc_id_matrix):
    log("Step 13: no-hypermutator sensitivity (count-based, SBS + ID)")
    _no_hyper_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL)
    _no_hyper_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL)


# --- Step 14 detail: ------------------------------------------------------
# Drop donors whose *unconstrained* SBS96 Assignment profile is dominated
# by MMR- or POLE-deficiency signatures, then rerun Extractor. PCAWG-
# consistent (mechanism-defined) alternative to step 13's count-based
# mean+2SD cutoff -- which conflates true MMR/POLE hypermutators with
# high-TRC donors. See docs/PLANNED_SENSITIVITY_STEPS.md section 6b
# for the full contract.

# Hypermutator-defining SBS signatures (COSMIC v3.5).
MMR_SIGS  = ["SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"]
POLE_SIGS = ["SBS10a", "SBS10b", "SBS10c", "SBS10d"]
HYPER_FRACTION_THRESHOLD = 0.30


def _identify_hypermutators_by_mechanism():
    """Donor IDs (e.g. 'pancreas__DO35442') whose unconstrained SBS96
    Assignment profile has MMR fraction >= HYPER_FRACTION_THRESHOLD OR
    POLE fraction >= HYPER_FRACTION_THRESHOLD. Returns None if the
    unconstrained Assignment activities file is missing.
    """
    activities_file = (
        WORKDIR / "assignment" / f"{PANCANCER_LABEL}_unconstrained"
        / "Assignment_Solution" / "Activities"
        / "Assignment_Solution_Activities.txt"
    )
    if not activities_file.exists():
        log(f"  ABORT: unconstrained Assignment not found: {activities_file}")
        return None
    act = pd.read_csv(activities_file, sep="\t", index_col=0)
    total = act.sum(axis=1).replace(0, pd.NA)
    mmr_cols  = [c for c in MMR_SIGS  if c in act.columns]
    pole_cols = [c for c in POLE_SIGS if c in act.columns]
    mmr_frac  = (act[mmr_cols].sum(axis=1)  / total).fillna(0.0)
    pole_frac = (act[pole_cols].sum(axis=1) / total).fillna(0.0)
    mmr_hits  = act.index[mmr_frac  >= HYPER_FRACTION_THRESHOLD]
    pole_hits = act.index[pole_frac >= HYPER_FRACTION_THRESHOLD]
    hyper = set(mmr_hits) | set(pole_hits)
    log(f"  unconstrained Assignment cohort: n={len(act)}")
    log(f"  MMR sigs present:  {mmr_cols}")
    log(f"  POLE sigs present: {pole_cols}")
    log(f"  MMR-dominant donors  (>= {HYPER_FRACTION_THRESHOLD:.0%}): "
        f"{len(mmr_hits)} -> {sorted(mmr_hits)}")
    log(f"  POLE-dominant donors (>= {HYPER_FRACTION_THRESHOLD:.0%}): "
        f"{len(pole_hits)} -> {sorted(pole_hits)}")
    log(f"  union to drop: {len(hyper)} donors")
    return hyper


def _no_mmr_pole_one(panc_matrix, kind, label, hyper_donors):
    if panc_matrix is None or not Path(panc_matrix).exists():
        log(f"  [{kind}] ABORT: pan-cancer matrix missing")
        return
    out_dir = WORKDIR / "sensitivity" / "pancancer_no_mmr_pole"
    maf_dir = out_dir / "maf"
    filt_matrix = maf_dir / f"{label}_no_mmr_pole.{KIND[kind]['matrix_ext']}"
    extractor_out = out_dir / f"extractor_{kind}"

    if cosmic_signatures_path(extractor_out, kind).exists():
        log(f"  [{kind}] exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, f"{label}_no_mmr_pole", kind)
        report_trc_cosine(extractor_out, f"{label}_no_mmr_pole", kind)
        return

    maf_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(panc_matrix, sep="\t", index_col=0)
    n_in = df.shape[1]
    keep = [c for c in df.columns if c not in hyper_donors]
    drop = sorted(set(df.columns) - set(keep))
    log(f"  [{kind}] cohort before: n={n_in}; "
        f"dropping {len(drop)} mechanism-defined hypermutator(s): {drop}")

    df[keep].to_csv(filt_matrix, sep="\t")
    log(f"  [{kind}] filtered matrix: {filt_matrix} ({len(keep)} samples)")

    run_extractor(filt_matrix, extractor_out, f"{label}_no_mmr_pole", kind)
    report_real_min_silhouette(extractor_out, f"{label}_no_mmr_pole", kind)
    report_trc_cosine(extractor_out, f"{label}_no_mmr_pole", kind)


def step14_no_mmr_pole(panc_sbs_matrix, panc_id_matrix):
    """Mechanism-based hypermutator removal sensitivity analysis.

    Defines hypermutators as donors with >= HYPER_FRACTION_THRESHOLD MMR
    (SBS6/14/15/20/21/26/44) OR >= HYPER_FRACTION_THRESHOLD POLE
    (SBS10a-d) attribution in the *unconstrained* pan-cancer SBS96
    Assignment. Applies the same donor exclusion list to both the SBS96
    and ID83 pan-cancer matrices and reruns Extractor.

    PCAWG-consistent (mechanism-defined) alternative to the count-based
    mean+2SD cutoff used in step13_no_hyper. Requires
    step12_unconstrained_assignment() to have produced
    `<PANCANCER_LABEL>_unconstrained` Assignment output beforehand.
    """
    log("Step 14: mechanism-based hypermutator removal (no-MMR/POLE, SBS + ID)")
    hyper = _identify_hypermutators_by_mechanism()
    if hyper is None:
        return
    _no_mmr_pole_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL, hyper)
    _no_mmr_pole_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL, hyper)


# --- Step 15: whole-genome SBS96 count hypermutator removal ---------------
# Drop donors whose unfiltered whole-genome SBS96 mutation count exceeds
# mean+2SD over the 658 PPP-HTG cohort. Single donor exclusion set is
# applied uniformly to BOTH the SBS96 and ID83 pan-cancer matrices --
# "hypermutator" is treated as a donor-level property and the SBS96
# whole-genome count is the more sensitive marker. Implementation
# decisions logged 2026-04-27 in docs/PLANNED_SENSITIVITY_STEPS.md §8.

def _load_whole_genome_sbs_counts():
    """pd.Series indexed by '<tumor>__<donor_id>', values = whole-genome
    `unfiltered_sbs_muts` from each per-tumor metadata file."""
    parts = []
    for tumor in TUMORS:
        meta = (ICGC_OUTPUTS / f"DCO__{TIMESTAMP_TAG}_{tumor}"
                / f"{tumor}_metadata.tsv")
        if not meta.exists():
            log(f"  [WGS-count] missing metadata for {tumor}: {meta}")
            continue
        m = pd.read_csv(meta, sep="\t")
        m = m.dropna(subset=["donor_id"])
        m["__key"] = m["donor_id"].apply(lambda d: f"{tumor}__{d}")
        parts.append(m.set_index("__key")["unfiltered_sbs_muts"])
    if not parts:
        return pd.Series(dtype=float)
    return pd.concat(parts)


def _identify_hypermutators_by_wgs_count(panc_sbs_matrix):
    """Donors whose whole-genome SBS96 count exceeds mean+2SD over the
    SBS96 pan-cancer matrix donors. Returns set or None on failure."""
    if panc_sbs_matrix is None or not Path(panc_sbs_matrix).exists():
        log(f"  ABORT: SBS96 pan-cancer matrix missing; "
            f"can't define WGS-hyper cohort")
        return None
    panc = pd.read_csv(panc_sbs_matrix, sep="\t", index_col=0, nrows=1)
    cohort = list(panc.columns)
    log(f"  cohort (from SBS96 pan-cancer matrix): n={len(cohort)}")
    counts_all = _load_whole_genome_sbs_counts()
    if counts_all.empty:
        log(f"  ABORT: no whole-genome SBS counts loaded from metadata files")
        return None
    counts = counts_all.reindex(cohort)
    missing = counts[counts.isna()].index.tolist()
    if missing:
        log(f"  WARNING: {len(missing)} cohort donors missing from metadata: "
            f"{missing}")
    counts = counts.dropna()
    mu, sd = float(counts.mean()), float(counts.std())
    cutoff = mu + 2 * sd
    hyper = set(counts[counts > cutoff].index)
    log(f"  WGS SBS96 cohort: n={len(counts)}, mean={mu:.0f}, sd={sd:.0f}, "
        f"cutoff={cutoff:.0f}")
    log(f"  WGS-hyper donors (> cutoff): {len(hyper)} -> {sorted(hyper)}")
    return hyper


def _no_wgs_hyper_one(panc_matrix, kind, label, hyper_donors):
    if panc_matrix is None or not Path(panc_matrix).exists():
        log(f"  [{kind}] ABORT: pan-cancer matrix missing")
        return
    out_dir = WORKDIR / "sensitivity" / "pancancer_no_wgs_hyper"
    maf_dir = out_dir / "maf"
    filt_matrix = maf_dir / f"{label}_no_wgs_hyper.{KIND[kind]['matrix_ext']}"
    extractor_out = out_dir / f"extractor_{kind}"

    if cosmic_signatures_path(extractor_out, kind).exists():
        log(f"  [{kind}] exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, f"{label}_no_wgs_hyper", kind)
        report_trc_cosine(extractor_out, f"{label}_no_wgs_hyper", kind)
        return

    maf_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(panc_matrix, sep="\t", index_col=0)
    n_in = df.shape[1]
    keep = [c for c in df.columns if c not in hyper_donors]
    drop = sorted(set(df.columns) - set(keep))
    log(f"  [{kind}] cohort before: n={n_in}; "
        f"dropping {len(drop)} WGS-hyper donor(s): {drop}")

    df[keep].to_csv(filt_matrix, sep="\t")
    log(f"  [{kind}] filtered matrix: {filt_matrix} ({len(keep)} samples)")

    run_extractor(filt_matrix, extractor_out, f"{label}_no_wgs_hyper", kind)
    report_real_min_silhouette(extractor_out, f"{label}_no_wgs_hyper", kind)
    report_trc_cosine(extractor_out, f"{label}_no_wgs_hyper", kind)


def step15_no_wgs_hyper(panc_sbs_matrix, panc_id_matrix):
    """Whole-genome SBS96 count hypermutator removal.

    Defines a single donor exclusion set from whole-genome unfiltered
    SBS96 counts (mean+2SD over the SBS96 pan-cancer cohort) and applies
    it uniformly to both the SBS96 and ID83 pan-cancer matrices, then
    reruns Extractor on each filtered matrix.
    """
    log("Step 15: whole-genome SBS96 count hypermutator removal (SBS + ID)")
    hyper = _identify_hypermutators_by_wgs_count(panc_sbs_matrix)
    if hyper is None:
        return
    _no_wgs_hyper_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL, hyper)
    _no_wgs_hyper_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL, hyper)


# --- Step 16: per-tumor SBS96 Extractor (PPP-HTG vs PPP-LTG validation) -------------

def step16_per_tumor_extractor():
    log("Step 16: per-tumor SBS96 Extractor (PPP-HTG vs PPP-LTG validation)")
    for tumor in TUMORS:
        for label in PER_TUMOR_GROUPS:
            project = f"{tumor}_{label}"
            matrix = matgen_matrix_path(project, "SBS96")
            if not matrix.exists():
                log(f"  skip {project}: no SBS96 matrix")
                continue
            out_dir = WORKDIR / "per_tumor" / project / "extractor"
            cosmic_db = cosmic_signatures_path(out_dir, "SBS96")
            if cosmic_db.exists():
                log(f"  exists: {project}")
                continue
            run_extractor(matrix, out_dir, project, "SBS96")
            report_real_min_silhouette(out_dir, project, "SBS96")


# --- Step 17: Ovary HRD exclusion ------------------------------------------
# In ovary, HRD/SBS3 dominates and may mask SBS96-TRC. Drop ovary donors
# with SBS3 fraction > HRD_SBS3_THRESHOLD, regenerate the ovary SBS96 matrix
# from a filtered MAF, and rerun Extractor.
#
# Source for SBS3 attribution: the per-tumor ovary Extractor's
# COSMIC-decomposed Activities file (produced by step 16). This file
# decomposes the ovary k=1 de novo signature into 5 COSMIC sigs
# (SBS1/3/5/98/102) and gives per-donor activities under that
# decomposition.
#
# Bug history (2026-04-27):
#   v1: read from CONSTRAINED Assignment -- failed silently because
#       constrained uses only the 10 COSMIC sigs decomposing the pan-cancer
#       Extractor signatures, a set that does NOT include SBS3. The
#       'SBS3 not in columns' branch always fired -> zero donors excluded.
#   v2: switched to UNCONSTRAINED pan-cancer Assignment -- also failed.
#       SigProfilerAssignment's NNLS-with-penalty (nnls_remove_penalty=0.01)
#       zeros out flat/diffuse signatures like SBS3 from per-donor refits
#       at low PPP-HTG mutation counts (median ~30 muts/ovary donor).
#       Verified empirically: 0/658 pan-cancer donors had SBS3 > 0 in the
#       unconstrained Assignment, so the threshold check excluded zero.
#   v3 (this version): read from the per-tumor ovary Extractor's
#       COSMIC-decomposed Activities. The decomposition is fit on the
#       cohort spectrum (~2,700 ovary mutations total) before the
#       per-donor projection, so SBS3 survives. Bimodal distribution:
#       40 donors at SBS3=0, 47 donors at SBS3>50%, only 2 between
#       30%-50%. Threshold is stable across 30-60% range.
# See docs/PLANNED_SENSITIVITY_STEPS.md §0 for the full incident log.

def step17_ovary_hrd_exclusion():
    log("Step 17: ovary HRD exclusion (remove SBS3-dominant donors, rerun Extractor)")
    ovary_project = f"ovary_{PRIMARY_GROUP}"

    # Per-tumor ovary Extractor's COSMIC-decomposed Activities (produced
    # by step 16). Donor IDs in this file are bare (e.g. 'DO35442') -- they
    # match the per-tumor matrix and the MAF Tumor_Sample_Barcode, so no
    # tumor-prefix stripping is needed for the downstream MAF filter.
    activities_file = (
        WORKDIR / "per_tumor" / ovary_project / "extractor" / "SBS96"
        / "Suggested_Solution" / "COSMIC_SBS96_Decomposed_Solution"
        / "Activities" / "COSMIC_SBS96_Activities.txt"
    )
    if not activities_file.exists():
        log(f"  ABORT: per-tumor ovary COSMIC-decomposed Activities not found "
            f"at {activities_file}. Has step 16 (per-tumor Extractor) run for ovary?")
        return

    out_root = WORKDIR / "ovary_hrd_excluded"
    excl_project = "ovary_promoter_high_no_hrd"
    excl_maf_dir = out_root / "maf" / excl_project
    excl_maf = excl_maf_dir / f"{excl_project}.maf"
    excl_extractor_out = out_root / "extractor"

    if cosmic_signatures_path(excl_extractor_out, "SBS96").exists():
        log(f"  exists: {excl_extractor_out}")
        report_real_min_silhouette(excl_extractor_out, excl_project, "SBS96")
        report_trc_cosine(excl_extractor_out, excl_project, "SBS96")
        return

    act = pd.read_csv(activities_file, sep="\t", index_col=0)
    if "SBS3" not in act.columns:
        # Per-tumor ovary k=1 decomposes to {SBS1, SBS3, SBS5, SBS98, SBS102}
        # by construction in the canonical run. If SBS3 is absent, the
        # canonical per-tumor extractor was rerun with a different result
        # and the HRD exclusion premise needs revisiting.
        log(f"  ABORT: SBS3 not in per-tumor decomposed Activities columns: "
            f"{list(act.columns)}. Has the canonical ovary extractor been "
            f"rerun with a different decomposition?")
        return
    total = act.sum(axis=1).replace(0, pd.NA)
    sbs3_frac = (act["SBS3"] / total).fillna(0.0)
    hrd_donors = set(sbs3_frac[sbs3_frac > HRD_SBS3_THRESHOLD].index)

    n_total = len(act)
    n_hrd   = len(hrd_donors)
    pct     = (100.0 * n_hrd / n_total) if n_total else 0.0
    log(f"  ovary cohort: {n_total} donors, {n_hrd} HRD ({pct:.1f}%) excluded")
    log(f"  HRD donors: {sorted(hrd_donors)}")

    # Filter the ovary SBS MAF.
    src_maf = sbs_maf_path(ovary_project)
    if not src_maf.exists():
        log(f"  ABORT: source ovary SBS MAF missing: {src_maf}")
        return
    excl_maf_dir.mkdir(parents=True, exist_ok=True)
    n_in = n_kept = 0
    with open(src_maf) as fin, open(excl_maf, "w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        header = next(reader)
        writer.writerow(header)
        sample_ix = header.index("Tumor_Sample_Barcode")
        for row in reader:
            if not row:
                continue
            n_in += 1
            if row[sample_ix] not in hrd_donors:
                writer.writerow(row)
                n_kept += 1
    log(f"  filtered MAF: kept {n_kept}/{n_in} mutation rows -> {excl_maf}")

    # Rerun MatrixGenerator on the filtered MAF.
    sbs96 = excl_maf_dir / "output" / "SBS" / f"{excl_project}.SBS96.all"
    if not sbs96.exists():
        from SigProfilerMatrixGenerator.scripts import (
            SigProfilerMatrixGeneratorFunc as matGen,
        )
        log(f"  matgen on {excl_project}")
        try:
            matGen.SigProfilerMatrixGeneratorFunc(
                excl_project, REFERENCE, str(excl_maf_dir),
                plot=False, exome=False, bed_file=None,
                chrom_based=False, tsb_stat=False, seqInfo=False,
            )
        except Exception as exc:
            log(f"    matgen FAILED: {exc}")
            return
    if not sbs96.exists():
        log(f"  ABORT: matgen did not produce {sbs96}")
        return

    # Rerun Extractor on the HRD-excluded ovary SBS96 matrix.
    run_extractor(sbs96, excl_extractor_out, excl_project, "SBS96")
    report_real_min_silhouette(excl_extractor_out, excl_project, "SBS96")
    report_trc_cosine(excl_extractor_out, excl_project, "SBS96")


# --- Step 18: Pan-cancer rerun with ovary HRD donors excluded --------------
# Step 17 reruns the Extractor on the *ovary-alone* HRD-excluded matrix
# (n=42), which gives a per-tumor sensitivity check on what the ovary
# signature looks like after HRD removal. Step 18 is different: it drops
# the same SBS3>50% ovary donors from the *pan-cancer* aggregate and
# reruns the pan-cancer SBS96 + ID83 Extractors. This is the
# manuscript-relevant test:
#   1. Robustness of pan-cancer SBS96-TRC -- cosine to the original
#      658-donor SBS96-TRC after dropping 47 HRD donors. Expected high
#      (>= 0.95) if SBS96-TRC is independent of HRD contamination.
#   2. Drops SBS3 from the COSMIC pairwise novelty list -- currently
#      SBS3 is rank-2 closest match (cosine 0.7027) because the HRD
#      donors pull the de novo signature toward SBS3's spectrum. Without
#      them, the rank-2 hit should change and the "no COSMIC sig
#      exceeds 0.80" novelty claim becomes cleaner.
#
# Reuses the same SBS3 source and threshold as step 17 -- reads donors
# from per-tumor ovary's COSMIC-decomposed Activities (step 16 output).

def step18_ovary_hrd_excluded():
    log("Step 18: ovary-HRD-excluded pan-cancer Extractor (SBS96 + ID83)")
    ovary_project = f"ovary_{PRIMARY_GROUP}"

    activities_file = (
        WORKDIR / "per_tumor" / ovary_project / "extractor" / "SBS96"
        / "Suggested_Solution" / "COSMIC_SBS96_Decomposed_Solution"
        / "Activities" / "COSMIC_SBS96_Activities.txt"
    )
    if not activities_file.exists():
        log(f"  ABORT: per-tumor ovary COSMIC-decomposed Activities not found "
            f"at {activities_file}. Has step 16 (per-tumor Extractor) run for ovary?")
        return

    act = pd.read_csv(activities_file, sep="\t", index_col=0)
    if "SBS3" not in act.columns:
        log(f"  ABORT: SBS3 not in per-tumor decomposed Activities columns: "
            f"{list(act.columns)}.")
        return
    total = act.sum(axis=1).replace(0, pd.NA)
    sbs3_frac = (act["SBS3"] / total).fillna(0.0)
    hrd_donor_ids = sorted(sbs3_frac[sbs3_frac > HRD_SBS3_THRESHOLD].index)
    hrd_columns = {f"ovary__{d}" for d in hrd_donor_ids}
    log(f"  ovary HRD donors to exclude (SBS3>{HRD_SBS3_THRESHOLD:.0%}): "
        f"{len(hrd_donor_ids)}")

    out_dir = WORKDIR / "sensitivity" / "pancancer_no_ovary_hrd"
    sbs_matrix = make_donor_excluded_aggregate(
        "SBS96", hrd_columns, out_dir,
        "pancancer_promoter_high_no_ovary_hrd",
    )
    id_matrix = make_donor_excluded_aggregate(
        "ID83", hrd_columns, out_dir,
        "pancancer_promoter_high_indel_no_ovary_hrd",
    )
    _run_sensitivity_extractor(sbs_matrix, out_dir,
        "pancancer_promoter_high_no_ovary_hrd", "SBS96")
    _run_sensitivity_extractor(id_matrix, out_dir,
        "pancancer_promoter_high_indel_no_ovary_hrd", "ID83")


# --- Step 19: Within-cohort non-promoter SBS96 Extractor -------------------
# Step 6 ran a non-promoter SBS96 Extractor on 1,785 samples — a different
# cohort from PPP-HTG (658 donors): includes 4 no-expression tumors
# (esophagus, gallbladder, kidney, PNET) and mutation-only donors that
# the all-muts switch made visible. So step 6's negative-control reading
# ("PPP-specificity holds: non-promoter does not yield SBS96-TRC") is
# *informative* but not a perfectly-matched cohort comparison.
#
# Step 19 fixes that: it subsets the non-promoter SBS96 matrix to ONLY
# the donor columns that also appear in the canonical PPP-HTG aggregate
# (~658 donors). Runs the Extractor on that subset. Same compartment
# difference, same donor pool — isolates the *compartment* effect from
# the *cohort* effect.
#
# Pass criterion: same as step 6 — closest cosine to canonical SBS96-TRC
# < 0.85. If a TRC-like component re-emerges in the matched cohort,
# the negative-control argument weakens.

def step19_nonpromoter_htg_donors_extractor():
    log("Step 19: within-cohort non-promoter SBS96 Extractor "
        "(non-promoter mutations restricted to PPP-HTG donor cohort)")
    out_dir   = WORKDIR / "sensitivity" / "pancancer_nonpromoter_htg_donors"
    out_label = "pancancer_non_promoter_htg_donors"
    out_path  = out_dir / "maf" / f"{out_label}.{KIND['SBS96']['matrix_ext']}"
    extractor_out = out_dir / "extractor_SBS96"

    if cosmic_signatures_path(extractor_out, "SBS96").exists():
        log(f"  exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, out_label, "SBS96")
        report_trc_cosine(extractor_out, out_label, "SBS96")
        return

    htg_src     = pancancer_matrix_path("SBS96")
    nonprom_src = (WORKDIR / "sensitivity" / "pancancer_nonpromoter" / "maf"
                   / f"pancancer_non_promoter.{KIND['SBS96']['matrix_ext']}")
    if not htg_src.exists():
        log(f"  ABORT: canonical PPP-HTG matrix missing: {htg_src}")
        return
    if not nonprom_src.exists():
        log(f"  ABORT: non-promoter aggregate matrix missing: {nonprom_src}. "
            f"Has step 6 (or its aggregator) run?")
        return

    if not out_path.exists():
        out_path.parent.mkdir(parents=True, exist_ok=True)
        htg_cols = set(pd.read_csv(htg_src, sep="\t", index_col=0, nrows=0).columns)
        nonprom  = pd.read_csv(nonprom_src, sep="\t", index_col=0)
        n_in     = nonprom.shape[1]
        keep     = [c for c in nonprom.columns if c in htg_cols]
        log(f"  intersect: {len(keep)}/{n_in} non-promoter samples are in PPP-HTG cohort")
        if len(keep) == 0:
            log(f"  ABORT: zero overlap between non-promoter and PPP-HTG cohorts")
            return
        nonprom[keep].to_csv(out_path, sep="\t")
        log(f"  wrote {out_path} ({len(keep)} samples)")

    _run_sensitivity_extractor(out_path, out_dir, out_label, "SBS96")


# --- Step 20: Within-cohort CFS SBS96 Extractor ----------------------------
# Step 7 ran the CFS Extractor on all 17 tumors (1,694 samples) including
# esophagus, gallbladder, kidney, PNET (no-expression tumors) and
# mutation-only donors visible after the all-muts switch. Partial result
# at k=4 (2026-04-29): the closest CFS component is SBS40a-dominated
# (cos 0.85 to SBS40a alone, cos 0.83 to canonical SBS96-TRC) -- the
# higher single-COSMIC cos means the component is best explained as
# SBS40a, not TRC re-emerging. Other k=4 components are clock/dMMR/SBS17b
# (esophageal-driven).
#
# Step 20 fixes the cohort confound: subset the CFS SBS96 matrix to ONLY
# the 658 PPP-HTG donor columns (intersection = 657; one breast donor
# has no CFS row), then rerun the Extractor on that within-cohort matrix.
# Same compartment difference as step 7, same donor pool as step 4 (PPP-
# HTG canonical) -- this is the cleanest unified-mechanism test at the
# spectrum level. Compares "spectrum of CFS mutations in PPP-HTG donors"
# to "spectrum of PPP-HTG mutations in PPP-HTG donors" (canonical SBS96-
# TRC).
#
# Pass criterion: closest cosine of any de novo component to canonical
# SBS96-TRC > 0.85 AND that cosine exceeds the closest single-COSMIC
# match for the same component (so it's not just SBS40a-driven). If
# both hold, the unified-mechanism claim is supported at spectrum level.

def step20_cfs_htg_donors_extractor():
    log("Step 20: within-cohort CFS SBS96 Extractor "
        "(CFS mutations restricted to PPP-HTG donor cohort)")
    out_dir   = WORKDIR / "sensitivity" / "pancancer_cfs_htg_donors"
    out_label = "pancancer_cfs_htg_donors"
    out_path  = out_dir / "maf" / f"{out_label}.{KIND['SBS96']['matrix_ext']}"
    extractor_out = out_dir / "extractor_SBS96"

    if cosmic_signatures_path(extractor_out, "SBS96").exists():
        log(f"  exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, out_label, "SBS96")
        report_trc_cosine(extractor_out, out_label, "SBS96")
        return

    htg_src = pancancer_matrix_path("SBS96")
    cfs_src = (WORKDIR / "sensitivity" / "pancancer_cfs" / "maf"
               / f"pancancer_cfs.{KIND['SBS96']['matrix_ext']}")
    if not htg_src.exists():
        log(f"  ABORT: canonical PPP-HTG matrix missing: {htg_src}")
        return
    if not cfs_src.exists():
        log(f"  ABORT: CFS aggregate matrix missing: {cfs_src}. "
            f"Has step 7 (or its aggregator) run?")
        return

    if not out_path.exists():
        out_path.parent.mkdir(parents=True, exist_ok=True)
        htg_cols = set(pd.read_csv(htg_src, sep="\t", index_col=0, nrows=0).columns)
        cfs      = pd.read_csv(cfs_src, sep="\t", index_col=0)
        n_in     = cfs.shape[1]
        keep     = [c for c in cfs.columns if c in htg_cols]
        log(f"  intersect: {len(keep)}/{n_in} CFS samples are in PPP-HTG cohort")
        if len(keep) == 0:
            log(f"  ABORT: zero overlap between CFS and PPP-HTG cohorts")
            return
        cfs[keep].to_csv(out_path, sep="\t")
        log(f"  wrote {out_path} ({len(keep)} samples)")

    _run_sensitivity_extractor(out_path, out_dir, out_label, "SBS96")


# --- Step 21: Pan-cancer non-promoter ID83 Extractor ----------------------
# Step 6 ran the non-promoter Extractor on SBS96 only — when the
# PPP-specificity control was pre-registered (legacy label 3c), it was
# scoped as a SBS-only test. This left a gap in Fig 2's compartment-
# specificity matrix: the non-promoter row has SBS-TRC and SBS-BG cells
# but the ID-TRC and ID-BG cells are blank.
#
# Step 21 closes that gap. Per-tumor non-promoter ID83 matrices already
# exist (step 2 MatrixGenerator produced them) — we just need to
# aggregate them into a pan-cancer non-promoter ID83 matrix and run an
# Extractor on it. Same pattern as step 7 (CFS Extractor) but on
# non-promoter regions.
#
# Pass criterion: same as step 6 SBS96 — the PPP-specificity control
# expects NO extracted ID83 component to reach cosine ≥ 0.85 to
# canonical pan-cancer ID83-TRC. If a TRC-like ID component DOES
# emerge, that's a finding (would suggest indel TRC biology operates
# more broadly than substitution TRC biology).

def step21_pancancer_nonpromoter_id_extractor():
    log("Step 21: pan-cancer non-promoter ID83 Extractor (PPP specificity, ID)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_nonpromoter"
    id_matrix = aggregate_nonpromoter("ID83")
    _run_sensitivity_extractor(id_matrix, out_dir,
                                "pancancer_non_promoter", "ID83")


# --- Step 22: Zou KO validation (constrained Assignment with TRC ref) -----
# Tests whether SBS96-TRC -- defined from ICGC PPP-HTG mutations -- is
# recovered in subclones knocked out for TRC-resolution helicases (RECQL5,
# SETX, ATP2B4) from the Zou et al. 2021 (Nik-Zainal) iPSC KO catalogue.
# Comparator panel: WRN/PIF1 (other helicases -- specificity control),
# TP53/ATM/NHEJ1/XRCC4/REV1 (low-burden non-TRC -- negative control),
# MSH6 (MMR positive sanity, expected to load on SBS6/14/etc, not TRC).
#
# Per-subclone constrained Assignment against COSMIC v3.5 + SBS96-TRC +
# SBS96B (the two pan-cancer step-4 de novo signatures appended to
# COSMIC). The standalone runner is `validation/run_validation.py`; this
# step thin-wraps it so the analysis is part of the pipeline graph.

def step22_validation_kos():
    log("Step 22: Zou KO validation (constrained Assignment with SBS96-TRC ref)")
    out_dir = WORKDIR / "validation_kos"
    final_tsv = out_dir / "validation_ko_activities.tsv"
    if final_tsv.exists():
        log(f"  exists: {final_tsv}")
        return
    import subprocess
    runner = REPO_ROOT / "validation" / "run_validation.py"
    if not runner.exists():
        log(f"  SKIP: runner missing at {runner}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"  invoking {runner}")
    subprocess.run(
        ["python3", str(runner)],
        check=True,
        cwd=str(REPO_ROOT),
    )


# --- main -------------------------------------------------------------------

def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    log(f"WORKDIR = {WORKDIR}")

    # --- Setup (steps 1-3) -----------------------------------------------
    step1_convert_all()
    step2_matrix_generator()
    sbs_panc, id_panc = step3_aggregate_high()
    _sanity_check_after_step3()

    # --- Pan-cancer Extractors -- compartment specificity (steps 4-7) ---
    sbs_cosmic_db, id_cosmic_db = step4_pancancer_extractor(sbs_panc, id_panc)
    step5_pancancer_low_extractor()          # PPP-LTG (PPP-HTG specificity)
    step6_pancancer_nonpromoter_extractor()  # non-promoter (PPP specificity)
    step7_pancancer_cfs_extractor()          # CFS (unified mechanism, all 17)

    # --- Pan-cancer Extractors -- compositional bias (steps 8-10) -------
    step8_liver_excluded()                   # drop liver
    step9_pancreas_excluded()                # drop pancreas
    step10_equal_weight()                    # cap=50, seed=42

    # --- Assignment (steps 11-12) ---------------------------------------
    step11_constrained_assignment(sbs_cosmic_db, id_cosmic_db)
    _sanity_check_after_step11()
    step12_unconstrained_assignment()        # required by step 14

    # --- Pan-cancer Extractors -- hypermutator subsets (steps 13-15) ----
    step13_no_hyper(sbs_panc, id_panc)       # count-based (mean+2SD)
    step14_no_mmr_pole(sbs_panc, id_panc)    # mechanism-based (needs step 12)
    step15_no_wgs_hyper(sbs_panc, id_panc)   # whole-genome SBS96 count-based

    # --- Per-tumor + subset (steps 16-20) -------------------------------
    step16_per_tumor_extractor()             # PPP-HTG vs PPP-LTG validation
    step17_ovary_hrd_exclusion()             # ovary-alone, uses step 16 Activities
    step18_ovary_hrd_excluded()              # pan-cancer rerun, drops 47 HRD donors
    step19_nonpromoter_htg_donors_extractor()  # within-cohort non-promoter (step 6 cleanup)
    step20_cfs_htg_donors_extractor()        # within-cohort CFS (step 7 cleanup)
    step21_pancancer_nonpromoter_id_extractor()  # non-promoter ID83 (Fig 2 gap fill)

    # --- KO validation (step 22) ----------------------------------------
    step22_validation_kos()                       # Zou 2021 KO subclones

    log("DONE")


if __name__ == "__main__":
    main()
