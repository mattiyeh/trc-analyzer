#!/usr/bin/env python3
"""End-to-end SigProfiler workflow for the TRC pan-cancer study.

Pipeline:
  1.  Convert each ICGC TSV (SBS + indel) -> minimal MAF.
  2.  Run SigProfilerMatrixGenerator on every (tumor x group) -> SBS96
      and ID83 matrices.
  3.  Aggregate the 17 tumor-type "promoter_high" SBS96 + ID83 matrices
      into pan-cancer cohorts.
  4.  Run SigProfilerExtractor de novo on the pan-cancer SBS96 and
      pan-cancer ID83 aggregates (separate runs).
  4b. Run SigProfilerExtractor independently on each tumor's promoter_high
      and promoter_low SBS96 matrices (per-tumor HTG vs LTG validation).
  5.  Run constrained SigProfilerAssignment (SBS + ID83) on every
      (tumor x group) using the Extractor-derived COSMIC signatures.
  7.  Unconstrained Assignment (full COSMIC database) on the pan-cancer
      SBS96 and ID83 aggregates -- the constrained-vs-unconstrained
      sensitivity comparison.
  6.  No-hypermutator sensitivity analysis: drop samples > mean + 2 SD
      and rerun Extractor (both kinds).

Sensitivity / control steps (added per docs/PLANNED_SENSITIVITY_STEPS.md):
  3b. Pan-cancer PPP-LTG aggregate Extractor (HTG specificity control).
  3c. Pan-cancer non-promoter SBS96 Extractor (PPP specificity control).
  3d. Liver-excluded pan-cancer Extractor (compositional bias).
  3e. Equal-weight pan-cancer Extractor, cap 50 donors per tumor seed=42.
  3f. Pancreas-excluded pan-cancer Extractor (compositional bias).
  3g. Pan-cancer CFS aggregate Extractor (unified mechanism, all 17 tumors).
  6b. Mechanism-based hypermutator removal: drop donors with MMR >=30%
      OR POLE >=30% attribution in the unconstrained SBS96 Assignment
      and rerun Extractor (PCAWG-consistent alternative to step 6).
  6c. Whole-genome SBS96 count hypermutator removal: drop donors whose
      whole-genome unfiltered SBS96 count > mean+2SD over the 658 PPP-HTG
      cohort. Single donor exclusion set applied to both SBS96 and ID83.
  Ovary HRD exclusion (FIXED): remove ovary donors with SBS3 attribution
      > 50% of total SBS in the *unconstrained* ovary Assignment, rerun
      MatrixGenerator + Extractor on the filtered ovary cohort.

Resumable: each step skips work whose outputs already exist.
"""

import csv
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


# --- Configuration ----------------------------------------------------------

REPO_ROOT     = Path("/data/research/projects/trc_signatures")
ICGC_OUTPUTS  = REPO_ROOT / "outputs/java/2026.03.15"
TIMESTAMP_TAG = "2026.03.15_00.10.17"
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

# Equal-weight (step3e) compositional rebalancing parameters.
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


# --- Step 4: Extractor de novo ---------------------------------------------

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


# --- Step 4b: per-tumor SBS96 extraction -----------------------------------

def step4b_per_tumor_extractor():
    log("Step 4b: per-tumor SBS96 Extractor (HTG vs LTG validation)")
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


# --- Step 5: constrained Assignment ----------------------------------------

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


def step5_assignment(sbs_cosmic_db, id_cosmic_db):
    log("Step 5: constrained Assignment (SBS + ID83)")

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


# --- Step 7: unconstrained Assignment (sensitivity comparison) -------------

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


def step7_unconstrained_comparison():
    log("Step 7: unconstrained Assignment on pan-cancer SBS96 + ID83")
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


# --- Step 6: no-hypermutator sensitivity analysis --------------------------

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


def step6_no_hyper(panc_sbs_matrix, panc_id_matrix):
    log("Step 6: no-hypermutator sensitivity analysis (SBS + ID)")
    _no_hyper_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL)
    _no_hyper_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL)


# --- Step 6b: mechanism-based hypermutator removal -------------------------
# Drop donors whose *unconstrained* SBS96 Assignment profile is dominated
# by MMR- or POLE-deficiency signatures, then rerun Extractor. PCAWG-
# consistent (mechanism-defined) alternative to step6_no_hyper's
# count-based mean+2SD cutoff -- which conflates true MMR/POLE
# hypermutators with high-TRC donors. See
# docs/PLANNED_SENSITIVITY_STEPS.md section 6b for the full contract.

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


def step6b_no_mmr_pole(panc_sbs_matrix, panc_id_matrix):
    """Mechanism-based hypermutator removal sensitivity analysis.

    Defines hypermutators as donors with >= HYPER_FRACTION_THRESHOLD MMR
    (SBS6/14/15/20/21/26/44) OR >= HYPER_FRACTION_THRESHOLD POLE
    (SBS10a-d) attribution in the *unconstrained* pan-cancer SBS96
    Assignment. Applies the same donor exclusion list to both the SBS96
    and ID83 pan-cancer matrices and reruns Extractor.

    PCAWG-consistent (mechanism-defined) alternative to the count-based
    mean+2SD cutoff used in step6_no_hyper. Requires
    step7_unconstrained_comparison() to have produced
    `<PANCANCER_LABEL>_unconstrained` Assignment output beforehand.
    """
    log("Step 6b: mechanism-based hypermutator removal (SBS + ID)")
    hyper = _identify_hypermutators_by_mechanism()
    if hyper is None:
        return
    _no_mmr_pole_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL, hyper)
    _no_mmr_pole_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL, hyper)


# --- Step 3 sensitivity family --------------------------------------------
# step3b: pan-cancer PPP-LTG aggregate Extractor (HTG specificity control)
# step3c: pan-cancer non-promoter SBS96 Extractor (PPP specificity control)
# step3d: liver-excluded pan-cancer Extractor (compositional bias)
# step3e: equal-weight pan-cancer Extractor (compositional bias, cap=50)
# step3f: pancreas-excluded pan-cancer Extractor (compositional bias)
# step3g: pan-cancer CFS aggregate Extractor (unified mechanism, all 17 tumors)
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


def step3b_aggregate_low():
    log("Step 3b: pan-cancer PPP-LTG aggregate Extractor (SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_low"
    sbs_matrix = aggregate_low("SBS96")
    id_matrix  = aggregate_low("ID83")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_promoter_low", "SBS96")
    _run_sensitivity_extractor(id_matrix,  out_dir, "pancancer_promoter_low", "ID83")


def step3c_pancancer_nonpromoter():
    log("Step 3c: pan-cancer non-promoter SBS96 Extractor (PPP specificity)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_nonpromoter"
    sbs_matrix = aggregate_nonpromoter("SBS96")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_non_promoter", "SBS96")


def step3d_liver_excluded():
    log("Step 3d: liver-excluded pan-cancer Extractor (SBS96 + ID83)")
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


def step3e_equal_weight():
    log(f"Step 3e: equal-weight pan-cancer Extractor "
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


def step3f_pancreas_excluded():
    log("Step 3f: pancreas-excluded pan-cancer Extractor (SBS96 + ID83)")
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


def step3g_pancancer_cfs_extractor():
    log(f"Step 3g: pan-cancer CFS Extractor (all {len(TUMORS)} tumors, SBS96 + ID83)")
    out_dir = WORKDIR / "sensitivity" / "pancancer_cfs"
    sbs_matrix = aggregate_cfs("SBS96")
    id_matrix  = aggregate_cfs("ID83")
    _run_sensitivity_extractor(sbs_matrix, out_dir, "pancancer_cfs", "SBS96")
    _run_sensitivity_extractor(id_matrix,  out_dir, "pancancer_cfs", "ID83")


# --- Step 6c: whole-genome SBS96 count hypermutator removal ---------------
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


def step6c_no_wgs_hyper(panc_sbs_matrix, panc_id_matrix):
    """Whole-genome SBS96 count hypermutator removal.

    Defines a single donor exclusion set from whole-genome unfiltered
    SBS96 counts (mean+2SD over the SBS96 pan-cancer cohort) and applies
    it uniformly to both the SBS96 and ID83 pan-cancer matrices, then
    reruns Extractor on each filtered matrix.
    """
    log("Step 6c: whole-genome SBS96 count hypermutator removal (SBS + ID)")
    hyper = _identify_hypermutators_by_wgs_count(panc_sbs_matrix)
    if hyper is None:
        return
    _no_wgs_hyper_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL, hyper)
    _no_wgs_hyper_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL, hyper)


# --- Ovary HRD exclusion ---------------------------------------------------
# In ovary, HRD/SBS3 dominates and may mask SBS96-TRC. Drop ovary donors
# with SBS3 fraction > HRD_SBS3_THRESHOLD in the *unconstrained* ovary
# Assignment, regenerate the ovary SBS96 matrix from a filtered MAF, and
# rerun Extractor.
#
# Bug history: previously this function read SBS3 from the *constrained*
# Assignment, which uses only the 10 COSMIC sigs that decompose to the
# pan-cancer Extractor-derived signatures -- a set that does NOT include
# SBS3. The "SBS3 not in columns" branch always fired and zero donors
# were excluded. Fixed 2026-04-27 to read from the unconstrained
# Assignment (full COSMIC v3.5, includes SBS3). See
# docs/PLANNED_SENSITIVITY_STEPS.md §0.

def step_ovary_hrd_exclusion():
    log("Ovary HRD exclusion: identify and remove SBS3-dominant donors")
    ovary_project = f"ovary_{PRIMARY_GROUP}"

    # Bootstrap the unconstrained ovary Assignment if it doesn't exist
    # (step7 only covers the pan-cancer aggregate by default).
    ovary_unconstrained_dir = (
        WORKDIR / "assignment" / f"{ovary_project}_unconstrained"
    )
    if not (ovary_unconstrained_dir / "Assignment_Solution").exists():
        log(f"  ovary unconstrained Assignment missing -- running it now")
        _unconstrained_one(
            ovary_project, matgen_matrix_path(ovary_project, "SBS96"),
            ovary_unconstrained_dir, "SBS96",
        )

    activities_file = (
        ovary_unconstrained_dir
        / "Assignment_Solution" / "Activities"
        / "Assignment_Solution_Activities.txt"
    )
    if not activities_file.exists():
        log(f"  ABORT: ovary unconstrained Assignment activities not found at "
            f"{activities_file}")
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
    if "SBS3" in act.columns:
        total = act.sum(axis=1).replace(0, pd.NA)
        sbs3_frac = (act["SBS3"] / total).fillna(0.0)
        hrd_donors = set(sbs3_frac[sbs3_frac > HRD_SBS3_THRESHOLD].index)
    else:
        # With unconstrained Assignment + COSMIC v3.5, SBS3 should ALWAYS
        # be one of the 97 columns. Reaching this branch indicates a
        # broken upstream Assignment or a COSMIC version change.
        log(f"  WARNING: SBS3 not in unconstrained columns: "
            f"{list(act.columns)}")
        hrd_donors = set()

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


# --- main -------------------------------------------------------------------

def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    log(f"WORKDIR = {WORKDIR}")

    # --- Canonical pipeline ----------------------------------------------
    step1_convert_all()
    step2_matrix_generator()
    sbs_panc, id_panc = step3_aggregate_high()

    sbs_extractor_out = WORKDIR / "pancancer" / PANCANCER_LABEL / "extractor"
    id_extractor_out  = WORKDIR / "pancancer" / PANCANCER_INDEL_LABEL / "extractor"
    sbs_cosmic_db = run_extractor(sbs_panc, sbs_extractor_out, PANCANCER_LABEL,       "SBS96")
    report_real_min_silhouette(sbs_extractor_out, PANCANCER_LABEL,       "SBS96")
    id_cosmic_db  = run_extractor(id_panc,  id_extractor_out,  PANCANCER_INDEL_LABEL, "ID83")
    report_real_min_silhouette(id_extractor_out,  PANCANCER_INDEL_LABEL, "ID83")

    step4b_per_tumor_extractor()
    step5_assignment(sbs_cosmic_db, id_cosmic_db)
    step7_unconstrained_comparison()       # produces pan-cancer unconstrained
                                           # activities used by step6b
    step6_no_hyper(sbs_panc, id_panc)      # PPP-HTG count-based hyper removal

    # --- Sensitivity / control steps (docs/PLANNED_SENSITIVITY_STEPS.md) -
    step3b_aggregate_low()                 # HTG specificity (pan-cancer LTG)
    step3c_pancancer_nonpromoter()         # PPP specificity (non-promoter)
    step3d_liver_excluded()                # compositional bias: drop liver
    step3e_equal_weight()                  # compositional bias: cap=50, seed=42
    step3f_pancreas_excluded()             # compositional bias: drop pancreas
    step3g_pancancer_cfs_extractor()       # unified mechanism (CFS, all 17)
    step6b_no_mmr_pole(sbs_panc, id_panc)  # mechanism-based hyper removal
    step6c_no_wgs_hyper(sbs_panc, id_panc) # whole-genome SBS-count hyper removal
    step_ovary_hrd_exclusion()             # FIXED: reads SBS3 from
                                           # ovary unconstrained Assignment

    log("DONE")


if __name__ == "__main__":
    main()
