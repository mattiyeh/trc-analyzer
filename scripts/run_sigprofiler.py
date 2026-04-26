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
  Ovary HRD exclusion: remove ovary donors with SBS3 attribution > 50%
      of total SBS, rerun MatrixGenerator + Extractor on the filtered
      ovary cohort.

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


def step6_no_hyper(panc_sbs_matrix, panc_id_matrix):
    log("Step 6: no-hypermutator sensitivity analysis (SBS + ID)")
    _no_hyper_one(panc_sbs_matrix, "SBS96", PANCANCER_LABEL)
    _no_hyper_one(panc_id_matrix,  "ID83",  PANCANCER_INDEL_LABEL)


# --- Ovary HRD exclusion ---------------------------------------------------
# In ovary, HRD/SBS3 dominates and may mask SBS96B. After per-tumor Assignment,
# drop ovary donors with SBS3 fraction > HRD_SBS3_THRESHOLD, regenerate the
# ovary SBS96 matrix from a filtered MAF, and rerun Extractor.

def step_ovary_hrd_exclusion():
    log("Ovary HRD exclusion: identify and remove SBS3-dominant donors")
    ovary_project = f"ovary_{PRIMARY_GROUP}"

    activities_file = (
        assignment_out_dir(ovary_project, "SBS96")
        / "Assignment_Solution" / "Activities"
        / "Assignment_Solution_Activities.txt"
    )
    if not activities_file.exists():
        log(f"  ABORT: ovary Assignment activities not found at {activities_file}")
        return

    out_root = WORKDIR / "ovary_hrd_excluded"
    excl_project = "ovary_promoter_high_no_hrd"
    excl_maf_dir = out_root / "maf" / excl_project
    excl_maf = excl_maf_dir / f"{excl_project}.maf"
    excl_extractor_out = out_root / "extractor"

    if cosmic_signatures_path(excl_extractor_out, "SBS96").exists():
        log(f"  exists: {excl_extractor_out}")
        report_real_min_silhouette(excl_extractor_out, excl_project, "SBS96")
        return

    act = pd.read_csv(activities_file, sep="\t", index_col=0)
    if "SBS3" in act.columns:
        total = act.sum(axis=1).replace(0, pd.NA)
        sbs3_frac = (act["SBS3"] / total).fillna(0.0)
        hrd_donors = set(sbs3_frac[sbs3_frac > HRD_SBS3_THRESHOLD].index)
    else:
        log(f"  SBS3 not in attribution columns: {list(act.columns)}")
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


# --- main -------------------------------------------------------------------

def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    log(f"WORKDIR = {WORKDIR}")

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
    step7_unconstrained_comparison()
    step6_no_hyper(sbs_panc, id_panc)
    step_ovary_hrd_exclusion()

    log("DONE")


if __name__ == "__main__":
    main()
