#!/usr/bin/env python3
"""End-to-end SigProfiler workflow for the TRC pan-cancer study.

Pipeline:
  1. Convert each ICGC SBS TSV -> minimal MAF (icgc_donor_id as sample ID).
  2. Run SigProfilerMatrixGenerator on every (tumor x group) -> SBS96 matrix.
  3. Aggregate the 17 tumor-type "promoter_high" SBS96 matrices into a single
     pan-cancer high cohort (the primary input for de novo extraction).
  4. Run SigProfilerExtractor de novo on the pan-cancer high aggregate.
  5. Run constrained SigProfilerAssignment on every other (tumor x group)
     using the Extractor-derived COSMIC SBS96 signature database.
  6. No-hypermutator sensitivity analysis: drop samples with mutation counts
     above mean + 2 SD and rerun Extractor on the trimmed cohort.

The script is resumable -- each step skips work whose outputs already exist.
Run it from anywhere; all paths are absolute.
"""

import csv
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


# --- Configuration ----------------------------------------------------------

REPO_ROOT     = Path("/data/research/projects/trc_signatures")
ICGC_OUTPUTS  = REPO_ROOT / "icgc-datasets-outputs__2026.03.15_00.10.17"
TIMESTAMP_TAG = "2026.03.15_00.10.17"
WORKDIR       = REPO_ROOT / "sigprofiler_run"
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

PRIMARY_GROUP   = "promoter_high"
PANCANCER_LABEL = "pancancer_promoter_high"

# Minimal MAF columns understood by SigProfilerMatrixGenerator.
MAF_COLUMNS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build",
    "Chromosome", "Start_position", "End_position", "Strand",
    "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
]

# Extractor parameters per the best-practices guide. KL divergence is the
# default objective in SigProfilerExtractor's underlying NMF and is not
# exposed as a top-level parameter.
EXTRACTOR_KW = dict(
    reference_genome=REFERENCE,
    opportunity_genome=REFERENCE,
    context_type="SBS96",
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


def convert_tsv_to_maf(tsv_path, maf_path):
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


def step1_convert_all():
    log("Step 1: TSV -> MAF")
    for tumor in TUMORS:
        for subdir, suffix, label in GROUPS:
            tsv = icgc_tsv_path(tumor, subdir, suffix)
            project = f"{tumor}_{label}"
            project_dir = WORKDIR / "maf" / project
            maf = project_dir / f"{project}.maf"

            if maf.exists():
                continue
            if not tsv.exists():
                log(f"  MISSING TSV: {tsv}")
                continue
            project_dir.mkdir(parents=True, exist_ok=True)
            n = convert_tsv_to_maf(tsv, maf)
            log(f"  {project}: {n} SBS rows -> {maf.name}")


# --- Step 2: SigProfilerMatrixGenerator -------------------------------------

def matgen_sbs96_path(project):
    return WORKDIR / "maf" / project / "output" / "SBS" / f"{project}.SBS96.all"


def step2_matrix_generator():
    log("Step 2: SigProfilerMatrixGenerator (SBS96)")
    from SigProfilerMatrixGenerator.scripts import (
        SigProfilerMatrixGeneratorFunc as matGen,
    )
    for tumor in TUMORS:
        for _, _, label in GROUPS:
            project = f"{tumor}_{label}"
            project_dir = WORKDIR / "maf" / project
            sbs96 = matgen_sbs96_path(project)
            if sbs96.exists():
                continue
            if not (project_dir / f"{project}.maf").exists():
                continue
            log(f"  matgen {project}")
            try:
                matGen.SigProfilerMatrixGeneratorFunc(
                    project, REFERENCE, str(project_dir),
                    plot=False, exome=False, bed_file=None,
                    chrom_based=False, tsb_stat=False, seqInfo=False,
                )
            except Exception as exc:
                log(f"    SKIP {project}: {exc}")


# --- Step 3: aggregate pan-cancer high --------------------------------------

def pancancer_matrix_path():
    return (
        WORKDIR / "pancancer" / PANCANCER_LABEL / f"{PANCANCER_LABEL}.SBS96.all"
    )


def step3_aggregate_high():
    log(f"Step 3: aggregate '{PRIMARY_GROUP}' across {len(TUMORS)} tumor types")
    out_path = pancancer_matrix_path()
    if out_path.exists():
        log(f"  exists: {out_path}")
        return out_path

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for tumor in TUMORS:
        project = f"{tumor}_{PRIMARY_GROUP}"
        sbs = matgen_sbs96_path(project)
        if not sbs.exists():
            log(f"  missing matrix for {project} -- skipping")
            continue
        df = pd.read_csv(sbs, sep="\t", index_col=0)
        if df.shape[1] == 0:
            log(f"  empty matrix for {project}")
            continue
        # ICGC donor IDs are globally unique, but prefix by tumor for clarity.
        df.columns = [f"{tumor}__{c}" for c in df.columns]
        frames.append((tumor, df))
        log(f"  {project}: {df.shape[1]} samples")

    if not frames:
        sys.exit("No high matrices found -- cannot aggregate")

    # All matrices must share the same 96 trinucleotide channels in the same
    # order; otherwise pd.concat would silently misalign rows.
    ref_tumor, ref_df = frames[0]
    for tumor, df in frames[1:]:
        if not df.index.equals(ref_df.index):
            sys.exit(
                f"SBS96 channel mismatch: '{tumor}' row index differs from "
                f"'{ref_tumor}'. Refusing to aggregate -- inputs were not "
                f"produced by a consistent MatrixGenerator run."
            )

    agg = pd.concat([df for _, df in frames], axis=1).fillna(0).astype(int)
    agg.index.name = "MutationType"
    agg.to_csv(out_path, sep="\t")
    log(f"  wrote {out_path} ({agg.shape[1]} samples, {agg.shape[0]} channels)")
    return out_path


# --- Step 4: Extractor de novo ---------------------------------------------

def cosmic_signatures_path(extractor_output):
    return (
        Path(extractor_output)
        / "SBS96" / "Suggested_Solution"
        / "COSMIC_SBS96_Decomposed_Solution" / "Signatures"
        / "COSMIC_SBS96_Signatures.txt"
    )


def run_extractor(matrix_path, output_dir, label):
    log(f"Extractor: de novo on {label}")
    cosmic_db = cosmic_signatures_path(output_dir)
    if cosmic_db.exists():
        log(f"  exists: {cosmic_db}")
        return cosmic_db
    from SigProfilerExtractor import sigpro as sig
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    log(f"  params: {EXTRACTOR_KW}")
    sig.sigProfilerExtractor(
        "matrix", str(output_dir), str(matrix_path),
        **EXTRACTOR_KW,
    )
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
    """Search a per-k directory for any tabular file that contains a column of
    per-signature silhouette / stability values, and return the minimum."""
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
                # A per-signature silhouette column has one row per signature
                # and values bounded in [-1, 1].
                if 0 < len(vals) <= k and vals.between(-1.0, 1.0).all():
                    m = float(vals.min())
                    if best is None or m < best:
                        best = m
    return best


def report_real_min_silhouette(extractor_output, label):
    log(f"Real per-signature minimum silhouette ({label}):")
    sbs_dir = Path(extractor_output) / "SBS96" / "All_Solutions"
    if not sbs_dir.exists():
        log("  (no All_Solutions/ directory)")
        return
    rows = []
    for k_dir in sorted(sbs_dir.glob("SBS96_*_Signatures"),
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


# --- Step 5: constrained Assignment ----------------------------------------

def step5_assignment(cosmic_db):
    log("Step 5: constrained Assignment")
    if not cosmic_db.exists():
        log(f"  ABORT: signature DB missing: {cosmic_db}")
        return
    from SigProfilerAssignment import Analyzer as Analyze

    def assign(label, matrix_path, out_dir):
        if not matrix_path.exists():
            return
        if (out_dir / "Assignment_Solution").exists():
            return
        log(f"  assign {label}")
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            Analyze.cosmic_fit(
                samples=str(matrix_path),
                output=str(out_dir),
                input_type="matrix",
                context_type="96",
                genome_build=REFERENCE,
                cosmic_version=COSMIC_VER,
                signature_database=str(cosmic_db),
                collapse_to_SBS96=True,
                nnls_add_penalty=0.05,
                nnls_remove_penalty=0.01,
            )
        except Exception as exc:
            log(f"    SKIP {label}: {exc}")

    # Pan-cancer aggregate (constrained refit on the same input the Extractor
    # ran on -- methodological transparency, per the best-practices guide).
    assign(PANCANCER_LABEL, pancancer_matrix_path(),
           WORKDIR / "assignment" / PANCANCER_LABEL)

    # Per-(tumor, group) cohorts, including each tumor's high group.
    for tumor in TUMORS:
        for _, _, label in GROUPS:
            project = f"{tumor}_{label}"
            assign(project, matgen_sbs96_path(project),
                   WORKDIR / "assignment" / project)


# --- Step 7: unconstrained Assignment (sensitivity comparison) -------------
# Constrained refit (Step 5) is the primary result. This adds the unconstrained
# fit against the full bundled COSMIC database, the comparison §7 of the
# best-practices guide expects every signatures paper to report.

def step7_unconstrained_comparison():
    log("Step 7: unconstrained Assignment on pan-cancer high (sensitivity)")
    matrix_path = pancancer_matrix_path()
    if not matrix_path.exists():
        log(f"  ABORT: pan-cancer matrix missing: {matrix_path}")
        return
    out_dir = WORKDIR / "assignment" / f"{PANCANCER_LABEL}_unconstrained"
    if (out_dir / "Assignment_Solution").exists():
        log(f"  exists: {out_dir}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    from SigProfilerAssignment import Analyzer as Analyze
    log(f"  assign {PANCANCER_LABEL} (full COSMIC database, signature_database=None)")
    Analyze.cosmic_fit(
        samples=str(matrix_path),
        output=str(out_dir),
        input_type="matrix",
        context_type="96",
        genome_build=REFERENCE,
        cosmic_version=COSMIC_VER,
        signature_database=None,
        collapse_to_SBS96=True,
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
    )


# --- Step 6: no-hypermutator sensitivity analysis --------------------------

def step6_no_hyper(panc_matrix):
    log("Step 6: no-hypermutator sensitivity analysis")
    no_hyper_dir    = WORKDIR / "pancancer" / f"{PANCANCER_LABEL}_no_hyper"
    no_hyper_matrix = no_hyper_dir / f"{PANCANCER_LABEL}_no_hyper.SBS96.all"
    extractor_out   = no_hyper_dir / "extractor"

    if cosmic_signatures_path(extractor_out).exists():
        log(f"  exists: {extractor_out}")
        report_real_min_silhouette(extractor_out, "pancancer_no_hyper")
        return

    no_hyper_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(panc_matrix, sep="\t", index_col=0)
    counts = df.sum(axis=0)
    mu, sd = float(counts.mean()), float(counts.std())
    cutoff = mu + 2 * sd
    keep = counts[counts <= cutoff].index
    drop = sorted(set(counts.index) - set(keep))
    log(f"  cohort: n={len(counts)}, mean={mu:.0f}, sd={sd:.0f}, "
        f"cutoff={cutoff:.0f}")
    log(f"  dropping {len(drop)} hypermutator(s): {drop}")

    df[keep].to_csv(no_hyper_matrix, sep="\t")
    log(f"  filtered matrix: {no_hyper_matrix} ({len(keep)} samples)")

    run_extractor(no_hyper_matrix, extractor_out, "pancancer_no_hyper")
    report_real_min_silhouette(extractor_out, "pancancer_no_hyper")


# --- main -------------------------------------------------------------------

def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    log(f"WORKDIR = {WORKDIR}")

    step1_convert_all()
    step2_matrix_generator()
    panc_matrix = step3_aggregate_high()

    extractor_out = WORKDIR / "pancancer" / PANCANCER_LABEL / "extractor"
    cosmic_db = run_extractor(panc_matrix, extractor_out, PANCANCER_LABEL)
    report_real_min_silhouette(extractor_out, PANCANCER_LABEL)

    step5_assignment(cosmic_db)
    step7_unconstrained_comparison()
    step6_no_hyper(panc_matrix)

    log("DONE")


if __name__ == "__main__":
    main()
