#!/usr/bin/env python3
"""
Refit a 4-signature basis {SBS96A (TRC), SBS96B (companion), SBS1, SBS5}
onto both the pan-cancer PPP-HTG and PPP-LTG SBS96 matrices.

Rationale: a 2-signature {SBS96A, SBS96B} refit forces low-burden LTG donors
into a binary TRC-or-companion fit that produces bimodal-artifact fractions
(median 0.000 vs HTG median 0.571 from a 51.5/48.5 direction split).
Adding SBS1 and SBS5 gives clock-like LTG mutations a natural home, so
SBS96A picks up only mutations that genuinely look like TRC.

Builds a merged signature_database file from:
  - SBS96A, SBS96B  -> de novo signatures from HTG Extractor
  - SBS1,   SBS5    -> COSMIC v3.5 decomposed signatures from same Extractor

Outputs:
  - assignment/pancancer_promoter_high_4sig/  (HTG side)
  - assignment/pancancer_promoter_low_4sig/   (LTG side)
"""
from pathlib import Path
import sys
import pandas as pd

REPO_ROOT  = Path("/data/research/projects/trc_signatures")
WORKDIR    = REPO_ROOT / "outputs/sigprofiler"

DENOVO_SIGS  = (WORKDIR / "pancancer/pancancer_promoter_high/extractor/SBS96"
                / "Suggested_Solution/SBS96_De-Novo_Solution/Signatures"
                / "SBS96_De-Novo_Signatures.txt")
COSMIC_SIGS  = (WORKDIR / "pancancer/pancancer_promoter_high/extractor/SBS96"
                / "Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures"
                / "COSMIC_SBS96_Signatures.txt")

HTG_MATRIX = WORKDIR / "pancancer/pancancer_promoter_high/pancancer_promoter_high.SBS96.all"
LTG_MATRIX = WORKDIR / "sensitivity/pancancer_low/maf/pancancer_promoter_low.SBS96.all"

MERGED_SIGS = WORKDIR / "assignment/_4sig_basis_signatures.txt"
OUT_HTG     = WORKDIR / "assignment/pancancer_promoter_high_4sig"
OUT_LTG     = WORKDIR / "assignment/pancancer_promoter_low_4sig"

REFERENCE  = "GRCh37"
COSMIC_VER = 3.5


def build_merged_sigs():
    if MERGED_SIGS.exists():
        print(f"merged sig file exists: {MERGED_SIGS}")
        return
    denovo = pd.read_csv(DENOVO_SIGS, sep="\t")
    cosmic = pd.read_csv(COSMIC_SIGS, sep="\t")
    if not (denovo["MutationType"] == cosmic["MutationType"]).all():
        sys.exit("ERROR: MutationType columns disagree between de novo and COSMIC files")
    merged = denovo[["MutationType", "SBS96A", "SBS96B"]].copy()
    merged["SBS1"] = cosmic["SBS1"].values
    merged["SBS5"] = cosmic["SBS5"].values
    MERGED_SIGS.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(MERGED_SIGS, sep="\t", index=False)
    print(f"wrote merged 4-sig file: {MERGED_SIGS}")
    print(f"  columns: {merged.columns.tolist()}")
    for c in ["SBS96A", "SBS96B", "SBS1", "SBS5"]:
        print(f"  {c} sums to {merged[c].sum():.4f}")


def run_assignment(matrix_path, out_dir, label):
    if (out_dir / "Assignment_Solution").exists():
        print(f"  already done: {label}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  refitting {label}: {matrix_path.name}  ->  {out_dir}")
    from SigProfilerAssignment import Analyzer as Analyze
    Analyze.cosmic_fit(
        samples=str(matrix_path),
        output=str(out_dir),
        input_type="matrix",
        context_type="96",
        genome_build=REFERENCE,
        cosmic_version=COSMIC_VER,
        signature_database=str(MERGED_SIGS),
        collapse_to_SBS96=False,
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
    )


if __name__ == "__main__":
    for path, lbl in [(DENOVO_SIGS, "de novo sigs"),
                      (COSMIC_SIGS, "cosmic sigs"),
                      (HTG_MATRIX,  "HTG matrix"),
                      (LTG_MATRIX,  "LTG matrix")]:
        if not path.exists():
            sys.exit(f"ERROR: {lbl} not found: {path}")

    build_merged_sigs()
    run_assignment(HTG_MATRIX, OUT_HTG, "HTG")
    run_assignment(LTG_MATRIX, OUT_LTG, "LTG")
    print("Done.")
