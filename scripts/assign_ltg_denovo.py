#!/usr/bin/env python3
"""
Refit the pan-cancer PPP-HTG Extractor's de novo signatures (SBS96A, SBS96B)
onto the pan-cancer PPP-LTG SBS96 matrix.

SBS96A = SBS96-TRC (per De_Novo_map_to_COSMIC: SBS40a 27.5% top component).
SBS96B = clock+repair companion.

For Panel B (HTG vs LTG per-donor SBS96-TRC fraction), compare:
  - HTG side: existing  pancancer/.../SBS96_De-Novo_Solution/Activities/
              SBS96_De-Novo_Activities_refit.txt  (658 donors)
  - LTG side: this script's output                                (~569 donors)

Output: outputs/sigprofiler/assignment/pancancer_promoter_low_denovo/
"""
from pathlib import Path
import sys

REPO_ROOT  = Path("/data/research/projects/trc_signatures")
WORKDIR    = REPO_ROOT / "outputs/sigprofiler"

LTG_MATRIX = WORKDIR / "sensitivity/pancancer_low/maf/pancancer_promoter_low.SBS96.all"
SIG_DB     = (WORKDIR / "pancancer/pancancer_promoter_high/extractor/SBS96"
              / "Suggested_Solution/SBS96_De-Novo_Solution/Signatures"
              / "SBS96_De-Novo_Signatures.txt")
OUT_DIR    = WORKDIR / "assignment/pancancer_promoter_low_denovo"

REFERENCE  = "GRCh37"
COSMIC_VER = 3.5

if __name__ == "__main__":
    for path, label in [(LTG_MATRIX, "LTG matrix"),
                        (SIG_DB,     "de novo signature DB")]:
        if not path.exists():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    if (OUT_DIR / "Assignment_Solution").exists():
        print(f"Already done: {OUT_DIR / 'Assignment_Solution'}")
        sys.exit(0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Refitting de novo HTG signatures onto LTG matrix")
    print(f"  matrix:   {LTG_MATRIX}")
    print(f"  sig DB:   {SIG_DB}")
    print(f"  output:   {OUT_DIR}")

    from SigProfilerAssignment import Analyzer as Analyze
    Analyze.cosmic_fit(
        samples=str(LTG_MATRIX),
        output=str(OUT_DIR),
        input_type="matrix",
        context_type="96",
        genome_build=REFERENCE,
        cosmic_version=COSMIC_VER,
        signature_database=str(SIG_DB),
        collapse_to_SBS96=False,
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
    )
    print("Done.")
