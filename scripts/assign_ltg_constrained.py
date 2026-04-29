#!/usr/bin/env python3
"""
Constrained SigProfilerAssignment on the pan-cancer PPP-LTG aggregate matrix
using the same constrained signature set produced by the pan-cancer PPP-HTG
Extractor. Feeds the HTG vs LTG per-donor SBS96-TRC fraction figure (Panel B
of the HTG-specificity figure).

Output: outputs/sigprofiler/assignment/pancancer_promoter_low_constrained/
"""
from pathlib import Path
import sys

REPO_ROOT  = Path("/data/research/projects/trc_signatures")
WORKDIR    = REPO_ROOT / "outputs/sigprofiler"

LTG_MATRIX = WORKDIR / "sensitivity/pancancer_low/maf/pancancer_promoter_low.SBS96.all"
SIG_DB     = (WORKDIR / "pancancer/pancancer_promoter_high/extractor/SBS96"
              / "Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures"
              / "COSMIC_SBS96_Signatures.txt")
OUT_DIR    = WORKDIR / "assignment/pancancer_promoter_low_constrained"

REFERENCE  = "GRCh37"
COSMIC_VER = 3.5

if __name__ == "__main__":
    for path, label in [(LTG_MATRIX, "LTG matrix"),
                        (SIG_DB,     "constrained signature DB")]:
        if not path.exists():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    if (OUT_DIR / "Assignment_Solution").exists():
        print(f"Already done: {OUT_DIR / 'Assignment_Solution'}")
        sys.exit(0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Running constrained Assignment")
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
