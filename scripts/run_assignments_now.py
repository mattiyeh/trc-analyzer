#!/usr/bin/env python3
"""
One-off helper: run pipeline steps 11 (constrained Assignment) and 12
(unconstrained Assignment) RIGHT NOW, in parallel with the main
run_sigprofiler.py that is currently grinding through Extractor steps.

This works because:
  - Step 11/12 only depend on step 4's Extractor-derived COSMIC DB,
    which already exists on disk.
  - Each Assignment call writes to its own output dir (no race with
    the running Extractor steps, which write under pancancer/sensitivity).
  - The pipeline is resumable: when it eventually reaches steps 11-12,
    `if (out_dir / "Assignment_Solution").exists(): return` skips them.

Output goes to outputs/sigprofiler/run_assignments_now.log via tee from
the calling shell.
"""
from __future__ import annotations
import sys
from pathlib import Path

# Add the pipeline scripts dir to sys.path so we can import the canonical
# helpers from run_sigprofiler instead of duplicating them.
_SCRIPTS = Path("/data/research/projects/trc_signatures/pipeline/trc-analyzer/scripts")
sys.path.insert(0, str(_SCRIPTS))

from run_sigprofiler import (  # noqa: E402
    log,
    cosmic_signatures_path,
    step11_constrained_assignment,
    step12_unconstrained_assignment,
    WORKDIR,
    PANCANCER_LABEL,
    PANCANCER_INDEL_LABEL,
)


def main():
    sbs_extractor_out = WORKDIR / "pancancer" / PANCANCER_LABEL       / "extractor"
    id_extractor_out  = WORKDIR / "pancancer" / PANCANCER_INDEL_LABEL / "extractor"

    sbs_cosmic_db = cosmic_signatures_path(sbs_extractor_out, "SBS96")
    id_cosmic_db  = cosmic_signatures_path(id_extractor_out,  "ID83")

    if not sbs_cosmic_db.exists():
        sys.exit(f"ERROR: SBS COSMIC DB missing: {sbs_cosmic_db}")
    if not id_cosmic_db.exists():
        sys.exit(f"ERROR: ID  COSMIC DB missing: {id_cosmic_db}")

    log(f"SBS COSMIC DB: {sbs_cosmic_db}")
    log(f"ID  COSMIC DB: {id_cosmic_db}")
    log("--- starting step 11 (constrained Assignment) ---")
    step11_constrained_assignment(sbs_cosmic_db, id_cosmic_db)
    log("--- starting step 12 (unconstrained Assignment) ---")
    step12_unconstrained_assignment()
    log("DONE")


if __name__ == "__main__":
    main()
