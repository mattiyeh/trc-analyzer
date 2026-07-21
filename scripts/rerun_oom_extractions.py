#!/usr/bin/env python3
"""
Re-run the two OOM-truncated Extractor runs so SigProfiler performs its OWN
selection (previously backfilled by run_sigprofiler.py:_stopped_k after the runs
were killed at high k). See docs/STABILITY_SWEEP_2026.md and
docs/STEP24_PREREGISTRATION.md-adjacent discussion.

  step 14 SBS96 : pancancer_no_mmr_pole   (dense; died at the k=9->10 peak)
  step 7  ID83  : pancancer_cfs           (sparse-per-channel; degenerate high-k)

DESIGN
------
- Runs SEQUENTIALLY. Each uses cpu=-1 (4 workers, ~44 GB peak). Two at once would
  need ~88 GB and OOM again. Do not parallelize.
- NON-DESTRUCTIVE. Writes to fresh *_rerun/ dirs. The truncated originals and
  their STOPPED.txt are left untouched for comparison. Refuses to clobber a
  non-empty rerun dir.
- IDENTICAL scientific params to the rest of the pipeline (EXTRACTOR_KW_BASE),
  so the certified selections are directly comparable and the Methods table stays
  uniform. cpu is the only resource knob and does not affect results.

REQUIRES the temporary SSD swap to be active (>= ~40 GB free swap). The script
checks and refuses to start otherwise.

PASS/FAIL is fixed in advance (stability-sweep rule): a run CERTIFIES iff
SigProfiler's own selected solution is k>=2, stable (avg silhouette >= 0.80,
min >= 0.20), and has a component with cosine >= 0.85 to canonical TRC. Whatever
it selects is reported, including step 14 dropping to k=1.

Run inside screen:
    screen -S trc_rerun
    cd /data/research/projects/trc_signatures
    python3 pipeline/trc-analyzer/scripts/rerun_oom_extractions.py
    # detach: Ctrl-a d   reattach: screen -r trc_rerun
"""
from __future__ import annotations

import glob
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/data/research/projects/trc_signatures")
SP = REPO / "outputs/sigprofiler"
REFERENCE = "GRCh37"
COSMIC_VER = 3.5

# Identical to run_sigprofiler.py:EXTRACTOR_KW_BASE (verified 2026-07-09).
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

MIN_FREE_SWAP_GB = 40
PASS_COSINE = 0.85

# (label, kind, context_type, matrix, fresh output dir, canonical sig file, canon col, class)
RUNS = [
    dict(
        label="step14_SBS96_no_mmr_pole",
        kind="SBS96",
        context_type="SBS96",
        matrix=SP / "sensitivity/pancancer_no_mmr_pole/maf/pancancer_promoter_high_no_mmr_pole.SBS96.all",
        out=SP / "sensitivity/pancancer_no_mmr_pole/extractor_SBS96_rerun",
        canon=SP / "pancancer/pancancer_promoter_high/extractor/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
        canon_col="SBS96A",
        cls="C>G",
    ),
    dict(
        label="step7_ID83_cfs",
        kind="ID83",
        context_type="ID83",
        matrix=SP / "sensitivity/pancancer_cfs/maf/pancancer_cfs.ID83.all",
        out=SP / "sensitivity/pancancer_cfs/extractor_ID83_rerun",
        canon=SP / "pancancer/pancancer_promoter_high_indel/extractor/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt",
        canon_col="ID83B",
        cls=None,
    ),
]


def log(msg: str) -> None:
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def free_swap_gb() -> float:
    with open("/proc/meminfo") as fh:
        m = {l.split(":")[0]: l.split()[1] for l in fh}
    return int(m.get("SwapFree", "0")) / 1024 / 1024


def cosine(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def post_check(run: dict) -> None:
    """Read the freshly SigProfiler-SELECTED solution and report the verdict."""
    try:
        sub = run["kind"]
        sig_glob = str(run["out"] / sub / "Suggested_Solution" /
                       f"{sub}_De-Novo_Solution" / "Signatures" / "*De-Novo_Signatures.txt")
        hits = glob.glob(sig_glob)
        if not hits:
            log(f"  {run['label']}: NO Suggested_Solution written -- did it truncate again? "
                f"check {run['out']}")
            return
        sel = pd.read_csv(hits[0], sep="\t")
        sel = sel.set_index(sel.columns[0])
        sel_k = sel.shape[1]

        canon = pd.read_csv(run["canon"], sep="\t")
        canon = canon.set_index(canon.columns[0])[run["canon_col"]].reindex(sel.index)
        best = max(sel.columns, key=lambda c: cosine(sel[c], canon))
        cos = cosine(sel[best], canon)

        stat = pd.read_csv(run["out"] / sub / "All_solutions_stat.csv")
        srow = stat[stat["Signatures"].astype(str).str.contains(r"\*")]
        avg = float(srow["Stability (Avg Silhouette)"].iloc[0]) if len(srow) else float("nan")
        mn = float(srow["Minimum Stability"].iloc[0]) if len(srow) else float("nan")

        stable = (avg >= 0.80) and (mn >= 0.20) and (avg + mn >= 1.00)
        if sel_k >= 2 and stable and cos >= PASS_COSINE:
            verdict = "CERTIFIED PASS"
        elif sel_k == 1:
            verdict = "k=1 blend (no conclusion)"
        elif not stable:
            verdict = f"selected k={sel_k} UNSTABLE"
        else:
            verdict = f"selected k={sel_k}, cosine {cos:.3f} < 0.85"
        log(f"  {run['label']} RESULT: selected k={sel_k}, avg {avg:.2f}/min {mn:.2f}, "
            f"cosine-to-canonical {cos:.4f}  ->  {verdict}")
    except Exception as exc:  # a post-check failure must not mask a successful run
        log(f"  {run['label']}: post-check error ({exc}); inspect {run['out']} by hand")


def main() -> None:
    fs = free_swap_gb()
    log(f"free swap = {fs:.1f} GB (need >= {MIN_FREE_SWAP_GB})")
    if fs < MIN_FREE_SWAP_GB:
        log("ABORT: insufficient free swap. Create /swapfile_trc first (see chat).")
        sys.exit(1)

    from SigProfilerExtractor import sigpro as sig

    for run in RUNS:
        if not run["matrix"].exists():
            log(f"ABORT: matrix missing: {run['matrix']}")
            sys.exit(1)
        if run["out"].exists() and any(run["out"].iterdir()):
            log(f"SKIP {run['label']}: output dir already non-empty ({run['out']}). "
                f"Remove it to force a re-run.")
            continue

        run["out"].mkdir(parents=True, exist_ok=True)
        kw = dict(EXTRACTOR_KW_BASE, context_type=run["context_type"])
        log(f"START {run['label']}")
        log(f"  matrix : {run['matrix']}")
        log(f"  out    : {run['out']}")
        log(f"  params : {kw}")
        t0 = datetime.now()
        try:
            sig.sigProfilerExtractor("matrix", str(run["out"]), str(run["matrix"]), **kw)
        except Exception as exc:
            log(f"  {run['label']} FAILED: {exc}")
            log("  stopping; the other run is not attempted after a failure.")
            sys.exit(1)
        log(f"DONE  {run['label']} in {datetime.now() - t0}")
        post_check(run)

    log("ALL RUNS COMPLETE. Remember to swapoff + rm /swapfile_trc when finished.")


if __name__ == "__main__":
    main()
