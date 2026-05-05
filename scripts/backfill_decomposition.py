#!/usr/bin/env python3
"""Backfill Suggested_Solution/ for Extractor runs OOM-truncated at high k.

Some sensitivity Extractor runs hit OOM at the high-k transition before
producing Suggested_Solution/, but the All_Solutions/ stability cliff was
already fully resolved at a lower k. We accept k_chosen via STOPPED.txt and
want the COSMIC decomposition without redoing the full k-sweep. This script
reruns Extractor with minimum_signatures=k_chosen, maximum_signatures=k_chosen
into a sibling temp dir, then copies the resulting Suggested_Solution/ into
the original extractor_out dir. Memory at k=2..4 is tiny -- no OOM risk.

The original All_Solutions/ entries from the truncated run are preserved;
only Suggested_Solution/ (which never existed) is added. Note that the
backfilled de novo signatures may differ from All_Solutions/<sub>_<k>_Signatures/
due to NMF stochasticity (different RNG draws across runs); both are valid
representations of the k=k solution.

Usage:
  python3 backfill_decomposition.py
"""
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from run_sigprofiler import (  # noqa: E402
    EXTRACTOR_KW_BASE, KIND, WORKDIR,
    _stopped_k, log, cosmic_signatures_path,
)

# (extractor_out, matrix_path, kind)
TARGETS = [
    (
        WORKDIR / "sensitivity" / "pancancer_no_mmr_pole" / "extractor_SBS96",
        WORKDIR / "sensitivity" / "pancancer_no_mmr_pole" / "maf"
        / "pancancer_promoter_high_no_mmr_pole.SBS96.all",
        "SBS96",
    ),
    (
        WORKDIR / "sensitivity" / "pancancer_cfs" / "extractor_ID83",
        WORKDIR / "sensitivity" / "pancancer_cfs" / "maf"
        / "pancancer_cfs.ID83.all",
        "ID83",
    ),
]


def backfill_one(extractor_out, matrix_path, kind):
    extractor_out = Path(extractor_out)
    matrix_path = Path(matrix_path)
    sub = KIND[kind]["extractor_subdir"]

    if cosmic_signatures_path(extractor_out, kind).exists():
        log(f"  SKIP {extractor_out}: Suggested_Solution already exists")
        return

    k = _stopped_k(extractor_out)
    if k is None:
        log(f"  SKIP {extractor_out}: no STOPPED.txt")
        return

    if not matrix_path.exists():
        log(f"  ABORT {extractor_out}: matrix missing at {matrix_path}")
        return

    temp_out = extractor_out.parent / f"{extractor_out.name}_decomp_temp"
    if temp_out.exists():
        shutil.rmtree(temp_out)
    temp_out.mkdir(parents=True)

    log(f"  [{kind}] backfilling at k={k}")
    log(f"  [{kind}] extractor_out: {extractor_out}")
    log(f"  [{kind}] matrix:        {matrix_path}")
    log(f"  [{kind}] temp:          {temp_out}")

    from SigProfilerExtractor import sigpro as sig
    kw = dict(
        EXTRACTOR_KW_BASE,
        context_type=KIND[kind]["context_type"],
        minimum_signatures=k,
        maximum_signatures=k,
    )
    log(f"  [{kind}] params: {kw}")
    sig.sigProfilerExtractor(
        input_type="matrix",
        output=str(temp_out),
        input_data=str(matrix_path),
        **kw,
    )

    src = temp_out / sub / "Suggested_Solution"
    dst = extractor_out / sub / "Suggested_Solution"
    if not src.exists():
        log(f"  ABORT {extractor_out}: temp run did not produce Suggested_Solution")
        return
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    log(f"  [{kind}] copied Suggested_Solution -> {dst}")
    shutil.rmtree(temp_out)
    log(f"  [{kind}] cleaned {temp_out}")


def main():
    for tgt in TARGETS:
        backfill_one(*tgt)
    log("DONE")


if __name__ == "__main__":
    main()
