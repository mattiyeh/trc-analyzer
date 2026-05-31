# SIGPROFILER.md — `run_sigprofiler.py` design notes

## Purpose

`run_sigprofiler.py` orchestrates the entire SigProfiler-based mutational
signature analysis for the TRC project: TSV→MAF conversion, matrix
generation, de novo extraction, constrained assignment, and the
sensitivity / exclusion analyses (no-hypermutator, ovary HRD exclusion).
It is the single entry point that produces every SigProfiler artifact
consumed by the manuscript figures.

Run from the repo root:

```
python3 pipeline/trc-analyzer/scripts/run_sigprofiler.py
```

The script is **idempotent and resumable**: each step checks for its
expected output before doing work, so a partial run can be restarted
without redoing finished steps.

## The KIND dispatch table

`KIND` is a top-level dict keyed by `"SBS96"` and `"ID83"` that
centralizes every parameter that differs between the two contexts —
matgen subdirectory, matrix file extension, extractor subdirectory and
k-directory prefix, COSMIC decomposed-solution paths, the
`context_type` strings expected by Extractor vs Assignment, the
`collapse_to_SBS96` flag, and `n_channels`.

**Why it exists:** the SigProfiler suite uses inconsistent argument
names and path conventions across its tools (Extractor wants
`context_type="ID83"`, Assignment wants `context_type="ID"`; Assignment
crashes if `collapse_to_SBS96=True` for indels; matgen writes to
`output/SBS/` for SBS but `output/ID/` for indels). Without KIND, every
step would carry a parallel SBS/INDEL branch and these mismatches drift
out of sync. With KIND, each step takes a `kind` argument and indexes
into one source of truth.

**How to use it:** any new step that touches matrices, extractor output,
or COSMIC decomposed solutions should accept `kind` and read parameters
via `KIND[kind][...]` rather than hard-coding paths or context strings.

## Indel MAF gotcha

When converting ICGC TSVs to MAF for indels, do **not** prepend `-` to
the ref/alt allele and do **not** decrement the start coordinate. The
SigProfilerMatrixGenerator MAF parser handles both transformations
internally. Doing them manually produces zero indel counts in every
ID83 matrix — silently, with no error.

This was the cause of the early empty-ID83-matrix bug. The fix is in
`convert_indel_tsv_to_maf`: emit the ICGC-style coordinates and alleles
unchanged.

## Constrained vs unconstrained Assignment

Step 5 runs `SigProfilerAssignment` against the **Extractor-derived
COSMIC signatures** (the `COSMIC_SBS96_Signatures.txt` from the
extractor's `Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/`),
not against the full COSMIC v3.5 catalog.

**Why constrained:** the full COSMIC database has ~80 SBS signatures,
many of which are biologically implausible in any given tumor type and
many of which are highly correlated. Unconstrained assignment overfits
— it splits TRC mutations across several look-alike signatures and
under-attributes SBS96B. Restricting the basis to the extractor's
decomposed solution forces every mutation to be explained by the
signatures the data actually selected, which gives stable per-donor
SBS96B attributions.

Step 7 (`step7_unconstrained_comparison`) runs the unconstrained
version as a sensitivity check only, so reviewers can see that the
overall conclusion does not depend on the constraint.

## The mislabeled silhouette column

`All_solutions_stat.csv` from SigProfilerExtractor has a column called
**"Minimum Stability"**, but its values are actually the **average**
per-signature silhouette across the k signatures, not the minimum. This
is a well-known SigProfiler labeling bug.

`report_real_min_silhouette` walks each `All_Solutions/SBS96_*` (or
`ID83_*`) directory, opens the per-k signature stats files, and
computes the true per-signature minimum silhouette. Always cite that
value — not the column from `All_solutions_stat.csv` — in the
manuscript and in any solution-selection discussion.

## Stale matgen cache fix

`SigProfilerMatrixGeneratorFunc` copies the input MAFs to
`<project_dir>/input/` on the first run and **reuses that cache** on
subsequent calls. If a MAF is added or modified after the first matgen
run, matgen silently ignores the new input.

`step2_matrix_generator` works around this by deleting
`<project_dir>/input/` before each matgen call **only when** the
expected matrix output is missing. This preserves resumability (we
don't re-matgen projects that already have valid matrices) while still
picking up any newly added MAF.

## Pan-cancer aggregate channel index assertion

`step3_aggregate_high` concatenates the 17 per-tumor SBS96 (or ID83)
matrices into one pan-cancer PPP-HTG matrix. Before concatenating, it
verifies that every tumor's matrix has the **same row index in the same
order**. If any matrix has a different channel ordering — which can
happen if matgen silently drops a context that has zero counts in a
tumor — concatenation by column would silently misalign mutation types
and corrupt the pan-cancer matrix.

If the assertion fails, the script exits with an error rather than
producing a corrupted aggregate. Do not weaken this check.

## Step ordering dependencies

Steps are numbered sequentially 1-23 (the renumbering table at the top of
`docs/SENSITIVITY_RESULTS_2026.md` is canonical). `main()` runs them in a
**results-first** order, not numeric order: the pan-cancer Extractor and
Assignment (which feed the figures) run as early as their dependencies
allow, and the slow / OOM-risky Extractors run last.

Actual `main()` execution order (verified against the script 2026-05-30):

    1 -> 2 -> 3 -> 4 -> 11 -> 12 -> 22 -> 10 -> 16 -> 17 -> 18 -> 13 ->
    14 -> 15 -> 5 -> 6 -> 7 -> 19 -> 20 -> 8 -> 9 -> 21 -> 23

Key dependencies that pin this order:

- **step 4 (pan-cancer Extractor) before step 11 (constrained Assignment):**
  Assignment refits against the step-4-derived COSMIC subset.
- **step 12 (unconstrained Assignment) before step 14 (MMR/POLE trim):**
  step 14 reads per-donor MMR/POLE attribution from the unconstrained fit.
- **step 16 (per-tumor Extractor) before steps 17/18 (ovary HRD exclusion):**
  the ovary HRD step reads SBS3 from ovary's *per-tumor COSMIC-decomposed
  Activities* (produced by step 16), NOT from constrained Assignment. An
  earlier version that read Assignment found 0 SBS3 and excluded no donors;
  see the bug-history comment in `run_sigprofiler.py` and
  `docs/PLANNED_SENSITIVITY_STEPS.md`.
- **step 21 (non-promoter ID83 Extractor) last:** the very sparse
  non-promoter ID83 matrix has the high-k OOM profile that crashed step 7
  ID83 on 2026-05-01; running it last means a crash doesn't block earlier
  results.

The in-code `main()` is authoritative; if it and this note ever disagree,
trust `main()`.

## Adding a new analysis step

When extending the script:

1. **Take a `kind` parameter** ("SBS96" or "ID83") and pull every
   context-specific value from `KIND[kind]`. Do not branch on `kind`
   inline; if you need a value KIND doesn't provide, add it to KIND.
2. **Maintain resumability.** Compute the expected output path first;
   if it exists and is non-empty, skip the work. Use the same
   "input MAF exists AND output matrix missing" pattern that step 2
   uses.
3. **Print timestamped progress** with `log(...)` at every meaningful
   boundary (step start, per-tumor iteration, completion). The script
   is run as an overnight job; without timestamps it's impossible to
   tell whether it's progressing or hung.
4. **Wire the new step into `main()`** in the position required by its
   data dependencies. If it depends on Assignment output, place it
   after step 5.

## Current run status

Started 2026-04-25 02:23, still running. This is the canonical 2026
SigProfiler run that will feed the updated SBS96B spectrum (Fig 3),
prevalence-by-tumor figure (Fig 4), HTG-vs-LTG comparison panels, and
the SBS96B+/SBS96B- donor classification used by the survival and
violin-plot analyses.
