#!/usr/bin/env bash
# diversity_parallel_batch.sh
# Usage:
#   diversity_parallel_batch.sh pop_lists outdir [N_JOBS]
#
# pop_lists : directory containing population bam-list .txt files
# outdir    : directory for results
# N_JOBS    : number of parallel jobs (default: 26)

set -u

POP_LIST_DIR="$1"
OUTDIR="$2"
JOBS="${3:-26}"

PIPELINE="../scripts/diversity_one_pop.sh"   # your single-population script

mkdir -p "$OUTDIR"

echo "Starting ANGSD pipeline"
echo "Pop list directory: $POP_LIST_DIR"
echo "Output directory:   $OUTDIR"
echo "Parallel jobs:      $JOBS"
echo ""

# -------------------------
# Run each population in parallel
# -------------------------
parallel -j "$JOBS" \
    "$PIPELINE" {} "$OUTDIR" \
    ::: "$POP_LIST_DIR"/*.txt

echo ""
echo "All populations finished"
echo "Combining per-population pi files..."
echo ""

# -------------------------
# Combine all *.pi.txt outputs
# Each .pi.txt is already formatted as:  pop,pi
# -------------------------
SUMMARY="$OUTDIR/population_pi_summary.csv"
echo "pop,pi" > "$SUMMARY"

cat "$OUTDIR"/*.pi.txt >> "$SUMMARY"
rm -f "$OUTDIR"/*.pi.txt

echo "Wrote combined summary:"
echo "  $SUMMARY"
echo "Finished at: $(date)"