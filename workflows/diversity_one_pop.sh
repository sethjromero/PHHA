#!/usr/bin/env bash
# pipeline for a single population, used as workhorse in diversity_parallel_batch.sh
# Run SAF -> realSFS -> saf2theta -> thetaStat -> extract mean pi for one population
# Usage: ./pipeline_for_one_population.sh /path/to/poplist.txt /abs/path/to/outdir

set -u
# --------- USER-CONFIGURE as needed ----------
ANGSD="${ANGSD:-$HOME/angsd}"   # path to angsd installation (ex: $HOME/angsd)
SITES="${SITES:-../vcf/pi_filtering/PHHA.bi.q30.i85.maxdp12.a70.kept.sites}"
ANC="${ANC:-../assembly/PHHA_ref.fa}"
MAXITER="${MAXITER:-200}"
FOLD="${FOLD:-1}"   # 1 for folded spectra (ancestral state unknown), 0 for unfolded (ancestral known)
P_ARG="${P_ARG:-1}"   # threads passed to angsd commands, most efficient at 1 if using de novo reference
KEEP_INTERMEDIATES="${KEEP_INTERMEDIATES:-0}"   # 0 -> remove intermediates, 1 -> keep them
# ----------------------------------------------

poplist="$1"          # path to the bam list (e.g. pop_lists/B0.BM_fmiss60.rare12.txt)
outdir="$2"           # absolute or relative output directory (will be created if missing)
logdir="${outdir}/logs"
mkdir -p "$outdir" "$logdir"

# Resolve popname: basename without .txt
popname="$(basename "$poplist" .txt)"

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

logfile="$logdir/${popname}.pipeline.log"
echo "[$(timestamp)] START $popname" >> "$logfile"

# If final pi already exists, skip
pi_out="${outdir}/${popname}.pi.txt"
if [ -f "$pi_out" ]; then
    echo "[$(timestamp)] SKIP $popname (already has $pi_out)" >> "$logfile"
    exit 0
fi

set -o pipefail

# ------------------------
# Step 1: ANGSD doSaf
# ------------------------
saf_prefix="${outdir}/${popname}"
echo "[$(timestamp)] Running ANGSD -doSaf for $popname" >> "$logfile"
"$ANGSD/angsd" -bam "$poplist" \
    -doSaf 1 \
    -anc "$ANC" \
    -GL 1 \
    -P "$P_ARG" \
    -sites "$SITES" \
    -out "$saf_prefix" \
    >> "$logfile" 2>&1
saf_idx="${saf_prefix}.saf.idx"
if [ ! -f "$saf_idx" ]; then
    echo "[$(timestamp)] ERROR: $saf_idx not created - aborting $popname" >> "$logfile"
    exit 2
fi

# ------------------------
# Step 2: realSFS 1D SFS with full convergence checking
# ------------------------
echo "[$(timestamp)] Running realSFS for $popname" >> "$logfile"
cnvlog="$outdir/${popname}.sfs.cnvrglog"
"$ANGSD/misc/realSFS" "$saf_idx" \
    -maxIter "$MAXITER" \
    -P "$P_ARG" \
    -fold "$FOLD" \
    > "${outdir}/${popname}.sfs" \
    2> "$cnvlog"
if [ ! -s "${outdir}/${popname}.sfs" ]; then
    echo "[$(timestamp)] ERROR: SFS not created for $popname" >> "$logfile"
    exit 3
fi

# Check for early EM-break
if grep -q "Breaking EM" "$cnvlog"; then
    em_line=$(grep "Breaking EM" "$cnvlog" | tail -n 1)
    em_iter=$(echo "$em_line" | grep -o "iter:[0-9]*" | cut -d: -f2)
    em_prev=$(grep "lik\[" "$cnvlog" | tail -n 1)
    em_diff=$(echo "$em_prev" | grep -o "diff=[0-9.eE+-]*" | cut -d= -f2)


echo "[$(timestamp)] Convergence reached early at iter=$em_iter with diff=$em_diff" >> "$logfile"
converged=1


else
    # Extract final diff
    last_diff=$(grep -o "diff=[0-9.eE+-]*" "$cnvlog" | tail -n 1 | cut -d= -f2)
    if [ -z "$last_diff" ]; then
        echo "[$(timestamp)] ERROR: Could not parse diff" >> "$logfile"
        exit 4
    fi
    echo "[$(timestamp)] Final SFS diff at maxIter was $last_diff" >> "$logfile"


    # Numeric threshold check
    awk -v v="$last_diff" 'BEGIN { exit !(v < 1e-6) }'
    if [ $? -eq 0 ]; then
        echo "[$(timestamp)] SFS convergence successful at maxIter (diff=$last_diff)" >> "$logfile"
        converged=1
    else
        echo "[$(timestamp)] ERROR: SFS failed to converge at maxIter (final diff=$last_diff)" >> "$logfile"
        exit 5
    fi
fi

# ------------------------
# Step 3: saf2theta
# ------------------------
echo "[$(timestamp)] Running saf2theta for $popname" >> "$logfile"
"$ANGSD/misc/realSFS" saf2theta "$saf_idx" \
    -sfs "${outdir}/${popname}.sfs" \
    -fold "$FOLD" \
    -outname "${outdir}/${popname}" \
    >> "$logfile" 2>&1
theta_idx="${outdir}/${popname}.thetas.idx"
if [ ! -f "$theta_idx" ]; then
    echo "[$(timestamp)] ERROR: $theta_idx not created" >> "$logfile"
    exit 6
fi

# ------------------------
# Step 4: thetaStat do_stat -> write to .thetas.idx.pestPG
# ------------------------
pestpg="${outdir}/${popname}.thetas.idx.pestPG"
echo "[$(timestamp)] Running thetaStat do_stat for $popname" >> "$logfile"
"$ANGSD/misc/thetaStat" do_stat "$theta_idx" -win 1 -step 1 > "$pestpg" 2>> "$logfile"
if [ ! -s "$pestpg" ]; then
    echo "[$(timestamp)] ERROR: $pestpg missing or empty" >> "$logfile"
    exit 7
fi

# ------------------------
# Step 5: extract mean pi (column 5)
# ------------------------
echo "[$(timestamp)] Extracting mean pi for $popname" >> "$logfile"
mean_pi=$(awk 'NR>1 {sum += $5; n++} END {if(n>0) printf "%.8f", sum/n; else print "NA"}' "$pestpg")
echo "$popname,$mean_pi" >> "$pi_out"
echo "[$(timestamp)] mean_pi=$mean_pi for $popname" >> "$logfile"

# ------------------------
# Step 6: cleanup intermediates if requested
# ------------------------
if [ "$KEEP_INTERMEDIATES" -eq 0 ]; then
    echo "[$(timestamp)] Cleaning up intermediates for $popname" >> "$logfile"
    rm -f "${saf_prefix}.saf" "${saf_prefix}.saf.idx" \
          "${saf_prefix}.sfs" "${saf_prefix}.sfs.cnvrglog" \
          "${saf_prefix}.thetas.gz" "${saf_prefix}.thetas.idx" \
          "${saf_prefix}.saf.gz" "${saf_prefix}.saf.pos.gz" "${saf_prefix}.arg" \
          "$pestpg"
fi

echo "[$(timestamp)] DONE $popname (pi=$mean_pi)" >> "$logfile"
exit 0
