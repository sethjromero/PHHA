#!/bin/bash
set -euo pipefail

# Defaults

THREADS=8

# Parse arguments

while getopts ":k:i:t:o:c:" opt; do
  case "$opt" in
    k) KVAL="$OPTARG" ;;
    i) IVAL="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    o) OUTFA="$OPTARG" ;;
    c) CLUSTER_PCT="$OPTARG" ;;
    *)
      echo "Usage: $0 -k <kval> -i <ival> -t <threads> -c <cluster_pct> -o <output.fa>"
      exit 1
      ;;
  esac
done

# Validate args

: "${KVAL:?Missing -k (k-value)}"
: "${IVAL:?Missing -i (i-value)}"
: "${CLUSTER_PCT:?Missing -c (cluster match percentage)}"
: "${OUTFA:?Missing -o (output fasta)}"

# Create temporary sequence set file

TMP_SEQS="k${KVAL}.i${IVAL}.c${CLUSTER_PCT}.seqs"

# Select contigs

echo "[INFO] Selecting contigs (k=$KVAL i=$IVAL threads=$THREADS)"

parallel --no-notice -j "$THREADS" \
	  mawk -v x="$KVAL" \''$1 >= x'\' ::: *.uniq.seqs \
	| cut -f2 \
	| perl -e '
		while (<>) {chomp; $z{$_}++;}
		while(($k,$v) = each(%z)) {print "$v\t$k\n";}
	  ' \
	| mawk -v x="$IVAL" '$1 >= x' \
	| cut -f2 \
	| mawk '{c= c + 1; print ">Contig_" c "\n" $1}' \
	| sed -e 's/NNNNNNNNNN/\t/g' \
	| cut -f1 \
	> "$TMP_SEQS"

# Cluster with cd-hit

echo "[INFO] Running cd-hit-est (c=$CLUSTER_PCT threads=$THREADS)"

cd-hit-est \
  -i "$TMP_SEQS" \
  -o "$OUTFA" \
  -M 0 \
  -T "$THREADS" \
  -c "$CLUSTER_PCT"

echo "[INFO] Done"

rm -f "$TMP_SEQS"
