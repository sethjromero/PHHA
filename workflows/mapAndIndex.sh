#!/bin/bash

while getopts ":t:r:f:" opt; do
    case ${opt} in
        t )
            THREADS="$OPTARG"
            ;;
        r )
            REFERENCE="$OPTARG"
            ;;
        f )
            FASTQS="$OPTARG"
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

# check that all required options were provided
if [ -z "$THREADS" ] || [ -z "$REFERENCE" ] || [ -z "$FASTQS" ]; then
    echo "Usage: $0 -t <threads> -r <reference> -f <fastq_glob>"
    exit 1
fi

ctr=1
fTotal=$(ls $FASTQS -1 | wc -l)

for file in $FASTQS
	do
	if test -f "$file"
	then
		fPrefix=$(echo "$file" | sed 's|.*/||' | cut -f1 -d ".")
		echo mapping and sorting individual "$ctr" of "$fTotal"
		bwa mem -t "$THREADS" "$REFERENCE" "$file" | \
		samtools view -u -b - | \
		samtools sort -l0 -@"$THREADS" -o "$fPrefix.sorted.bam"
		((ctr++))
	fi
done
for sBam in *.sorted.bam
	do
	if test -f "$sBam"
	then
		samtools index -@"$THREADS" "$sBam"
	fi
done
