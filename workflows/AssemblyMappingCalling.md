# Assembly, Mapping, & Variant Calling

## Table of Contents

1. [De novo reference assembly](#1-de-novo-reference-assembly)  
    1a. [Characterizing dataset](#1a-characterizing-dataset)  
    1b. [Generating unique sequences](#1b-generating-unique-sequence-sets)  
    1c. [Clustering contigs for reference assembly](#1c-clustering-contigs-for-reference-assembly)  
    1d. [Indexing final assembly](#1d-indexing-final-assembly)  
2. [Mapping](#2-mapping)
3. [Variant calling](#3-variant-calling)


## 1. De novo reference assembly

### 1a. Characterizing dataset

All directories and paths in this markdown will be relative to:  
`/working/romero/PHHA/`

Starting in `./fastq`:

```sh
cd fastq
```

Checking number of total individuals

```sh
ls *.fastq.gz -1 | wc -l
```

Checking number of total populations (**NOTE**: sensitive to naming structure)

```sh
ls *.fastq.gz | cut -d'_' -f2 | sort -u | wc -l
```

PHHA dataset:

+ **256** individuals
+ **26** populations

Checking distribution of fastq file sizes

```sh
ls -l *.fastq.gz | awk '
{
  s=$5
  if (s >= 1e9)        b=">1GB"
  else if (s >= 5e8)   b="500MB–1GB"
  else if (s >= 1e8)   b="100MB-500MB"
  else if (s >= 1e7)   b="10–100MB"
  else if (s >= 1e6)   b="1–10MB"
  else                 b="<1MB"
  c[b]++
}
END {
  bins[1]=">1GB"
  bins[2]="500MB–1GB"
  bins[3]="100MB-500MB"
  bins[4]="10–100MB"
  bins[5]="1–10MB"
  bins[6]="<1MB"
  for (i=1; i<=6; i++)
    printf "%-12s %d\n", bins[i], c[bins[i]]+0
}'
```

```
>1GB         0
500MB–1GB    0
100MB-500MB  12
10–100MB     238
1–10MB       4
<1MB         2
```

Quick look at high/low extreme individuals

```sh
ls -lhS | awk 'NR>1 {print $5, $9}' | head
```

```
151M PH_PC_3.fastq.gz
111M PH_PC_4.fastq.gz
103M PH_FV_7.fastq.gz
103M PH_BA_7.fastq.gz
101M PH_PC_7.fastq.gz
100M PH_PC_9.fastq.gz
99M PH_AS_1.fastq.gz
99M PH_BC_10.fastq.gz
98M PH_BB_9.fastq.gz
97M PH_PC_6.fastq.gz
```

```sh
ls -lhS | awk '{print $5, $9}' | tail
```

```
34M PH_SL_7.fastq.gz
34M PH_SH_10.fastq.gz
25M PH_BC_12.fastq.gz
15M PH_BA_12.fastq.gz
4.4M PH_FV_3.fastq.gz
4.4M PH_PR_10.fastq.gz
3.5M PH_BA_11.fastq.gz
1.6M PH_BA_6.fastq.gz
968K PH_PM_5.fastq.gz
62K PH_WC_4.fastq.gz
```

Going ahead and building assembly on all individuals (too few outliers to bother tossing for now).

### 1b. Generating unique sequence sets

Making a list of individual IDs from the .fastq.gz files
```sh
ls *.fastq.gz | sed -e 's/.fastq.gz//g' > namelist
```

For each individual, creating a file with only the unique reads from that individual (and a count of their occurances). This will speed up future steps. Defining variables to use with awk and perl:

```sh
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'
```

Running awk and perl commands:

```sh
nohup parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs" :::: namelist &>/dev/null &
```

### 1c. Clustering contigs for reference assembly

Using the script `generateAssembly.sh` to automate this

```sh
module load cd-hit/4.6
```

From within `PHHA/fastq/`

```sh
nohup bash ../scripts/generateAssembly.sh \
    -k 2 \
    -i 2 \
    -c 0.96 \
    -t 16 \
    -o ../assembly/PHHA.k2.i2.c96.fa \
    > ../assembly/k2.i2.c96.log 2>&1 &
```

- `-k` minimum number of times a sequence must occur within a single sample to be included in clustering
- `-i` minimum number of individuals a sequence must occur in to be included in clustering
- `-c` clustering match percentage used by **cd-hit**
- `-t` number of threads to use
- `-o` output name for assembly file (contigs, FASTA format)

Ran `cd-hit` on sequence sets of `k=2` and `i=2` with a range of `c` values (**0.89 - 0.97**)

Summarizing number of contigs for each assembly based on clustering match percentage

```sh
{
  echo -e "cluster_pct\tcontig_count"
  awk '
    BEGIN{OFS="\t"}
    FNR==1{match(FILENAME,/c([0-9]{2})/,m);c=m[1];n=0}
    /^>/{n++}
    ENDFILE{print c,n}' PHHA*.fa
} > assembly_contig_counts.txt
```

```sh
cluster_pct	    contig_count
89	            763499
90	            815035
91	            878248
92	            958586
93	            1016270
94	            1072715
95	            1232428
96	            1504062
97	            1940870
```

Continuing with the assembly based on **c = 0.94**

```sh
mv PHHA.k2.i2.c94.fa PHHA_ref.fa
mv PHHA.k2.i2.c94.fa.clstr PHHA_ref.fa.clstr
```

### 1d. Indexing final assembly

Final steps of assembly

```sh
module load bwa/0.7.17-r1188
```

```sh
bwa index -p PHHA_ref PHHA_ref.fa
```


## 2. Mapping

Straightforward. Automating via the script `mapAndIndex.sh`. Look in there for details.

```sh
module load bwa/0.7.17-r1188
module load samtools/1.10
```

From within `bams` directory:

```sh
nohup bash ../scripts/mapAndIndex.sh -t 28 -r ../assembly/PHHA_ref -f "../fastq/*.fastq.gz" > mapping.log 2>/dev/null &
```

For checking ongoing status:

```sh
tail -n 1 mapping.log
```

## 3. Variant calling

From within `bams/`:

```sh
ls *.sorted.bam > bam_list.txt
```

```sh
module load bcftools/1.9
```

```sh
nohup bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f ../assembly/PHHA_ref.fa -q 30 -Q 20 -I -b bam_list.txt -o ../vcf/PHHA.bcf > ../vcf/mpileup.log 2>&1 &
```

Producing a vcf ready for filtering

```sh
nohup bcftools call -v -m -f GQ PHHA.bcf -O z -o PHHA.vcf.gz &>/dev/null &
```