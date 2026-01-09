# Assembly, Mapping, & Variant Calling

## Table of Contents

1. [De novo reference assembly](#1-de-novo-reference-assembly)  
    1a. [Defining dataset](#1a-defining-dataset)  
    1b. [Generating unique sequence sets](#1b-generating-unique-sequence-sets)



## 1. De novo reference assembly

### 1a. Defining dataset

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

Checking number of total populations (**NOTE**: sensitive to naming convention)

```sh
ls *.fastq.gz | cut -d'_' -f2 | sort -u | wc -l
```

PHHA dataset:

+ **256** individuals
+ **26** populations

Checking distribution of fastq file sizes via:

```sh
ls -l | awk '
> NR>1 {
>   s=$5
>   if (s >= 1e9)        b=">1GB"
>   else if (s >= 1e8)   b="100MB–1GB"
>   else if (s >= 1e7)   b="10–100MB"
>   else if (s >= 1e6)   b="1–10MB"
>   else                 b="<1MB"
>   c[b]++
> }
> END {
>   bins[1]=">1GB"
>   bins[2]="100MB–1GB"
>   bins[3]="10–100MB"
>   bins[4]="1–10MB"
>   bins[5]="<1MB"
>   for (i=1; i<=5; i++)
>     printf "%-12s %d\n", bins[i], c[bins[i]]+0
> }'
```

```
>1GB         0
100MB–1GB    12
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
Running the awk and perl commands:

```sh
nohup cat namelist | parallel --no-notice -j 8 "zcat {}.fastq | mawk '$AWK1' | mawk '$AWK2' | perl -e '$PERLT' > {}.uniq.seqs" 2> /dev/null &
```
