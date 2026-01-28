Initially filtering on 3 things:

+ SNPs only (no indels) `--remove-indels`
+ biallelic sites only `--min-alleles 2` and `--max-alleles 2`
+ 1 SNP per contig (helps account for LD, our reads are typically <100 BP length) `--thin 100`
+ **OPTIONAL (but I typically do)** - filter on call quality. Our call qualities are typically very good so it doesn't get rid of much but might as well cut out the worst stuff `--minQ 30`


```sh
vcftools \
--gzvcf PHHA.vcf.gz \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--thin 100 \
--minQ 30 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out PHHA.bithin.q30
```
Produces file **PHHA.bithin.q30.recode.vcf**

Make a summary file to look at individual missingness

```sh
vcftools \
--vcf PHHA.bithin.q30.recode.vcf \
--missing-indv \
--out PHHA.bithin.q30
```
Produces file called **PHHA.bithin.q30.imiss**

Based on missingness, removing individuals with >88% missigness on the minimally-filtered SNP set defined above)

```sh
awk '$5 > 0.88 {print $1}' PHHA.bithin.q30.imiss | tail -n +2 > indmiss88.txt
```

Now filter out those identified individuals and make a downsized VCF via the following:

```sh
vcftools \
--vcf PHHA.bithin.q30.recode.vcf \
--remove indmiss88.txt \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out PHHA.bithin.q30.i88
```

Next taking a look at depth, particularly to set a maximum loci depth. 

```sh
vcftools \
--vcf PHHA.bithin.q30.i88.recode.vcf \
--site-mean-depth \
--out PHHA.bithin.q30.i88
```
Produces file called **PHHA.bithin.q30.i88.ldepth.mean**

```sh
vcftools \
--vcf PHHA.bithin.q30.i88.recode.vcf \
--positions sites_maxdp12.tsv \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out PHHA.bithin.q30.i88.maxdp12
```