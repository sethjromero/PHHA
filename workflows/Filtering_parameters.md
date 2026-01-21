Initially filter on 3 things:

+ SNPs only (no indels) `--remove-indels`
+ biallelic sites only `--min-alleles 2` and `--max-alleles 2`
+ 1 SNP per contig (helps account for LD, our reads are typically >100 BP length) `--thin 100`
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
Produces file called **PHHA.bithin.q30.recode.vcf**

Make a summary file to look at individual missingness

```sh
vcftools \
--vcf PHHA.bithin.q30.recode.vcf \
--missing-indv \
--out PHHA.bithin.q30
```
Produces file called **PHHA.bithin.q30.imiss**