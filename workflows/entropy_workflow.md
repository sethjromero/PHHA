# Running entropy

##### Notes / things to talk about:

1) number of PCs to use (27->10->5)

2) currently in `/working/emcmullen/phha_temp`(no write access to romero?)

3) will add plotting code soon

## Prepare files for entropy

In `/working/romero/PHHA/vcf/final_set/` we have our filtered vcf, `PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.vcf`, along with the 012 filesFirst we need to run 2 perl scripts

```
perl vcf2mpgl_universal.pl PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.vcf
# this produces 3 files
# header_mpglPHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.txt
# indiv_ids.txt
# PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.mpgl

perl gl2genest_universal.pl PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.mpgl mean
# this produces pntest_mean_PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.txt

# the first script pulls the Phred-scaled genotype likelihoods out of the vcf file so the .mpgl file has one line per marker with the marker name and 3 values per individual (1 for each genotype) 

# the second script converts these triplets into a genotype likelihood between 0 and 2 (dosage)
```

Now copy over the .012.indv and pntest_mean file to local

In R:

```
# librarieslibrary(dplyr)library(tidyr)library(readr)library(stringr)library(MASS)library(LEA)# path to the .012.indv fileindiv_file= 'PHHA.bithin.q30.i88.maxdp12.a70.maf03.012.indv'# path to the pntest_mean filepntest_file = 'pntest_mean_PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.txt'# define function to pull population ID out of filenames
# XX_YY_NNN -> pop YY sample num NNNmakePopId <- function(fileIndv){  PopIDdf = read.table(fileIndv, sep="\t") %>%    as.data.frame() %>%    rename(All = V1) %>%    mutate(Population = str_split_i(All, '_', 2),           ID = str_split_i(All, '_', 3))  return(PopIDdf)}# define function to do pca on geno matrix
# g is ind x loci
# pulls mean dosage of each locus, calc allele freq from that
# normalize / scale all genos by mean and variance of the site
# missing vals -> 0 (mean geno)
# do PCAPCA_entropy <- function(g){  colmean = apply(g, 2, mean, na.rm = T)  normalize = matrix(nrow = nrow(g), ncol = ncol(g))  af = colmean/2  for (m in 1:length(af)){    nr = g[,m]-colmean[m]    dn = sqrt(af[m]*(1-af[m]))    normalize[,m] = nr/dn  }  normalize[is.na(normalize)] = 0  method1 = prcomp(normalize, scale. = F,center = F)  pca_df = method1$x[,1:27]  return(pca_df)}# read indiv and geno filesPopID <- makePopId(indiv_file)g <- read.table(pntest_file, header = F)# pull first 10 PCs
# transpose because pntest was loci x indv
# tack on pop IDspca_df <- PCA_entropy(t(g)) %>%  .[,1:10] %>%  cbind(PopID)# define function to get initial pop assignments
# kmeans clustering on pca space
# lda to find pcs that split clusters
# use the lda model to get prob of membership in each clusterwriteLDAfile <- function(pcaDF, k){  kCluster = kmeans(pcaDF[,1:5], k, iter.max = 10, nstart = 10, algorithm = 'Hartigan-Wong')  ldaOut = lda(x = pcaDF[,1:5], grouping = kCluster$cluster, CV = T)  write.table(round(ldaOut$posterior, 5),              file = paste('ldak', as.character(k), '.txt', sep = ''),              quote = F, row.names = F, col.names = F)}# write a file for each num pops you want to run entropy withfor (i in 2:10){  writeLDAfile(pca_df, i)}# make header for .mpgl
# pull YY_NNN back together for indiv ID
# this should be 2nd line of final .mpgl filePopID_list <- paste(PopID$Pop, PopID$ID, sep = '_')# need n indiv and n loci as the 1st lineheader <- data.frame(dims = NA, PopID_list)df <- t(header)dims <- paste(dim(g)[2], dim(g)[1], sep = " ")df[1,1] <- dimsf = paste('entropy_header.txt', sep = '')# write 2 line headerwrite.table(df,f,            sep = " ", na = "",            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Now we have ldak2.txt - ldak10.txt and entropy_header.txt locally. Copy these files back to ponderosa.

Combine the header file with the mpgl file

```
cat entropy_header.txt PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.mpgl > PHHA_entropy.mpgl
```

## Run entropy

Then run entropy for each level of k

```
conda activate entropy

nohup entropy -i PHHA_entropy.mpgl -m 1 -n 2 -k 2 -q ldak2.txt -Q 1 -l 60000 -b 10000 -t 10 -o PHHAentropyout_k2_c1.hdf5 2> temp.txt &

nohup entropy -i PHHA_entropy.mpgl -m 1 -n 2 -k 3 -q ldak3.txt -Q 1 -l 60000 -b 10000 -t 10 -o PHHAentropyout_k3_c1.hdf5 2> temp.txt &

...

nohup entropy -i PHHA_entropy.mpgl -m 1 -n 2 -k 10 -q ldak10.txt -Q 1 -l 60000 -b 10000 -t 10 -o PHHAentropyout_k10_c1.hdf5 2> temp.txt &

```

## Run estpost

Run the python script in the directory with the .hdf5 files

```
conda activate entropy_post
python estpost_entropy.py
```

Contents of estpost.py:
be sure to update the max k, path to estpost, etc

```
# assumes that hdf5 files are named XXXX_kN_cN.hdf5 where XXXX is any text
# not containing _ and then the number of populations and chain number is given

# packages
import sys
import pandas as pd
import numpy as np
import scipy as sp
import glob
import re
import random
import subprocess
import os

# display more digits
np.set_printoptions(precision=8)
pd.set_option("display.precision", 8)

# find all the hdf5 files in the directory
result = subprocess.run(["find", ".", "-name", "*hdf5"], capture_output=True, text=True, check=True)
hdf5_files = result.stdout.strip().split("\n")
hdf5_files = sorted(hdf5_files)
hdf5_files

# path to estpost (look in your conda environment)
estpost = '~/miniconda3/envs/entropy2/bin/estpost.entropy'
estpost = os.path.expanduser(estpost)

# for each hdf5 file, pull the DIC from it and write to a file
for i in range(0,len(hdf5_files)):
	f = hdf5_files[i]
	k = f.split('_')[1] #set this 
	c = f.split('_')[2].split('.hdf5')[0]
	print(k,c)
	dic = "DIC_%s_%s.txt" % (k,c)
	print(dic)
	with open(dic, "w") as out:
		subprocess.run([estpost, f, "-s", "3", "-p", "deviance"], stdout=out, check=True)

# now find all the new DIC files
dic_files = subprocess.run(["find", '.', '-name', 'DIC*'], capture_output=True, text=True, check=True).stdout.strip().split('\n')
for d in dic_files:
	subprocess.run(["cat", d], check=True)
	print("\n")

# combine DIC files to compare DIC by k and c
dic_list = []
for d in dic_files:
	k = d.split('_k')[1].split('_')[0]
	c = d.split('_c')[1].split('.txt')[0]
	result = subprocess.run(["grep", "DIC", d], capture_output=True, text=True, check=True)
	# directly extract number and convert to float
	dic = float(re.search(r"(\d+\.\d+)", result.stdout).group(0))
	dic_list.append([k, dic, c])
dic_df = pd.DataFrame(dic_list,columns=['k','DIC','chain'])
dic_df.head()

# save DIC comparison to file
dic_df.to_csv('dic_list.csv')
dic_sum = dic_df.groupby('k').describe().DIC
dic_sum.sort_values('mean')
dic_sum.to_csv('dic_sum.csv')

# for each level of k find all chains run at that level of k
# pull q for those files
for k in range(2,10):
	try:
		result=subprocess.run(['find','.','-name',f'*k{k}*.hdf5'], capture_output=True, text=True, check=True)
		hdf5_files=result.stdout.strip().split('\n')
	except subprocess.CalledProcessError:
		hdf_files=[]
	hdf5_files = [f for f in hdf5_files if f]
	if not hdf5_files:
		print(f"No .hdf5 files found for k={k}, skipping...")
		continue
	outfile=f'q{k}.txt'
	subprocess.run([estpost, *hdf5_files, "-p", "q", "-s", "0", "-o", outfile], check=True)
	print(f'wrote {outfile} for k = {k}')

# for each level of k find all chains run at that level of k
# pull Q for those files
for k in range(2,10):
	try:
		result=subprocess.run(['find','.','-name',f'*k{k}*.hdf5'], capture_output=True, text=True, check=True)
		hdf5_files=result.stdout.strip().split('\n')
	except subprocess.CalledProcessError:
		hdf_files=[]
	hdf5_files = [f for f in hdf5_files if f]
	if not hdf5_files:
		print(f"No .hdf5 files found for k={k}, skipping...")
		continue
	outfile=f'Qhet{k}.txt'
	subprocess.run([estpost, *hdf5_files, "-p", "Q", "-s", "0", "-o", outfile], check=True)
	print(f'wrote {outfile} for k = {k}')
	
# for each level of k find all chains run at that level of k
# pull MCMC diagnostics for those files
for k in range(2, 7):
	# Use find to get all matching .hdf5 files for this k
	try:
		result = subprocess.run(["find", ".", "-name", f"*k{k}*.hdf5"], capture_output=True, text=True, check=True)
		hdf5_files = result.stdout.strip().split("\n")
	except subprocess.CalledProcessError:
		hdf5_files = []
	hdf5_files = [f for f in hdf5_files if f]
	if not hdf5_files:
		print(f"No .hdf5 files found for k={k}, skipping...")
		continue
	outfile = f"MCMC_k{k}.txt"
	subprocess.run([estpost, *hdf5_files, "-p", "q", "-s", "4", "-o", outfile],check=True)
	print(f'wrote {outfile} for k = {k}')


# for each level of k find all chains run at that level of k
# pull gprobs for those files
for k in range(2,10):
	try:
		result = subprocess.run(['find','.','-name',f'*k{k}*.hdf5'], capture_output=True, text=True, check=True)
		hdf5_files=result.stdout.strip().split('\n')
	except subprocess.CalledProcessError:
		hdf5_files=[]
	hdf5_files = [f for f in hdf5_files if f]
	if not hdf5_files:
		print(f"No .hdf5 files found for k={k}, skipping...")
		continue
	outfile=f'gprob_k{k}.txt'
	subprocess.run([estpost, *hdf5_files, '-p', 'gprob', '-s', '0', '-o', outfile], check=True)
	print(f'wrote {outfile} for k = {k}')

# combine gprobs across levels of k
try:
	result=subprocess.run(['find', '.', '-name', '*.hdf5'], capture_output=True, text=True, check=True)
	hdf5_files=result.stdout.strip().split('\n')
except subprocess.CalledProcessError:
	hdf5_files=[]
hdf5_files = [f for f in hdf5_files if f]
if not hdf5_files:
	raise FileNotFoundError('no .hdf5 files found')
outfile = 'gprob_all_combined.txt'
subprocess.run([estpost,*hdf5_files,'-p','gprob','-s','0', '-o', outfile], check=True)
print(f'wrote combined genotype probabilities to {outfile}')
```

## Plot results

Pull q and Q files back to local

```
# add code here
```

# Install conda and entropy

In your home directory, get the latest miniconda installer & install

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash ~/Miniconda3-latest-Linux-x86_64.sh
```

Follow the interactive prompts to install. More documentation can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

Once you have conda installed / activated, make an environment for entropy

```
conda create -n entropy
conda activate entropy
conda config --append channels bioconda
conda install popgen-entropy
conda install gsl=2.6
# I had to specify gsl version or I would get errors

```
Separate environment for estpost

```
conda create -n entropy_post
conda activate entropy_post
conda install python
conda install pandas
conda install scipy

```