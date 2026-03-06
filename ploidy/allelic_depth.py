# set up conda env in terminal
# conda create -n phha
# conda activate phha
# conda install python pandas matplotlib seaborn cyvcf2

# in r console
#install.packages("reticulate")
#library(reticulate)
#use_condaenv(phha) # name of conda env
#reticulate::py_config() # check that it is the right one

# in terminal copy over the current vcf and the vcf with the 2 removed pops


# running first command also runs reticulate::repl_python() so you can use python in rstudio
# exit will return to r

# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cyvcf2 import VCF
import os

os.getcwd()
os.chdir('Documents/GitHub/PHHA/ploidy')

# set up paths
vers = 'current'

if vers == 'current':
  vcf_file = 'PHHA.bithin.q30.i85.maxdp12.a70.maf03.recode.vcf'
  pop_file = 'PHHA.bithin.q30.i85.maxdp12.a70.maf03.012.indv'
elif vers == 'old':
  vcf_file = 'PHHA.bithin.q30.i88.maxdp12.a70.maf03.recode.vcf'
  pop_file = 'PHHA.bithin.q30.i88.maxdp12.a70.maf03.012.indv'
else:
  print("invalid version")

csv_file = f'allele_fraction_by_population_{vers}.csv'

# read in pop data
pop_df = pd.read_csv(pop_file, sep='\t', header = None)
pop_df.head()
pop_df['pop'] = [x.split('_')[1] for x in pop_df[0]]
#pop_df['sample'] = [x.split('.')[0] for x in pop_df[0]]
pop_df.head()
sample_to_pop = dict(zip(pop_df[0], pop_df['pop']))
sample_to_pop

# read in vcf
vcf = VCF(vcf_file)
samples = vcf.samples
sample_pop = [sample_to_pop.get(s, 'Unknown') for s in samples]

# calc allele fractions for hets
allele_data = []

# for every snp
for variant in vcf:
  # skip multiallelic sites
  if len(variant.ALT) > 1:
    continue
  types = variant.gt_types # 0 = 0/0, 1 = het, 2 = unknown, 3 = 1/1 for each indiv
  ads = variant.format('AD') # num ref reads, num alt reads for each indiv
  dps = variant.format('DP') # total num reads for each indiv
  # for each individual
  for i in range(len(types)):
    gt = types[i]
    alt_dep = ads[i][1]
    tot_dep = dps[i]
    # qc
    if tot_dep > 0:
      alt_frac = alt_dep/tot_dep
    else:
      alt_frac = np.nan
    if gt != 1:
      alt_frac = np.nan
    # save that individuals results at that snp
    allele_data.append({
      'sample': samples[i],
      'population': sample_pop[i],
      'allele_fraction': alt_frac,
      'depth': tot_dep,
      'snp_id': f"{variant.CHROM}_{variant.POS}"})
      
df = pd.DataFrame(allele_data)    
df.head()
df['allele_fraction'] = [x[0] if not np.isnan(x) else np.nan for x in df['allele_fraction']]
df.head()
df.to_csv(csv_file, index=False)

min_depth = 1
plot_file = f'allele_fraction_by_population_{vers}_min_dp_{min_depth}.png'
# aggregate by population and plot
pops = df['population'].unique()
pops

df2 = df

for p in pops:
  in_p = df2['population'] == p
  pdf = df2[in_p.values]
  pop_dfs.append(pdf)
  plt.figure()
  sns.kdeplot(pdf['allele_fraction'], fill=True, color='purple', alpha=0.5)
  plt.vlines([0, 0.25, 0.5, 0.75, 1], 0, 2, 'k')
  plt.title(f'{p}')
  plt.show()
  plt.savefig(f'allele_fraction_pop_{p}_{vers}_min_dp_{min_depth}.png')
