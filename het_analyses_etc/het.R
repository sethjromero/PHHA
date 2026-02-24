#### packages ####
library(stringr)
library(adegenet)
library(hierfstat)
library(pheatmap)
library(vegan) 
library(pegas)
library(poppr)

par_defaults = par(no.readonly = TRUE)

#### data prep ####
# read in gprob file (continuous dosage)
gprob_file = "entropy/estpost_files/gprob_all_combined.txt"
g = read.csv(gprob_file)
dim(g)
# 222 indivs, 7759 snps
ind_file = "entropy/PHHA.bithin.q30.i85.maxdp12.a70.maf03.012.indv"
i = read.table(ind_file)
i$id = str_extract(i$V1, "(?<=_)\\d+(?=\\.sorted\\.bam)")
i$pop = str_extract(i$V1, "(?<=_)\\w+(?=_)")
g[,1] = paste0(i$pop, "_", i$id)
pop = i$pop
rm(i, gprob_file, ind_file)
n_pops = length(unique(pop))

#### convert to genotypes ####
# options are round or strict
# round = all gprobs go to 0/1/2
# strict = 0-0.3 = 0, 0.7-1.3 = 1, 1-7-2 = 2, else NA
version = 'round'

if (version == 'round'){
  # round to nearest value
  g_mat = as.matrix(g[, 2:ncol(g)])
  g_mat = round(g_mat)
  sum(is.na(g_mat)) # should be 0
  # convert to gt format
  gt_mat = ifelse(g_mat==0, 11,
                  ifelse(g_mat==1, 12,
                         ifelse(g_mat==2, 22,
                                NA)))
  sum(is.na(gt_mat)) # should be 0
} else{
  if (version == 'strict'){
    # convert to 0/1/2 if within limits
    g_mat = as.matrix(g[, 2:ncol(g)])
    g_mat = ifelse(g_mat >= 0 & g_mat <= 0.3, 0,
                    ifelse(g_mat >= 0.7 & g_mat <= 1.3, 1,
                           ifelse(g_mat >= 1.7 & g_mat <= 2, 2,
                                  NA)))
    sum(is.na(g_mat)) # should be greater than 0
    gt_mat = ifelse(g_mat==0, 11,
                    ifelse(g_mat==1, 12,
                           ifelse(g_mat==2, 22,
                                  NA)))
    sum(is.na(gt_mat)) # should be greater than 0
  } else {
    print('invalid version')
  }
}

#### genind format ####
genind = df2genind(gt_mat, ploidy=2, pop=pop, sep="") 
summary(genind)

#### genpop format ####
genpop = genind2genpop(genind)
summary(genpop)

#### hierfstat format ####
hfs = genind2hierfstat(genind)

#### heterozygosity - over all individuals/populations ####
h_obs = summary(genind)$Hobs
h_exp = summary(genind)$Hexp
plot(NA, type = "n", xlim = c(0,1), 
     ylim = range(c(density(h_exp)$y, density(h_obs)$y)), 
     main = str_glue("Observed vs Expected Heterozygosity - {version}"), xlab = "Value", ylab = "Density")
polygon(density(h_exp), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(density(h_obs), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
legend("topright",
       c("Expected", "Observed"),
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
       border = c("red", "blue"))
# save plot
png(str_glue("het_analyses_etc/Ho_He_{version}.png"))
plot(NA, type = "n", xlim = c(0,1), 
     ylim = range(c(density(h_exp)$y, density(h_obs)$y)), 
     main = str_glue("Observed vs Expected Heterozygosity - {version}"), xlab = "Value", ylab = "Density")
polygon(density(h_exp), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(density(h_obs), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
legend("topright",
       c("Expected", "Observed"),
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
       border = c("red", "blue"))
dev.off()

#### heterozygosity - per population ####
stats = basic.stats(hfs)
h_obs_pop = stats$Ho
h_exp_pop = stats$Hs
# plot
par(mfrow = c(4,3), oma = c(0,0,2,0))
for (p in 1:n_pops){
  plot(NA, type = "n", xlim = c(0,1),
       main = str_glue("Population {colnames(h_exp_pop)[p]}"),
       ylim = range(c(density(h_obs_pop[,p])$y, density(h_exp_pop[,p])$y)), xlab = "Value", ylab = "Density")
  polygon(density(h_exp_pop[,p]), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  polygon(density(h_obs_pop[,p]), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
  legend("topright",
         c("Expected", "Observed"),
         fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
         border = c("red", "blue"))
  if (p%%12 == 0){
    mtext(str_glue("Observed vs Expected Heterozygosity by Population - {version}"), outer = TRUE)
  }
}
# save
ctr = 1
png(str_glue("het_analyses_etc/Ho_He_by_pop_{version}_pt{ctr}.png"))
par(mfrow = c(4,3), oma = c(0,0,2,0))
for (p in 1:n_pops){
  plot(NA, type = "n", xlim = c(0,1),
       main = str_glue("Population {colnames(h_exp_pop)[p]}"),
       ylim = range(c(density(h_obs_pop[,p])$y, density(h_exp_pop[,p])$y)), xlab = "Value", ylab = "Density")
  polygon(density(h_exp_pop[,p]), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  polygon(density(h_obs_pop[,p]), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
  legend("topright",
         c("Expected", "Observed"),
         fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
         border = c("red", "blue"))
  if (p%%12 == 0){
    mtext(str_glue("Observed vs Expected Heterozygosity by Population - {version}"), outer = TRUE)
    dev.off()
    ctr = ctr + 1
    png(str_glue("het_analyses_etc/Ho_He_by_pop_{version}_pt{ctr}.png"))
    par(mfrow = c(4,3), oma = c(0,0,2,0))
  }
}
dev.off()  

#### allele freqs ####
stats
df = do.call(rbind, lapply(names(stats$pop.freq), function(locus) {
  freqs = stats$pop.freq[[locus]][1, ]   # frequency of allele 1
  minor_freqs <- pmin(freqs, 1 - freqs)  # take smaller allele frequency
  data.frame(t(minor_freqs), row.names = locus)
}))
colnames(df) = colnames(stats$n.ind.samp)
write.csv(df, str_glue("het_analyses_etc/allele_freqs_{version}.csv"))
# get overall allele freqs as weighted avg of pop allele freqs
ps = summary(genind)$n.by.pop
ps = p[colnames(df)]
af = as.numeric(as.matrix(df) %*% ps) / sum(ps)


#### SFS - over all individuals/populations ####
par(par_defaults)
plot(density(af), type = "n", xlim = c(0,0.5),
     main = str_glue("SFS - {version}"),
     ylim = range(density(af)$y), xlab = "Value", ylab = "Density")
polygon(density(af), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
# save
png(str_glue("het_analyses_etc/sfs_{version}.png"))
plot(density(af), type = "n", xlim = c(0,0.5),
     main = str_glue("SFS - {version}"),
     ylim = range(density(af)$y), xlab = "Value", ylab = "Density")
polygon(density(af), col = rgb(1, 0, 0, 0.5), border = "red")
dev.off()

#### SFS - by population ####
# plot
par(mfrow = c(4,3), oma = c(0,0,2,0))
for (p in 1:n_pops){
  plot(NA, type = "n", xlim = c(0,0.5),
       main = str_glue("Population {names(ps)[p]}"),
       ylim = range(density(af)$y), xlab = "Value", ylab = "Density")
  polygon(density(df[,p]), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  if (p%%12 == 0){
    mtext(str_glue("SFS by Population - {version}"), outer = TRUE)
  }
}
# save
ctr = 1
png(str_glue("het_analyses_etc/sfs_{version}_pt{ctr}.png"))
par(mfrow = c(4,3), oma = c(0,0,2,0))
for (p in 1:n_pops){
  plot(NA, type = "n", xlim = c(0,0.5),
       main = str_glue("Population {names(ps)[p]}"),
       ylim = range(density(af)$y), xlab = "Value", ylab = "Density")
  polygon(density(df[,p]), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  if (p%%12 == 0){
    mtext(str_glue("SFS by Population - {version}"), outer = TRUE)
    dev.off()
    ctr = ctr + 1
    png(str_glue("het_analyses_etc/sfs_{version}_pt{ctr}.png"))
    par(mfrow = c(4,3), oma = c(0,0,2,0))
  }
}
dev.off() 

#### Pairwise FST ####
pwfst = pairwise.fst.dosage(g_mat, pop)
print(pwfst)
# save
write.csv(pwfst, str_glue("het_analyses_etc/pairwise_fst_{version}.csv"))
# compare to overall
print(stats$overall[7])
print(pwfst>stats$overall[7])
View(pwfst>stats$overall[7])
# plot heatmap
pheatmap(pwfst, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST",
         silent=FALSE)
# reorder by the order we did pops in for entropy
pop_order = c('ML', 'CR', 'KC', 'BA', 'WH', 'FC', 'CA', 'SH', 'MC', 'WC', 'SB',
              'PC', 'BB', 'OF', 'FV', 'SC', 'PT', 'GB', 'WF', 'BC', 'TB', 'AS',
              'SL', 'PR')
pwfst_reordered = pwfst[pop_order, pop_order]
pwfst_reordered
pheatmap(pwfst_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST",
         silent=FALSE)
# the reorder shows that the order by eye from entropy matches the clustering by pairwise fst
# we only ever see 2 pop swaps that are the closest to each other
# may be a good way to order pops by k/q quicker
# save
png(str_glue("het_analyses_etc/fst_heatmap_{version}.png"))
pheatmap(pwfst_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST",
         silent=FALSE)
dev.off()



#### Plot template ####
png(str_glue("het_analyses_etc/_{version}.png"))
dev.off()
