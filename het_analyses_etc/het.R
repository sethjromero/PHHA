#### packages ####
library(stringr)
library(adegenet)
library(hierfstat)
library(pheatmap)
library(vegan) 
library(pegas)
library(poppr)
library(sf)
library(ggplot2)
library(dplyr)
library(maps)
library(ggmap) 
library(tidyr)

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
version = 'strict'

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
     ylim = range(c(density(h_exp, na.rm = TRUE)$y, density(h_obs, na.rm = TRUE)$y)), 
     main = str_glue("Observed vs Expected Heterozygosity - {version}"), xlab = "Value", ylab = "Density")
polygon(density(h_exp, na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(density(h_obs, na.rm = TRUE), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
legend("topright",
       c("Expected", "Observed"),
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
       border = c("red", "blue"))
# save plot
png(str_glue("het_analyses_etc/Ho_He_{version}.png"))
plot(NA, type = "n", xlim = c(0,1), 
     ylim = range(c(density(h_exp, na.rm = TRUE)$y, density(h_obs, na.rm = TRUE)$y)), 
     main = str_glue("Observed vs Expected Heterozygosity - {version}"), xlab = "Value", ylab = "Density")
polygon(density(h_exp, na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(density(h_obs, na.rm = TRUE), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
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
       ylim = range(c(density(h_obs_pop[,p], na.rm = TRUE)$y, density(h_exp_pop[,p], na.rm = TRUE)$y)), xlab = "Value", ylab = "Density")
  polygon(density(h_exp_pop[,p], na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  polygon(density(h_obs_pop[,p], na.rm = TRUE), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
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
       ylim = range(c(density(h_obs_pop[,p], na.rm = TRUE)$y, density(h_exp_pop[,p], na.rm = TRUE)$y)), xlab = "Value", ylab = "Density")
  polygon(density(h_exp_pop[,p], na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
  polygon(density(h_obs_pop[,p], na.rm = TRUE), col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
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
ps = ps[colnames(df)]
af = as.numeric(as.matrix(df) %*% ps) / sum(ps)


#### SFS - over all individuals/populations ####
par(par_defaults)
plot(density(af, na.rm = TRUE), type = "n", xlim = c(0,0.5),
     main = str_glue("SFS - {version}"),
     ylim = range(density(af, na.rm = TRUE)$y), xlab = "Value", ylab = "Density")
polygon(density(af, na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
# save
png(str_glue("het_analyses_etc/sfs_{version}.png"))
plot(density(af, na.rm = TRUE), type = "n", xlim = c(0,0.5),
     main = str_glue("SFS - {version}"),
     ylim = range(density(af, na.rm = TRUE)$y), xlab = "Value", ylab = "Density")
polygon(density(af, na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")
dev.off()

#### SFS - by population ####
# plot
par(mfrow = c(4,3), oma = c(0,0,2,0))
for (p in 1:n_pops){
  plot(NA, type = "n", xlim = c(0,0.5),
       main = str_glue("Population {names(ps)[p]}"),
       ylim = range(density(af, na.rm = TRUE)$y), xlab = "Value", ylab = "Density")
  polygon(density(df[,p], na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
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
       ylim = range(density(af, na.rm = TRUE)$y), xlab = "Value", ylab = "Density")
  polygon(density(df[,p], na.rm = TRUE), col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
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

#### Nei's D ####
rownames(genpop@tab) # order of the pops for dist matrix
gen_dist = as.data.frame(as.matrix(dist.genpop(genpop)))
pheatmap(gen_dist, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Nei's D",
         silent=FALSE)
gen_dist_reordered = gen_dist[pop_order, pop_order]
pheatmap(gen_dist_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Nei's D",
         silent=FALSE) # doesn't match fst and entropy as well
# save
png(str_glue("het_analyses_etc/neisd_heatmap_{version}.png"))
pheatmap(gen_dist_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Nei's D",
         silent=FALSE)
dev.off()

#### Mantel Test ####
# pull in coords
coords_df = read.csv("environ/PHHA_latlong.csv")
coords_df = coords_df[order(coords_df$Two.letter.abbrev),]
coords_df = coords_df[!coords_df$Two.letter.abbrev %in% c("DA", "PM"),]
coords_df$Two.letter.abbrev == rownames(genpop@tab)
coords = coords_df[,c("Long", "Lat")]
geo_dist = dist(coords)
geo_dist = as.data.frame(as.matrix(geo_dist))
colnames(geo_dist) = pop_order
rownames(geo_dist) = pop_order

mantel_pwfst = mantel(pwfst, geo_dist, method = "pearson")
print(mantel_pwfst)

mantel_neisd = mantel(gen_dist, geo_dist, method = "pearson")
print(mantel_neisd)

#### Plot pairwise fst by geo dist ####
pwfst_vec = as.vector(pwfst[lower.tri(pwfst)])
geo_vec = as.vector(geo_dist[lower.tri(geo_dist)])
# plot results
plot(geo_vec, pwfst_vec, xlab = "Geographic distance", 
     ylab = "Fst", pch = 16, cex = 0.9)
abline(lm(pwfst_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_pwfst$statistic, 3),
                       "\nP = ", mantel_pwfst$signif), bty = "n")
# highlight SL on plot
pairs_mat = t(combn(unique(pop), 2))
colors = ifelse(pairs_mat[,1] == "SL" | pairs_mat[,2] == "SL", "red", "black")
plot(geo_vec, pwfst_vec, xlab = "Geographic distance", 
     ylab = "Fst", pch = 16, cex = 0.9, col = colors)
abline(lm(pwfst_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_pwfst$statistic, 3),
                       "\nP = ", mantel_pwfst$signif,
                       "\nRed = SL population"), bty = "n")
# save
png(str_glue("het_analyses_etc/ibd_fst_{version}.png"))
plot(geo_vec, pwfst_vec, xlab = "Geographic distance", 
     ylab = "Fst", pch = 16, cex = 0.9)
abline(lm(pwfst_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_pwfst$statistic, 3),
                       "\nP = ", mantel_pwfst$signif), bty = "n")
dev.off()

#### Plot Nei's d by geo dist ####
neisd_vec = as.vector(gen_dist[lower.tri(gen_dist)])
# plot results
plot(geo_vec, neisd_vec, xlab = "Geographic distance", 
     ylab = "Nei's D", pch = 16, cex = 0.9)
abline(lm(neisd_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_neisd$statistic, 3),
                       "\nP = ", mantel_neisd$signif), bty = "n")
# highlight SL on plot
plot(geo_vec, neisd_vec, xlab = "Geographic distance", 
     ylab = "Nei's D", pch = 16, cex = 0.9, col = colors)
abline(lm(neisd_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_neisd$statistic, 3),
                       "\nP = ", mantel_neisd$signif,
                       "\nRed = SL population"), bty = "n")
# save
png(str_glue("het_analyses_etc/ibd_neisd_{version}.png"))
plot(geo_vec, neisd_vec, xlab = "Geographic distance", 
     ylab = "Nei's D", pch = 16, cex = 0.9)
abline(lm(neisd_vec ~ geo_vec), lwd = 2)
legend("topleft",
       legend = paste0("Mantel r = ", round(mantel_neisd$statistic, 3),
                       "\nP = ", mantel_neisd$signif), bty = "n")
dev.off()

#### map divergent pops - fst ####
coords_df
coords_sf = st_as_sf(coords_df, coords = c("Long", "Lat"), crs = 4326)
fst_df = as.data.frame(pairs_mat)
fst_df$fst = pwfst_vec
colnames(fst_df) = c('pop1', 'pop2', 'fst')
coords_clean = coords_df %>%
  rename(
    pop = Two.letter.abbrev,
    lat = Lat,
    lon = Long
  ) %>%
  select(pop, lat, lon)
fst_map = fst_df %>%
  left_join(coords_clean, by = c("pop1" = "pop")) %>%
  rename(lat1 = lat, lon1 = lon) %>%
  left_join(coords_clean, by = c("pop2" = "pop")) %>%
  rename(lat2 = lat, lon2 = lon)
threshold = quantile(fst_map$fst, 0.9, na.rm = TRUE)
fst_map_sig = fst_map %>%
  filter(fst >= threshold)
lon_range = range(coords_clean$lon)
lat_range = range(coords_clean$lat)
lon_buffer = 0.5
lat_buffer = 0.5
# map
bbox = c(
  left = lon_range[1] - lon_buffer,
  bottom = lat_range[1] - lat_buffer,
  right = lon_range[2] + lon_buffer,
  top = lat_range[2] + lat_buffer
)
# ggmap::register_stadiamaps(key = '')
base_map = get_stadiamap(
  bbox = bbox,
  maptype = "stamen_terrain",   # natural terrain background
  zoom = 7
)
coords_clean2 = coords_clean
coords_clean2$lat = coords_clean2$lat + 0.2
coords_clean2$lon = coords_clean2$lon + 0.2
ggmap(base_map) +
  # Draw FST lines for highly divergent pairs
  geom_segment(data = fst_map_sig,
               aes(x = lon1, y = lat1,
                   xend = lon2, yend = lat2,
                   size = fst),
               color = "red",
               alpha = 0.5) +
  # Population points
  geom_point(data = coords_clean,
             aes(x = lon, y = lat),
             size = 3, color = "black") +
  # Optional labels
  geom_text(data = coords_clean2, aes(x = lon, y = lat, label = pop), nudge_y = 0.1) +
  scale_size(range = c(0.5, 3)) +
  ggtitle("Top 10% FST between Populations") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/fst_outliers_{version}.png"))

#### map divergent pops - neis d ####
neis_df = as.data.frame(pairs_mat)
neis_df$neis = neisd_vec
colnames(neis_df) = c('pop1', 'pop2', 'neis')
neis_map = neis_df %>%
  left_join(coords_clean, by = c("pop1" = "pop")) %>%
  rename(lat1 = lat, lon1 = lon) %>%
  left_join(coords_clean, by = c("pop2" = "pop")) %>%
  rename(lat2 = lat, lon2 = lon)
threshold = quantile(neis_map$neis, 0.9, na.rm = TRUE)
neis_map_sig = neis_map %>%
  filter(neis >= threshold)
ggmap(base_map) +
  # Draw FST lines for highly divergent pairs
  geom_segment(data = neis_map_sig,
               aes(x = lon1, y = lat1,
                   xend = lon2, yend = lat2,
                   size = neis),
               color = "red",
               alpha = 0.5) +
  # Population points
  geom_point(data = coords_clean,
             aes(x = lon, y = lat),
             size = 3, color = "black") +
  # Optional labels
  geom_text(data = coords_clean2, aes(x = lon, y = lat, label = pop), nudge_y = 0.1) +
  scale_size(range = c(0.5, 3)) +
  ggtitle("Top 10% Nei's D between Populations") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/neisd_outliers_{version}.png"))

#### HWE ####
hw_results = hw.test(genind, B = 0)
plot(density(hw_results[,3], na.rm = TRUE), type = "n", xlim = c(0,1),
     main = str_glue("HWE p-values - {version}"),
     ylim = range(density(hw_results[,3], na.rm = TRUE)$y),
     xlab = "Value", ylab = "Density")
polygon(density(hw_results[,3], na.rm = TRUE),
        col = rgb(1, 0, 0, 0.5), border = "red")
# save csv and plot
write.csv(hw_results,
          str_glue("het_analyses_etc/hwe_results_{version}.csv"),
          row.names = FALSE)
png(str_glue("het_analyses_etc/hew_results_{version}.png"))
plot(density(hw_results[,3], na.rm = TRUE), type = "n", xlim = c(0,1),
     main = str_glue("HWE p-values - {version}"),
     ylim = range(density(hw_results[,3], na.rm = TRUE)$y),
     xlab = "Value", ylab = "Density")
polygon(density(hw_results[,3], na.rm = TRUE),
        col = rgb(1, 0, 0, 0.5), border = "red")
dev.off()

#### FIS ####
fis_results = stats$Fis  # hierfstat calculates Fis
print(fis_results)
fis_col_means = colMeans(fis_results, na.rm = TRUE)
print(fis_col_means)
fis_long = fis_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "locus") %>%
  pivot_longer(-locus, names_to = "pop", values_to = "FIS") %>%
  filter(!is.na(FIS) & !is.nan(FIS))
ggplot(fis_long, aes(x = FIS)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  xlab(expression(F[IS])) +
  ylab("Density") +
  ggtitle("Overall FIS Distribution") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/fis_distr_{version}.png"))
ggplot(fis_long, aes(x = FIS, fill = pop)) +
  geom_density(alpha = 0.4) +
  xlab(expression(F[IS])) +
  ylab("Density") +
  ggtitle("FIS Distribution by Population") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(str_glue("het_analyses_etc/fis_by_pop_stacked_{version}.png"))
ggplot(fis_long, aes(x = FIS)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~pop, scales = "free") +
  xlab(expression(F[IS])) +
  ylab("Density") +
  ggtitle("FIS Distribution by Population") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/fis_by_pop_facet_{version}.png"))

#### allelic richness ####
allelic_richness = allelic.richness(hfs)
allelic_richness$Ar
# save
write.csv(allelic_richness$Ar,
          str_glue("het_analyses_etc/allelic_richness_{version}.csv"),
          row.names = FALSE)
ar_means = colMeans(allelic_richness$Ar, na.rm = TRUE)
# plot
ar_long = allelic_richness$Ar %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "locus") %>%
  pivot_longer(-locus, names_to = "pop", values_to = "Ar") %>%
  filter(!is.na(Ar) & !is.nan(Ar))
# overall
ggplot(ar_long, aes(x = Ar)) +
  geom_density(fill = "forestgreen", alpha = 0.5) +
  xlab("Allelic richness (Ar)") +
  ylab("Density") +
  ggtitle("Overall Allelic Richness Distribution") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/allelic_richness_{version}.png"))
ggplot(ar_long, aes(x = Ar, fill = pop)) +
  geom_density(alpha = 0.4) +
  xlab("Allelic richness (Ar)") +
  ylab("Density") +
  ggtitle("Allelic Richness Distribution by Population") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(str_glue("het_analyses_etc/allelic_richness_stacked_{version}.png"))
ggplot(ar_long, aes(x = Ar)) +
  geom_density(fill = "forestgreen", alpha = 0.5) +
  facet_wrap(~pop, scales = "free") +
  xlab("Allelic richness (Ar)") +
  ylab("Density") +
  ggtitle("Allelic Richness Distribution by Population") +
  theme_minimal()
ggsave(str_glue("het_analyses_etc/allelic_richness_facet_{version}.png"))

#### private alleles ####
private_alleles = private_alleles(genind)
print(private_alleles) # none

#### AMOVA ####
#strata(genind) = data.frame(Pop = pop(genind))
#amova_res = poppr.amova(genind, ~Pop)
#print(amova_res)
