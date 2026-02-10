#### packages ####
library(stringr)
library(adegenet)
library(hierfstat)
library(pheatmap)
library(vegan) 
library(pegas)
library(poppr)

#### data prep ####
# convert gprob file (&ind file) to genind format
# assume wd = the PHHA dir
gprob_file = "entropy/estpost_files/gprob_all_combined.txt"
g = read.csv(gprob_file)
ind_file = "entropy/PHHA.bithin.q30.i88.maxdp12.a70.maf03.012.indv"
i = read.table(ind_file)
i$id = str_extract(i$V1, "(?<=_)\\d+(?=\\.sorted\\.bam)")
i$pop = str_extract(i$V1, "(?<=_)\\w+(?=_)")
g[,1] = paste0(i$pop, "_", i$id)
pop = i$pop
rm(i, gprob_file, ind_file)
# 245 inds, 6809 snps

#### call gts from gprobs ####
# try 0/1/2 by rounding or by 0-0.3, 0.7-1.3, 1.6-2.0 with all else NA (missing)
g_mat = as.matrix(g[, 2:ncol(g)])

g_mat_round = round(g_mat)
sum(is.na(g_mat_round))
g_mat_strict = ifelse(g_mat >= 0 & g_mat <= 0.3, 0,
                      ifelse(g_mat >= 0.7 & g_mat <= 1.3, 1,
                             ifelse(g_mat >= 1.7 & g_mat <= 2, 2,
                                    NA)))
sum(is.na(g_mat_strict))
# convert dosage to gt
gt_mat_round = ifelse(g_mat_round==0, 11,
                      ifelse(g_mat_round==1, 12,
                             ifelse(g_mat_round==2, 22,
                                    NA)))
sum(is.na(gt_mat_round))
gt_mat_strict = ifelse(g_mat_strict==0, 11,
                      ifelse(g_mat_strict==1, 12,
                             ifelse(g_mat_strict==2, 22,
                                    NA)))
sum(is.na(gt_mat_strict))


#### genind ####
genind_round = df2genind(gt_mat_round, ploidy=2, pop=pop, sep="") 
genind_strict = df2genind(gt_mat_strict, ploidy=2, pop=pop, sep="") 

# summarize
summary(genind_round)
summary(genind_strict)

#### heterozygosity ####
# observed het per locus (all individuals)
hobs_round = summary(genind_round)$Hobs
hobs_strict = summary(genind_strict)$Hobs
# plot comparison between versions
dens1 = density(hobs_round)
dens2 = density(hobs_strict)
png("het_analyses_etc/Ho.png")
plot(dens1, type = "n", xlim = c(0,1), 
     ylim = range(c(dens1$y, dens2$y)), 
     main = "Observed heterozygosity", xlab = "Value", ylab = "Density")
polygon(dens1, col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(dens2, col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
dev.off()

# expected het per locus (all individuals/pops)
hexp_round = summary(genind_round)$Hexp
hexp_strict = summary(genind_strict)$Hexp
# plot comparison between versions
dens1 = density(hexp_round)
dens2 = density(hexp_strict)
png("het_analyses_etc/He.png")
plot(dens1, type = "n", xlim = c(0,1), 
     ylim = range(c(dens1$y, dens2$y)), 
     main = "Expected heterozygosity", xlab = "Value", ylab = "Density")
polygon(dens1, col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(dens2, col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
# higher het than expected at some loci
dev.off()

# compare observed and expected for each version
dens1 = density(hobs_round)
dens2 = density(hexp_round)
png("het_analyses_etc/Ho_He_round.png")
plot(dens1, type = "n", xlim = c(0,1), 
     ylim = range(c(dens1$y, dens2$y)), 
     main = "Observed vs expected heterozygosity - round", xlab = "Value", ylab = "Density")
polygon(dens1, col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(dens2, col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
dev.off()
# compare observed and expected for each version
dens1 = density(hobs_strict)
dens2 = density(hexp_strict)
png("het_analyses_etc/Ho_He_strict.png")
plot(dens1, type = "n", xlim = c(0,1), 
     ylim = range(c(dens1$y, dens2$y)), 
     main = "Observed vs expected heterozygosity - strict", xlab = "Value", ylab = "Density")
polygon(dens1, col = rgb(1, 0, 0, 0.5), border = "red")  # red with 50% transparency
polygon(dens2, col = rgb(0, 0, 1, 0.5), border = "blue") # blue with 50% transparency
dev.off()

# convert to hierfstat format
hfs_round = genind2hierfstat(genind_round)
hfs_strict = genind2hierfstat(genind_strict)

#### fstats ####
# stats
bs_round = basic.stats(hfs_round)
bs_round
df = do.call(rbind, lapply(names(bs_round$pop.freq), function(locus) {
  freqs = bs_round$pop.freq[[locus]][1, ]  # first row = first allele
  data.frame(t(freqs), row.names = locus)
}))
colnames(df) = colnames(bs_round$n.ind.samp)
write.csv(df, "het_analyses_etc/allele_freqs_round.csv")
bs_strict = basic.stats(hfs_strict)
bs_strict
df = do.call(rbind, lapply(names(bs_strict$pop.freq), function(locus) {
  freqs = bs_strict$pop.freq[[locus]][1, ]  # first row = first allele
  data.frame(t(freqs), row.names = locus)
}))
colnames(df) = colnames(bs_strict$n.ind.samp)
write.csv(df, "het_analyses_etc/allele_freqs_strict.csv")

# calculate pairwise FST on dosage matrices
pairwise_fst_round = pairwise.fst.dosage(g_mat_round, pop)
print(pairwise_fst_round)
write.csv(pairwise_fst_round, "het_analyses_etc/pairwise_fst_round.csv")
# compare to overall
print(bs_round$overall[7])

# calculate pairwise FST
pairwise_fst_strict = pairwise.fst.dosage(g_mat_strict, pop)
print(pairwise_fst_strict)
write.csv(pairwise_fst_strict, "het_analyses_etc/pairwise_fst_strict.csv")
print(bs_strict$overall[7])

# plot comparison
png("het_analyses_etc/fst_heatmap_round.png")
pheatmap(pairwise_fst_round, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST", silent=FALSE)
dev.off()
png("het_analyses_etc/fst_heatmap_strict.png")
pheatmap(pairwise_fst_strict, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST", silent=FALSE)
dev.off()

# reorder to compare better to the entropy results
pop_order = c("DA", "PM", "GB", "SL", "PR", "AS", "TB",
              "BC", "WF", "PT", "SC", "OF", "BB",
              "WC", "PC", "FV", "SB", "BA", "MC", "SH", "FC",
              "WH","CA", "CR", "KC", "ML")
pairwise_fst_strict_reordered = pairwise_fst_strict[pop_order, pop_order]
png("het_analyses_etc/fst_heatmap_round.png")
pheatmap(pairwise_fst_strict_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST", silent=FALSE)
dev.off()
pairwise_fst_round_reordered = pairwise_fst_round[pop_order, pop_order]
png("het_analyses_etc/fst_heatmap_strict.png")
pheatmap(pairwise_fst_round_reordered, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Pairwise FST", silent=FALSE)
dev.off()

#### Nei's distanct ####
# genetic vs geo distance
# genetic distance
genpop_round = genind2genpop(genind_round)
rownames(genpop_round@tab) # order of the pops for dist matrix
gen_dist_round = dist.genpop(genpop_round)
genpop_strict = genind2genpop(genind_strict)
rownames(genpop_strict@tab) # order of the pops for dist matrix
gen_dist_strict = dist.genpop(genpop_strict)
# order is the same
rownames(genpop_round@tab) == rownames(genpop_strict@tab) 

# read lat/lon data
coords_df = read.csv("environ/PHHA_latlong.csv")
# reorder 
coords_df = coords_df[order(coords_df$Two.letter.abbrev),]
coords_df$Two.letter.abbrev == rownames(genpop_strict@tab)

# pull just lon and lat
coords = coords_df[,c("Long", "Lat")]
geo_dist = dist(coords)

#### mantel test ####
mantel_res_round = mantel(gen_dist_round, geo_dist, method = "pearson")
print(mantel_res_round)

mantel_res_strict = mantel(gen_dist_strict, geo_dist, method = "pearson")
print(mantel_res_strict)

#### maps of gen dist ####
# plot gen vs geo dist
# convert distance objects / matrices to vectors
gen_vec_round = as.vector(gen_dist_round)
gen_vec_strict = as.vector(gen_dist_strict)
geo_vec = as.vector(geo_dist)

plot(
  geo_vec, gen_vec_round,
  xlab = "Geographic distance",
  ylab = "Genetic distance",
  pch = 16,
  cex = 0.9
)
abline(lm(gen_vec_round ~ geo_vec), lwd = 2)
legend("topright",
       legend = paste0("Mantel r = ", round(mantel_res_round$statistic, 3),
                       "\nP = ", mantel_res_round$signif), bty = "n")

plot(
  geo_vec, gen_vec_strict,
  xlab = "Geographic distance",
  ylab = "Genetic distance",
  pch = 16,
  cex = 0.9
)
abline(lm(gen_vec_strict ~ geo_vec), lwd = 2)
legend("topright",
       legend = paste0("Mantel r = ", round(mantel_res_strict$statistic, 3),
                       "\nP = ", mantel_res_strict$signif), bty = "n")

identical(labels(gen_dist_round), coords_df[labels(geo_dist),]$Two.letter.abbrev)
identical(labels(gen_dist_strict), coords_df[labels(geo_dist),]$Two.letter.abbrev)

# redo scatterplots and color all points that contain DA or PM differenly
d = as.dist(gen_dist_round)
d
gen_vec_round = as.vector(d)
gen_vec_round
pairs_mat = t(combn(unique(pop), 2))

colors = ifelse(pairs_mat[,1] == "DA" | pairs_mat[,2] == "DA", "lightgreen", 
                ifelse(pairs_mat[,1] == "PM" | pairs_mat[,2] == "PM", "thistle", "black"))
colors[pairs_mat[,1]=="DA" & pairs_mat[,2]=="PM"] = "deeppink"

png("het_analyses_etc/mantel_round.png")
plot(
  geo_vec, gen_vec_round,
  xlab = "Geographic distance",
  ylab = "Genetic distance",
  pch = 16,
  cex = 0.9,
  col = colors
)
abline(lm(gen_vec_round ~ geo_vec), lwd = 2)
legend("topright",
       legend = paste0("Mantel r = ", round(mantel_res_round$statistic, 3),
                       "\nP = ", mantel_res_round$signif), bty = "n")
dev.off()


d = as.dist(gen_dist_strict)
d
gen_vec_strict = as.vector(d)
gen_vec_strict
pairs_mat = t(combn(unique(pop), 2))

colors = ifelse(pairs_mat[,1] == "DA" | pairs_mat[,2] == "DA", "lightgreen", 
                ifelse(pairs_mat[,1] == "PM" | pairs_mat[,2] == "PM", "thistle", "black"))
colors[pairs_mat[,1]=="DA" & pairs_mat[,2]=="PM"] = "deeppink"

png("het_analyses_etc/mantel_strict.png")
plot(
  geo_vec, gen_vec_strict,
  xlab = "Geographic distance",
  ylab = "Genetic distance",
  pch = 16,
  cex = 0.9,
  col = colors
)
abline(lm(gen_vec_strict ~ geo_vec), lwd = 2)
legend("topright",
       legend = paste0("Mantel r = ", round(mantel_res_strict$statistic, 3),
                       "\nP = ", mantel_res_strict$signif), bty = "n")
dev.off()


png("het_analyses_etc/gen_dist_map_round.png")
# plot map
plot(NA,
     xlab = "Longitude",
     ylab = "Latitude",
     asp = 1,
     xlim = c(-121, -109),
     ylim = c(39, 50)
)
# connect pops with lines with thickness proportional to fst / gen dist

# coords for each pair
coord1 = coords_df[match(pairs_mat[,1], coords_df$Two.letter.abbrev), c("Long","Lat")]
coord2 = coords_df[match(pairs_mat[,2], coords_df$Two.letter.abbrev), c("Long","Lat")]
# scale gen dist
gen_dist_scaled_round = (gen_vec_round - min(gen_vec_round)) / (max(gen_vec_round) - min(gen_vec_round))
color_fun = colorRampPalette(c("red","blue"))
colors = color_fun(100)[as.numeric(cut(gen_dist_scaled_round, breaks = 100))]
colors_trans = adjustcolor(colors, alpha.f = 0.5)
for(i in seq_along(gen_vec_round)) {
  segments(coord1$Long[i], coord1$Lat[i], coord2$Long[i], coord2$Lat[i],
           col = colors_trans[i], lwd = 2)
}
points(
  coords_df$Long, coords_df$Lat,
  pch = 19,
  col = "black",
  cex = 1.5,
  xlab = "Longitude",
  ylab = "Latitude",
  asp = 1,
  xlim = c(-121, -109),
  ylim = c(39, 50)
)
text(
  coords_df$Long, coords_df$Lat,
  labels = coords_df$Two.letter.abbrev,
  pos = 3,
  cex = 1.5
)
dev.off()

png("het_analyses_etc/gen_dist_map_strict.png")
# plot map
plot(NA,
     xlab = "Longitude",
     ylab = "Latitude",
     asp = 1,
     xlim = c(-121, -109),
     ylim = c(39, 50)
)
# connect pops with lines with thickness proportional to fst / gen dist

# coords for each pair
coord1 = coords_df[match(pairs_mat[,1], coords_df$Two.letter.abbrev), c("Long","Lat")]
coord2 = coords_df[match(pairs_mat[,2], coords_df$Two.letter.abbrev), c("Long","Lat")]
# scale gen dist
gen_dist_scaled_strict = (gen_vec_strict - min(gen_vec_strict)) / (max(gen_vec_strict) - min(gen_vec_strict))
color_fun = colorRampPalette(c("red","blue"))
colors = color_fun(100)[as.numeric(cut(gen_dist_scaled_strict, breaks = 100))]
colors_trans = adjustcolor(colors, alpha.f = 0.5)
for(i in seq_along(gen_vec_strict)) {
  segments(coord1$Long[i], coord1$Lat[i], coord2$Long[i], coord2$Lat[i],
           col = colors_trans[i], lwd = 2)
}
points(
  coords_df$Long, coords_df$Lat,
  pch = 19,
  col = "black",
  cex = 1.5,
  xlab = "Longitude",
  ylab = "Latitude",
  asp = 1,
  xlim = c(-121, -109),
  ylim = c(39, 50)
)
text(
  coords_df$Long, coords_df$Lat,
  labels = coords_df$Two.letter.abbrev,
  pos = 3,
  cex = 1.5
)
dev.off()

#### maps of fst ####
# do with fst instead of nei's dist
fst_vec_strict = as.vector(pairwise_fst_strict[upper.tri(pairwise_fst_strict)])
fst_scaled = (fst_vec_strict - min(fst_vec_strict)) / (max(fst_vec_strict) - min(fst_vec_strict))
color_fun = colorRampPalette(c("blue", "red"))
colors = color_fun(100)[as.numeric(cut(fst_scaled, breaks = 100))]
colors_trans = adjustcolor(colors, alpha.f = 0.5)

png("het_analyses_etc/fst_map_strict.png")
plot(NA,
     xlab = "Longitude",
     ylab = "Latitude",
     asp = 1,
     xlim = c(-121, -109),
     ylim = c(39, 50)
)
for(i in seq_along(fst_vec_strict)) {
  segments(coord1$Long[i], coord1$Lat[i],
           coord2$Long[i], coord2$Lat[i],
           col = colors_trans[i], lwd = 2)
}
points(
  coords_df$Long, coords_df$Lat,
  pch = 19,
  col = "black",
  cex = 1.5,
  xlab = "Longitude",
  ylab = "Latitude",
  asp = 1,
  xlim = c(-121, -109),
  ylim = c(39, 50)
)
text(
  coords_df$Long, coords_df$Lat,
  labels = coords_df$Two.letter.abbrev,
  pos = 3,
  cex = 1.5
)
dev.off()

# do with fst instead of nei's dist
fst_vec_round = as.vector(pairwise_fst_round[upper.tri(pairwise_fst_round)])
fst_scaled = (fst_vec_round - min(fst_vec_round)) / (max(fst_vec_round) - min(fst_vec_round))
color_fun = colorRampPalette(c("blue", "red"))
colors = color_fun(100)[as.numeric(cut(fst_scaled, breaks = 100))]
colors_trans = adjustcolor(colors, alpha.f = 0.5)

png("het_analyses_etc/fst_map_round.png")
plot(NA,
     xlab = "Longitude",
     ylab = "Latitude",
     asp = 1,
     xlim = c(-121, -109),
     ylim = c(39, 50)
)
for(i in seq_along(fst_vec_round)) {
  segments(coord1$Long[i], coord1$Lat[i],
           coord2$Long[i], coord2$Lat[i],
           col = colors_trans[i], lwd = 2)
}
points(
  coords_df$Long, coords_df$Lat,
  pch = 19,
  col = "black",
  cex = 1.5,
  xlab = "Longitude",
  ylab = "Latitude",
  asp = 1,
  xlim = c(-121, -109),
  ylim = c(39, 50)
)
text(
  coords_df$Long, coords_df$Lat,
  labels = coords_df$Two.letter.abbrev,
  pos = 3,
  cex = 1.5
)
dev.off()

#### hwe ####
hw_results = hw.test(genind_round, B = 0)
hist(hw_results[,3])

hw_results = hw.test(genind_strict, B = 0)
hist(hw_results[,3])

#### FIS ####
fis_results = bs_round$Fis  # hierfstat calculates Fis
print(fis_results)
fis_col_means = colMeans(fis_results, na.rm = TRUE)
print(fis_col_means)

fis_results = bs_strict$Fis  # hierfstat calculates Fis
print(fis_results)
fis_col_means = colMeans(fis_results, na.rm = TRUE)
print(fis_col_means)

#### alleleic richness ####
allelic_richness = allelic.richness(hfs_round)
allelic_richness$Ar
ar_means = colMeans(allelic_richness$Ar, na.rm = TRUE)
ar_means
private_alleles = private_alleles(genind_round)
print(private_alleles)

allelic_richness = allelic.richness(hfs_strict)
allelic_richness$Ar
ar_means = colMeans(allelic_richness$Ar, na.rm = TRUE)
ar_means
private_alleles = private_alleles(genind_strict)
print(private_alleles)


#### amova ####
strata(genind_round) = data.frame(Pop = pop(genind_round))
amova_res = poppr.amova(genind_round, ~Pop)
print(amova_res)

strata(genind_strict) = data.frame(Pop = pop(genind_strict))
amova_res = poppr.amova(genind_strict, ~Pop)
print(amova_res)
