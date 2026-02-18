# libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(MASS)
library(LEA)
setwd('entropy')

# path to the .012.indv file
indiv_file= 'PHHA.bithin.q30.i85.maxdp12.a70.maf03.012.indv'

# path to the pntest_mean file
pntest_file = 'pntest_mean_PHHA.bithin.q30.i85.maxdp12.a70.maf03.recode.txt'

# define function to pull population ID out of filenames
# XX_YY_NNN -> pop YY sample num NNN
makePopId <- function(fileIndv){
  PopIDdf = read.table(fileIndv, sep="\t") %>%
    as.data.frame() %>%
    rename(All = V1) %>%
    mutate(Population = str_split_i(All, '_', 2),
           ID = str_split_i(All, '_', 3))
  return(PopIDdf)
}

# define function to do pca on geno matrix
# g is ind x loci
# pulls mean dosage of each locus, calc allele freq from that
# normalize / scale all genos by mean and variance of the site
# missing vals -> 0 (mean geno)
# do PCA
PCA_entropy <- function(g){
  colmean = apply(g, 2, mean, na.rm = T)
  normalize = matrix(nrow = nrow(g), ncol = ncol(g))
  af = colmean/2
  for (m in 1:length(af)){
    nr = g[,m]-colmean[m]
    dn = sqrt(af[m]*(1-af[m]))
    normalize[,m] = nr/dn
  }
  normalize[is.na(normalize)] = 0
  method1 = prcomp(normalize, scale. = F,center = F)
  pca_df = method1$x[,1:27]
  return(pca_df)
}

# read indiv and geno files
PopID <- makePopId(indiv_file)
g <- read.table(pntest_file, header = F)


# pull first 10 PCs
# transpose because pntest was loci x indv
# tack on pop IDs
pca_df <- PCA_entropy(t(g)) %>%
  .[,1:10] %>%
  cbind(PopID)

# define function to get initial pop assignments
# kmeans clustering on pca space
# lda to find pcs that split clusters
# use the lda model to get prob of membership in each cluster
writeLDAfile <- function(pcaDF, k){
  kCluster = kmeans(pcaDF[,1:5], k, iter.max = 10, nstart = 10, algorithm = 'Hartigan-Wong')
  ldaOut = lda(x = pcaDF[,1:5], grouping = kCluster$cluster, CV = T)
  write.table(round(ldaOut$posterior, 5),
              file = paste('ldak', as.character(k), '.txt', sep = ''),
              quote = F, row.names = F, col.names = F)
}

# write a file for each num pops you want to run entropy with
for (i in 2:10){
  writeLDAfile(pca_df, i)
}


# make header for .mpgl
# pull YY_NNN back together for indiv ID
# this should be 2nd line of final .mpgl file
PopID_list <- paste(PopID$Pop, PopID$ID, sep = '_')
# need n indiv and n loci as the 1st line
header <- data.frame(dims = NA, PopID_list)
df <- t(header)
dims <- paste(dim(g)[2], dim(g)[1], sep = " ")
df[1,1] <- dims
f = paste('entropy_header.txt', sep = '')

# write 2 line header
write.table(df,f,
            sep = " ", na = "",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
