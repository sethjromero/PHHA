library(rhdf5)

# update filename for each run you want to examine
#f = 'hdf5_files/PHHAentropyout_k2_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k3_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k4_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k5_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k6_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k7_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k8_c1.hdf5'
#f = 'hdf5_files/PHHAentropyout_k9_c1.hdf5'
f = 'hdf5_files/PHHAentropyout_k10_c1.hdf5'


# see items in f
h5ls(f)
# should be Q, alpha, args, deviance, fst, gamma, gprob, q, and zprob

# look at q first
w.q<-h5read(f,"q")
# dim should be n iters x n pops x n indivs
dim(w.q)
# look at 4 random indivs at once
par(mfrow=c(2,2))
# random pop
# random individual
# plot all iters
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='',ylim=c(0,1))
mtext("Trace plots for admixture estimates", outer = TRUE, side=3, line=-2)

# look at Q
w.Q<-h5read(f,"Q")
# dim should be n iters x n pop combos (00 ,01/10, 11 for 2 pops) x n indivs
dim(w.Q)
# look at 4 random indivs at once
par(mfrow=c(2,2))
# random Q (11,12,22)
# random individual
# plot all iters
plot(w.Q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.Q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.Q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.Q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.Q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.Q)[3],1)],type="l",ylab='',ylim=c(0,1))
plot(w.Q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.Q)[3],1)],type="l",ylab='',ylim=c(0,1))
mtext("Trace plots for admixture estimates", outer = TRUE, side=3, line=-2)

# alpha = the genetic var in the ancestral pop
# (prior on pi = allele freqs in ancestral pop)
w.a<-h5read(f,"alpha")
dim(w.a)
par(mfrow=c(1,1))
plot(w.a,type='l')

# deviance of model
w.d<-h5read(f,"deviance")
dim(w.d)
par(mfrow=c(1,1))
plot(w.d,type='l')

# fst of each pop from ancestor
w.fst<-h5read(f,"fst")
# dim is n iters x n pops
dim(w.fst)
par(mfrow=c(dim(w.fst)[2],1))
for (i in 1:dim(w.fst)[2]){
  plot(w.fst[,i],type='l')
}
# look at each separately too
par(mfrow=c(1,1))
for (i in 1:dim(w.fst)[2]){
  plot(w.fst[,i],type='l')
}

# prior membership in each pop
w.g<-h5read(f,"gamma")
dim(w.g)
par(mfrow=c(1,1))
plot(w.g,type='l')

#genotype probabilities
w.G<-h5read(f,"gprob")
# dim is n genotypes x n indivs x n loci
# why not also over iters? i think it just saves the final state
dim(w.G)
par(mfrow=c(2,2))
# prob of ref ref
# indiv1
# random snp
plot(w.G[3,1,sample(1:dim(w.G)[3],1)],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,1,sample(1:dim(w.G)[3],1)],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,1,sample(1:dim(w.G)[3],1)],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,1,sample(1:dim(w.G)[3],1)],ylab='',ylim=c(-0.1,1.1))
# or prop across SNPs
# prob of ref ref
# random indiv
# across SNPS (x)
plot(w.G[3,sample(1:dim(w.G)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,sample(1:dim(w.G)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,sample(1:dim(w.G)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.G[3,sample(1:dim(w.G)[2],1),],ylab='',ylim=c(-0.1,1.1))

#now instead of geno, source of allele
w.z<-h5read(f,"zprob")
# n pops x n indivs x n snps
dim(w.z)
par(mfrow=c(2,2))
#prob of pop0 pop0
#indiv1
#random snp
plot(w.z[1,1,sample(1:dim(w.z)[3],1)],ylab='',ylim=c(0,1))
plot(w.z[1,1,sample(1:dim(w.z)[3],1)],ylab='',ylim=c(0,1))
plot(w.z[1,1,sample(1:dim(w.z)[3],1)],ylab='',ylim=c(0,1))
plot(w.G[3,1,sample(1:dim(w.G)[3],1)],ylab='',ylim=c(0,1))
# or prop across SNPs
# prob of ref ref
# random indiv
# across SNPS (x)
plot(w.z[1,sample(1:dim(w.z)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.z[1,sample(1:dim(w.z)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.z[1,sample(1:dim(w.z)[2],1),],ylab='',ylim=c(-0.1,1.1))
plot(w.z[1,sample(1:dim(w.z)[2],1),],ylab='',ylim=c(-0.1,1.1))

