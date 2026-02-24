library(stringr)

# set k for each run through
k = 7
q_file = str_glue("estpost_files/q{k}.txt")
qhet_file = str_glue("estpost_files/Qhet{k}.txt")
q_admix = read.table(q_file, header = TRUE, sep = ',')
q_het = read.table(qhet_file, header = TRUE, sep = ',')

# pull in indiv file
inds_df = read.table('PHHA.bithin.q30.i85.maxdp12.a70.maf03.012.indv')
# split pop and id numbers
inds_df$individual_id = str_extract(inds_df$V1, "(?<=_)\\d+(?=\\.sorted\\.bam)")
inds_df$population = str_extract(inds_df$V1, "(?<=_)\\w+(?=_)")
n_inds = nrow(inds_df)

# pull in the q and Qhet values for each of k populations and (k 2) population pairs
# for each population
for (i in 1:k){
  # q_admix files 1:n_inds rows are for pop1 the next n_inds rows are pop2.. etc
  q_admix_k_rows = q_admix[((i-1)*n_inds+1):(i*n_inds),]
  colname = str_glue("q_{i}")
  inds_df[[colname]] = q_admix_k_rows$mean
  colname = str_glue("lb_q_{i}")
  inds_df[[colname]] = q_admix_k_rows$ci_0.950_LB
  colname = str_glue("ub_q_{i}")
  inds_df[[colname]] = q_admix_k_rows$ci_0.950_UB
}

# example ordering for 4 pops
# 1-1
# 2-1
# 2-2
# 3-1
# 3-2
# 3-3
# 4-1
# 4-2
# 4-3
# 4-4
ctr = 0
for (p1 in 1:k){
  for (p2 in 1:k){
    if (p2 <= p1){
      ctr = ctr+1
      q_het_k1_k2_rows = q_het[((ctr-1)*n_inds+1):(ctr*n_inds),]
      colname = str_glue("Q_{p1}_{p2}")
      inds_df[[colname]] = q_het_k1_k2_rows$mean
      colname = str_glue("lb_Q_{p1}_{p2}")
      inds_df[[colname]] = q_het_k1_k2_rows$ci_0.950_LB
      colname = str_glue("ub_Q_{p1}_{p2}")
      inds_df[[colname]] = q_het_k1_k2_rows$ci_0.950_UB
    }
  }
}
colnames(inds_df)
# inds_df now has all you need for plotting -> save and repeat for next k
write.csv(inds_df, str_glue("PlottingData_k{k}.csv"), row.names = FALSE)

