# libraries
library(stringr)

# params
max_k = 7

# setup
avg_gprob_f = 'estpost_files/gprob_all_combined.txt'
avg_gprob = as.matrix(read.csv(avg_gprob_f)[,-1])

# visualize avg
heatmap(avg_gprob, NA, NA, main = 'Average')
# for each k, visualize that version
for (k in 2:max_k){
  gprob_k_f = str_glue("estpost_files/gprob_k{k}.txt")
  gprob_k = as.matrix(read.csv(gprob_k_f)[,-1])
  heatmap(gprob_k, NA, NA, main = str_glue("k = {k}"))
  # also plot dif between that version and the avg
  heatmap((gprob_k - avg_gprob), NA, NA,
          main = str_glue('Difference from average, k = {k}'),
          col = colorRampPalette(c("blue", "red"))(100))
}
