library(stringr)
library(dplyr)
par_defaults = par(no.readonly = TRUE)
# take files named PlottingData_k2.csv, PlottingData_k3.csv, etc and plot 
# bars and triangle

max_k = 7
kn_data = list()
for (i in 2:max_k){
  df = read.csv(str_glue("PlottingData_k{i}.csv"))
  kn_data = append(kn_data,list(df))
}
admix_palette_full = c(
  "#1b9e77",  # teal
  "#d95f02",  # orange
  "#7570b3",  # purple
  "#e7298a",  # pink
  "#a6761d",  # brown
  "#1f78b4",  # blue
  "#b2df8a"   # light green
)


plot_one_df = function(df, draw_labels = TRUE, pop_order = NULL){
  # if pop_order is provided, reorder df so that the pops plot in that order
  # along the x axis
  if (!is.null(pop_order)) {
    if (!"population" %in% colnames(df)) {
      stop("df must contain a 'population' column when using pop_order")
    }
    # enforce order; populations not in pop_order go last
    df$population = factor(df$population, levels = pop_order)
    df = df[order(df$population), , drop = FALSE]
  }
  
  n = nrow(df)
  # pull out needed columns
  q_cols  = grep("^q_\\d+$", colnames(df), value = TRUE)
  lb_cols = paste0("lb_", q_cols)
  ub_cols = paste0("ub_", q_cols)
  
  k = length(q_cols)
  pop_cols = admix_palette[seq_len(k)]
  
  # get proportions and error
  admix_mat = t(as.matrix(df[, q_cols]))
  error_mat = as.matrix(df[, ub_cols] - df[, lb_cols])
  error_mat = t(error_mat)  # rows = populations, cols = individuals
  
  # normalize error
  err_range = range(error_mat, na.rm = TRUE)
  err_scaled = (error_mat - err_range[1]) / diff(err_range)
  err_scaled = pmin(pmax(err_scaled, 0), 1)
  
  # create color palette for error
  err_cols = colorRampPalette(c("white", "red"))(100)
  total_error_height = k*0.05
  
  # make room for error cells below bar plot
  if (draw_labels){
    par(mar = c(4, 4, 6, 1), xpd = TRUE)
  } else {
    par(mar = c(1, 4, 6, 1), xpd = TRUE)
  }
  
  
  bp = barplot(
    admix_mat,
    col = pop_cols,
    border = NA,
    space = 0,
    xaxt = "n",
    ylim = c(-total_error_height, 1),
    ylab = "Ancestry proportion"
  )
  
  # put the error cells below
  square_height = 0.04
  
  # for each pop
  for (p in seq_len(k)) {
    # set the y range of where to plot the cell
    y_bottom = -p * square_height
    y_top = y_bottom + square_height
    
    # for each individual
    for (i in seq_len(n)) {
      # pull which error color bin it should be in
      col_idx = ceiling(err_scaled[p, i] * 99) + 1
      # plot the error cell
      rect(
        xleft   = bp[i] - 0.45,
        xright  = bp[i] + 0.45,
        ybottom = y_bottom,
        ytop    = y_top,
        col     = err_cols[col_idx],
        border  = NA
      )
    }
  }
  
  # add population (locality) labels below individuals
  # find contiguous blocks on indivs from same pop/loc
  if (draw_labels){
    loc = as.character(df$population)
    r = rle(loc)
    # start and end indices for each locality block
    ends = cumsum(r$lengths)
    starts = ends - r$lengths + 1
    # find where each pop/loc label should go
    label_x = (bp[starts] + bp[ends]) / 2
    label_y = -total_error_height - 0.02
    text(
      x = label_x,
      y = label_y,
      labels = r$values,
      xpd = TRUE,
      cex = 0.8
    )
    # bracket off which indiv came from that pop/loc
    bracket_y = label_y + 0.01
    segments(
      x0 = bp[starts],
      x1 = bp[ends],
      y0 = bracket_y,
      y1 = bracket_y,
      xpd = TRUE
    )
    segments(
      x0 = bp[starts],
      x1 = bp[starts],
      y0 = bracket_y,
      y1 = bracket_y + 0.01,
      xpd = TRUE
    )
    segments(
      x0 = bp[ends],
      x1 = bp[ends],
      y0 = bracket_y,
      y1 = bracket_y + 0.01,
      xpd = TRUE
    )
  }
  
  # add legend for poulations (admixture)
  legend(
    "top",
    inset = c(0, -0.15),
    legend = q_cols,
    fill = pop_cols,
    horiz = TRUE,
    bty = "n",
    cex = 0.9
  )
  
  
  # add colorramp for error legend
  # ramp parameters
  if (draw_labels){
    ramp_n = 50
    ramp_cols = colorRampPalette(c("white", "red"))(ramp_n)
    
    xleft  = par("usr")[1]
    xright = par("usr")[2]
    ybot = 1.05
    ytop = 1.10
    
    xs = seq(xleft, xright, length.out = ramp_n + 1)
    
    for (i in seq_len(ramp_n)) {
      rect(
        xs[i], ybot,
        xs[i + 1], ytop,
        col = ramp_cols[i],
        border = NA
      )
    }
    text(xleft,  ybot - 0.01, "Low error",  adj = c(0, 1), cex = 0.8)
    text(xright, ybot - 0.01, "High error", adj = c(1, 1), cex = 0.8)
  }
  
}

# plot pop averages for better visibility
plot_pop_means = function(df, pop_order = NULL){
  
  if (!"population" %in% colnames(df)) {
    stop("df must contain a 'population' column")
  }
  
  # admixture columns
  q_cols = grep("^q_\\d+$", colnames(df), value = TRUE)
  k = length(q_cols)
  pop_cols = admix_palette[seq_len(k)]
  
  ## ---- enforce population order ----
  if (!is.null(pop_order)) {
    df$population = factor(df$population, levels = pop_order)
  } else {
    df$population = factor(df$population)
  }
  
  # compute mean admixture per population
  mean_df = aggregate(
    df[, q_cols],
    by = list(population = df$population),
    FUN = mean,
    na.rm = TRUE
  )
  
  # drop unused populations (e.g. not in pop_order)
  mean_df = mean_df[!is.na(mean_df$population), ]
  mean_df$population = droplevels(mean_df$population)
  
  # matrix for barplot
  admix_mat = t(as.matrix(mean_df[, q_cols]))
  
  # margins
  par(mar = c(6, 4, 4, 1))
  
  bp = barplot(
    admix_mat,
    col = pop_cols,
    border = NA,
    space = 0.2,
    ylim = c(0, 1),
    ylab = "Mean ancestry proportion",
    xaxt = "n"
  )
  
  # population labels
  axis(
    side = 1,
    at = bp,
    labels = mean_df$population,
    las = 2,
    tick = FALSE
  )
  
  # admixture legend
  legend(
    "top",
    inset = c(0, -0.15),
    legend = q_cols,
    fill = pop_cols,
    horiz = TRUE,
    bty = "n",
    cex = 0.9
  )
}
# from prev runs
my_order = c("GB", "SL", "PR", "AS", "TB",
             "BC", "WF", "PT", "SC", "OF", "BB",
             "WC", "PC", "FV", "SB", "MC", "SH", "FC",
             "WH", "CA", "BA", "CR", "KC", "ML")
# find order by decreasing q =1
df = kn_data[1][[1]]
q_cols = grep("^q_\\d+$", colnames(df), value = TRUE)
mean_df = aggregate(
  df[, q_cols],
  by = list(population = df$population),
  FUN = mean,
  na.rm = TRUE
)
# view mean_df by q1 descending for k = 2
my_order = c('ML', 'CA', 'CR', 'KC',
             'WH', 'FC', 'BA', 'SH', # very teal ends, orange begins
             'MC', 'GB', 'WF', 'WC',
             'SB', 'PC', 'BB', 'OF',
             'PT', 'BC', 'FV', 'TB',
             'SC', 'SL', 'AS', 'PR')

for (i in 2:max_k){
  admix_palette = admix_palette_full
  df = kn_data[i-1][[1]]
  plot_one_df(df, pop_order = my_order)
}
for (i in 2:max_k){
  admix_palette = admix_palette_full
  df = kn_data[i-1][[1]]
  plot_pop_means(df, pop_order = my_order)
}
# now adjust to also consider k = 3
# we want color palette to be teal, purple, orange (purple in middle of extremes)
admix_palette_list = list(
  c("#1b9e77",
    "#d95f02"), # k = 2
  c("#1b9e77",
    "#7570b3",
    "#d95f02"), # k = 3
  admix_palette_full, # k = 4
  admix_palette_full, # k = 5
  admix_palette_full, # k = 6
  admix_palette_full # k = 7
)
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  #print(k)
  #print(admix_palette)
  df = kn_data[i-1][[1]]
  plot_one_df(df, pop_order = my_order)
}
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  df = kn_data[i-1][[1]]
  plot_pop_means(df, pop_order = my_order)
}
# adjust due to q=3, decreasing q=3
df = kn_data[2][[1]]
q_cols = grep("^q_\\d+$", colnames(df), value = TRUE)
mean_df = aggregate(
  df[, q_cols],
  by = list(population = df$population),
  FUN = mean,
  na.rm = TRUE
)
my_order = c('ML', 'CA', 'CR', 'KC',
             'WH', 'FC', 'BA', # very teal ends, purple begins
             'SH',
             'WC', 'SB', 'PC', 'BB',
             'OF', 'FV', 'SC', 'MC', # purple ends, very orange begins
             'PT', 'GB', 'WF', 'BC',
             'TB', 'AS', 'SL', 'PR')
admix_palette_list = list(
  c("#1b9e77",
    "#d95f02"), # k = 2
  c("#1b9e77",
    "#7570b3",
    "#d95f02"), # k = 3
  c("#1b9e77",
    "#e7298a",
    "#d95f02",
    "#7570b3"
    ), # k = 4
  admix_palette_full, # k = 5
  admix_palette_full, # k = 6
  admix_palette_full # k = 7
)
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  #print(k)
  #print(admix_palette)
  df = kn_data[i-1][[1]]
  plot_one_df(df, pop_order = my_order)
}
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  df = kn_data[i-1][[1]]
  plot_pop_means(df, pop_order = my_order)
}
# adjust due to q=4, decreasing q=2
df = kn_data[3][[1]]
q_cols = grep("^q_\\d+$", colnames(df), value = TRUE)
mean_df = aggregate(
  df[, q_cols],
  by = list(population = df$population),
  FUN = mean,
  na.rm = TRUE
)
my_order = c('ML', 'CR', 'KC', 'BA', # very teal ends, purple begins
             'WH', 'FC', 'CA', 'SH', 'MC',
             'WC', 'SB', 'PC', 'BB',
             'OF', 'FV', 'SC', # purple ends, very orange begins
             'PT', 'GB', 'WF', 'BC',
             'TB', 'AS', 'SL', 'PR')
admix_palette_list = list(
  c("#1b9e77",
    "#d95f02"), # k = 2
  c("#1b9e77",
    "#7570b3",
    "#d95f02"), # k = 3
  c("#1b9e77",
    "#e7298a",
    "#d95f02",
    "#7570b3"
  ), # k = 4
  admix_palette_full, # k = 5
  admix_palette_full, # k = 6
  admix_palette_full # k = 7
)
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  #print(k)
  #print(admix_palette)
  df = kn_data[i-1][[1]]
  plot_one_df(df, pop_order = my_order)
}
for (i in 2:max_k){
  admix_palette = admix_palette_list[i-1][[1]]
  df = kn_data[i-1][[1]]
  plot_pop_means(df, pop_order = my_order)
}

# keep to k = 4
# higher k just splits ancestry too much

# instead of giving a different palette each time
# update q ordering
sort(colnames(kn_data[3][[1]]))

# chatgpt function
relabel_q_columns = function(df, permutation){
  # permutation: named integer vector
  # names = new q index
  # values = old q index
  #
  # example for swapping 1<->2 and 3<->4:
  # c("1"=2, "2"=1, "3"=4, "4"=3)
  
  for (new_i in names(permutation)) {
    old_i = permutation[[new_i]]
    
    # temporary names to avoid overwriting
    tmp_q  = paste0(".__tmp_q_",  new_i)
    tmp_lb = paste0(".__tmp_lb_", new_i)
    tmp_ub = paste0(".__tmp_ub_", new_i)
    
    df[[tmp_q]]  = df[[paste0("q_", old_i)]]
    df[[tmp_lb]] = df[[paste0("lb_q_", old_i)]]
    df[[tmp_ub]] = df[[paste0("ub_q_", old_i)]]
  }
  
  # overwrite originals
  for (new_i in names(permutation)) {
    df[[paste0("q_", new_i)]]     = df[[paste0(".__tmp_q_",  new_i)]]
    df[[paste0("lb_q_", new_i)]]  = df[[paste0(".__tmp_lb_", new_i)]]
    df[[paste0("ub_q_", new_i)]]  = df[[paste0(".__tmp_ub_", new_i)]]
  }
  
  # remove temp columns
  df[, grep("^\\.__tmp_", colnames(df))] = NULL
  
  df
}
perm_k4 = c(
  "1" = 1,
  "2" = 3,
  "3" = 4,
  "4" = 2
)
kn_data[[3]] = relabel_q_columns(kn_data[[3]], perm_k4)
perm_k3 = c(
  "1" = 1,
  "2" = 3,
  "3" = 2
)
kn_data[[2]] = relabel_q_columns(kn_data[[2]], perm_k3)

admix_palette_k_constant = c("#1b9e77",
                             "#d95f02",
                             "#7570b3",
                             "#e7298a"
)

# now plot with one palette, only k 2-3
for (i in 2:4){
  admix_palette = admix_palette_k_constant
  df = kn_data[i-1][[1]]
  plot_pop_means(df, pop_order = my_order)
}
for (i in 2:4){
  admix_palette = admix_palette_k_constant
  df = kn_data[i - 1][[1]]
  plot_one_df(df, pop_order = my_order)
}

for (i in 2:4){
  admix_palette = admix_palette_k_constant
  df = kn_data[i-1][[1]]
  png(
    filename = str_glue("plots/Admixture_k{i}_Indivs.png"),
    width = 2000,
    height = 1600,
    res = 200
  )
  plot_one_df(df, pop_order = my_order)
  dev.off()
}
# plot each k result in separate plot (save)
for (i in 2:4){
  admix_palette = admix_palette_k_constant
  df = kn_data[i-1][[1]]
  png(
    filename = str_glue("plots/Admixture_k{i}_Pops.png"),
    width = 2000,
    height = 1600,
    res = 200
  )
  plot_pop_means(df, pop_order = my_order)
  dev.off()
}
