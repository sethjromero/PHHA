library(stringr)
par_defaults = par(no.readonly = TRUE)
# take files named PlottingData_k2.csv, PlottingData_k3.csv, etc and plot bars and triangle

max_k = 7
kn_data = list()
for (i in 2:max_k){
  df = read.csv(str_glue("PlottingData_k{i}.csv"))
  kn_data = append(kn_data,list(df))
}

admix_palette = c(
  "#1b9e77",  # teal
  "#d95f02",  # orange
  "#7570b3",  # purple
  "#e7298a",  # pink
  "#a6761d",  # brown
  "#1f78b4",  # blue
  "#b2df8a"   # light green
)


plot_one_df = function(df, draw_labels = TRUE){
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
    loc = df$population
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
# plot each k result in separate plot (view)
par(par_defaults)
for (i in 2:max_k){
  df = kn_data[i-1][[1]]
  plot_one_df(df)
}

# plot each k result in separate plot (save)
par(par_defaults)
for (i in 2:max_k){
  df = kn_data[i-1][[1]]
  png(
    filename = str_glue("plots/Admixture_k{i}.png"),
    width = 2000,
    height = 1600,
    res = 200
  )
  plot_one_df(df)
  dev.off()
}



