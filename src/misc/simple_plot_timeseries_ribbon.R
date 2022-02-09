### TODO: convert into a gneric function, i.e. not specific to our dataset
### time series plotting functions for spatially-explicit stochastic suppression gene drive simulation
#' @export
q_W_shredder = function(q, c=1){(8*q - 3*c*q - c*q^2) / (8 - 4*c)}
q_X_shredder = function(q, c=1){(2*q) / (2 - c + c*q)}
N = function(N, R_max, N_star, c, q, area=50*50){
  n = N/area
  P_female = (2 - c - c*q) / (4 - 2*c)
  b = R_max / ( 1 + (n * ((R_max-2) / (2*N_star))) )
  return(N * P_female * b)
}
simple_plot_timeseries_ribbon = function(R_max, DF_ALL, c=1, area=2500, vec_t=1:50, q_intro=0.01, vec_ribbon_colours=c(rgb(127/256,205/256,187/256,alpha=0.8), rgb(44/256,127/256,184/256,alpha=0.8)), vec_line_type=1:2, show.legend=TRUE){
  ### test #############
  ### fixed parameters
  # area = 50 * 50
  # c = 1
  # vec_t = 1:50
  # q_intro = 0.01
  # sigma = 10
  # vec_ribbon_colours = c(rgb(127/256,205/256,187/256,alpha=0.8), rgb(44/256,127/256,184/256,alpha=0.8))
  # vec_line_type = c(1, 2)
  # ### variable
  # R_max = 5
  #######################
  #################
  ### simulated ###
  idx = (DF_ALL$Rmax==R_max)
  subdat = DF_ALL[idx, ]
  df_n_min = aggregate(n ~ gen + drive_type, data=subdat, FUN=min, na.rm=TRUE)
  df_n_max = aggregate(n ~ gen + drive_type, data=subdat, FUN=max, na.rm=TRUE)
  df_q_min = aggregate(q ~ gen + drive_type, data=subdat, FUN=min, na.rm=TRUE)
  df_q_max = aggregate(q ~ gen + drive_type, data=subdat, FUN=max, na.rm=TRUE)
  #####################
  ### deterministic ###
  vec_q_W_shredder = vec_q_X_shredder = c(q_intro)
  ### population density
  N_star = mean(subdat$n[subdat$gen==1], na.rm=TRUE)/area ### assumming we are at N_star at g=0 to account for the effect of the absorbing boundaries
  vec_N_W_shredder = vec_N_X_shredder = c(N_star * area)
  for (i in vec_t){
    vec_q_W_shredder = c(vec_q_W_shredder, q_W_shredder(q=tail(vec_q_W_shredder,1), c=c))
    vec_q_X_shredder = c(vec_q_X_shredder, q_X_shredder(q=tail(vec_q_X_shredder,1), c=c))
    vec_N_W_shredder = c(vec_N_W_shredder, N(N=tail(vec_N_W_shredder,1), R_max=R_max, N_star=N_star, c=c, q=tail(vec_q_W_shredder,1), area=area))
    vec_N_X_shredder = c(vec_N_X_shredder, N(N=tail(vec_N_X_shredder,1), R_max=R_max, N_star=N_star, c=c, q=tail(vec_q_X_shredder,1), area=area))
  }
  ############
  ### Plot ###
  vec_drive_types = c("W_shredder", "X_shredder")
  vec_drive_types_labels = c("W-shredder", "X-shredder")
  ### population density
  xlim = c(0, max(vec_t))
  ylim = c(0, N_star+max(df_n_max$n))
  plot(0, type="n", xlim=xlim, ylim=ylim, xlab="Generation", ylab="Population size", main=paste0("R = ", R_max))
  for (i in 1:length(vec_drive_types)){
    # i = 2
    drive_type = vec_drive_types[i]
    ribbon_colour = vec_ribbon_colours[i]
    line_type = vec_line_type[i]
    df_min_sub = df_n_min[df_n_min$drive_type==drive_type, ]
    df_max_sub = df_n_max[df_n_max$drive_type==drive_type, ]
    vec_x = c(df_min_sub$gen, rev(df_min_sub$gen))
    vec_y = c(df_min_sub$n, rev(df_max_sub$n))
    polygon(x=vec_x, y=vec_y, col=ribbon_colour, border=NA)
    lines(x=vec_t, y=head(eval(parse(text=paste0("vec_N_", drive_type))), -1), lty=line_type, lwd=2)
  }
  grid()
  if (show.legend){
    legend("topright", legend=c(paste0(vec_drive_types_labels, ": simulated range"), paste0(vec_drive_types_labels, ": deterministic equation")),
          col=c(vec_ribbon_colours, "black", "black"),
          lty=c(1,1,1,2),
          lwd=c(10,10,2,2))
  }
  ### suppression gene drive allele frequency
  xlim = c(0, max(vec_t))
  ylim = c(0, 1)
  plot(0, type="n", xlim=xlim, ylim=ylim, xlab="Generation", ylab="Suppression gene drive allele frequency", main=paste0("R = ", R_max))
  for (i in 1:length(vec_drive_types)){
    # i = 2
    drive_type = vec_drive_types[i]
    ribbon_colour = vec_ribbon_colours[i]
    line_type = vec_line_type[i]
    df_min_sub = df_q_min[df_q_min$drive_type==drive_type, ]
    df_max_sub = df_q_max[df_q_max$drive_type==drive_type, ]
    vec_x = c(df_min_sub$gen, rev(df_min_sub$gen))
    vec_y = c(df_min_sub$q, rev(df_max_sub$q))
    polygon(x=vec_x, y=vec_y, col=ribbon_colour, border=NA)
    lines(x=vec_t, y=head(eval(parse(text=paste0("vec_q_", drive_type))), -1), lty=line_type, lwd=2)
  }
  grid()
  if (show.legend){
    legend("bottomright", legend=c(paste0(vec_drive_types_labels, ": simulated range"), paste0(vec_drive_types_labels, ": deterministic equation")),
          col=c(vec_ribbon_colours, "black", "black"),
          lty=c(1,1,1,2),
          lwd=c(10,10,2,2))
  }
}