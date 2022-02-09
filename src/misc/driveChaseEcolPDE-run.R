###################################################################################################
###                                                                                             ###
### Numerical solutions to 1D PDEs for W-shredder, X-shredder, and TADS suppression gene drives ###
###                                                                                             ###
###################################################################################################

#################################################################
### Load packages for solving the PDEs and parallel execution ###
#################################################################
library(deSolve)
library(doParallel)
library(doSNOW)

#########################################
### Load PDEs and execution functions ###
#########################################
source("driveChaseEcolPDE-functions.R")

##################
### Parameters ###
##################
nx = 1000 ### number of x-axis windows
nt = 200  ### number of generations
dx = 1.00 ### distance between x-axis windows
K = 5 ### maximum population size in each interval, i.e. carrying capacity
threshold = 0.001 ### the minimum population density at which q becomes zero, i.e. the population has essentially crashed
coef_K = 1 ### fraction of the carrying capacity at which we start the WT population at
q_intro = 0.5 ### introduce drive heterozygotes at the middle of the landscape, i.e. all individuals at that point are heterozygote drive-carriers
intro_pop = FALSE
intro_drive = TRUE
vec_x = (c(1:nx)-round(nx/2))*dx ### x-axis values
vec_colours_u = colorRampPalette(c("#a1d99b", "#2c7fb8"))(nt+1)
vec_colours_q = colorRampPalette(c("#feb24c", "#de2d26"))(nt+1)
plot.out = FALSE

#################
### Variables ###
#################
vec_pde   = c("pdq_pdt", "pdn_pdt")
vec_D     = round(10^(seq(0, 4, length=21)))
vec_r     = round(seq(0,K,length=21),2) ### r<= 0.25 and r<=1 also r>1 are interesting based on the theoretical results
vec_c     = seq(0.0,1.0,length=21)
# vec_dType = c("W_shredder", "X_shredder", "TADS")
# ### tests
# vec_pde   = c("pdq_pdt", "pdn_pdt")
# vec_D     = 100
# vec_r     = seq(0, 4,length=81)
# vec_c     = 1
vec_dType = c("W_shredder", "X_shredder")

########################
### Output dataframe ###
########################
# ODE_OUT = array(NA, dim=c(length(vec_pde), length(vec_D), length(vec_r), length(vec_c), length(vec_dType), nt+1, nx))
df = expand.grid(D=vec_D, r=vec_r, c=vec_c, dType=vec_dType, period=0, chase=FALSE, crash=FALSE, q_max=NA, drive_steepness=NA)
### datafram of input variable indices
df_idx = expand.grid(i=1:length(vec_D), j=1:length(vec_r), k=1:length(vec_c), l=1:length(vec_dType))

##################################
### Numerically solve the PDEs ###
##################################
nCores = parallel::detectCores()
cl = makeSOCKcluster(nCores-1) ### set aside 1 free core
registerDoSNOW(cl)
pb = txtProgressBar(min=0, max=length(vec_D)*length(vec_r)*length(vec_c)*length(vec_dType), char="=", style=3); counter=1
progress = function(i) setTxtProgressBar(pb, i)
opts = list(progress=progress)
foreach(idx = 1:nrow(df_idx), .options.snow=opts, .combine="c") %dopar% {
  library(deSolve)
  i = df_idx$i[idx]
  j = df_idx$j[idx]
  k = df_idx$k[idx]
  l = df_idx$l[idx]
  ### test: i=j=k=l=1
  D     = vec_D[i]
  r     = vec_r[j]
  c     = vec_c[k]
  dType = vec_dType[l]
  ### test: D=10; r=3; c=1; dType="W_shredder"; plot.out=TRUE; intro_pop=FALSE; intro_drive=TRUE
  ### numerically solve the PDEs
  ode_out = func_run(nt=nt, nx=nx, dx=dx, K=K, D=D, r=r, dType=dType, c=c,
                    coef_K=coef_K, q_intro=q_intro, intro_pop=intro_pop, intro_drive=intro_drive,
                    threshold=threshold, progress.bar=TRUE)
  ### extract the output arrays
  array_q = ode_out$array_q
  array_n = ode_out$array_n
  ### save output for merging into ODE_OUT array after parallel execution
  fname_out = paste0(i, "-", j, "-",  k, "-",  l, ".rds")
  saveRDS(array_q, file=paste0("q-", fname_out))
  saveRDS(array_n, file=paste0("n-", fname_out))
  ### retrun NULL to save mmemory
  NULL
}
close(pb)
stopCluster(cl)

############################
### Summarise the output ###
############################
### summarise individual numerical solution output into a dataframe and detect chasing
# system("mkdir PDE_output")
plot.show = FALSE
pb = txtProgressBar(min=0, max=nrow(df_idx), char="=", style=3)
for (idx in 1:nrow(df_idx)){
  i = df_idx$i[idx]
  j = df_idx$j[idx]
  k = df_idx$k[idx]
  l = df_idx$l[idx]
  fname_out = paste0(i, "-", j, "-",  k, "-",  l, ".rds")
  array_q = readRDS(paste0("q-", fname_out))
  array_n = readRDS(paste0("n-", fname_out))
  # ODE_OUT[1, i, j, k, l, , ] = array_q
  # ODE_OUT[2, i, j, k, l, , ] = array_n
  ### use the middle of the landscape
  idx_mid = round(ncol(array_n)/2)
  vec_n = array_n[,idx_mid]
  vec_q = array_q[,idx_mid]
  ### remove the first 10 generations since there is an intrinsic osciliation after drive introduction
  vec_n = vec_n[-1:-10]
  ### detect the maximum period
  spec = tryCatch(spectrum(vec_n, log="no", plot=FALSE), error=function(e){0})
  # spec = spectrum(vec_n, log="no", plot=FALSE)
  period = tryCatch(spec$freq[spec$spec==max(spec$spec)][1]*nt, error=function(e){0})
  # period = spec$freq[spec$spec==max(spec$spec)]*nt
  ### visualise
  if (plot.show) {
    plot_spectrum_densities_and_frequencies(array_n, array_q, nt, nx, K, spec, period)
  }
  ### chasing occured if period > 2
  if (period > 2) {
    df$period[idx] = period
    df$chase[idx] = TRUE
  }
  ### the population crashed if the population density across the landscape at the final generation in zero
  if (sum(array_n[nrow(array_n),],na.rm=TRUE) < ncol(array_n)*threshold) {
    df$crash[idx] = TRUE
  }
  ### drive allele wave statistics
  df$q_max[idx] = func_q_wave_stats(array_q=array_q, K=K, nx=nx)$height
  if(df$dType[idx]=="W_shredder"){
    df$steepness[idx] = (1/2)*sqrt(1/df$D[idx])
  } else if (df$dType[idx]=="X_shredder"){
    df$steepness[idx] =       sqrt(1/df$D[idx])
  }
  # system(paste0("mv *-", fname_out, " PDE_output/"))
  setTxtProgressBar(pb, idx)
}
close(pb)

### define the 4 outcomes
df$CC = (df$chase==TRUE) & (df$crash==TRUE)
df$CN = (df$chase==TRUE) & (df$crash==FALSE)
df$NC = (df$chase==FALSE) & (df$crash==TRUE)
df$NN = (df$chase==FALSE) & (df$crash==FALSE)

### save the output array
saveRDS(df, file="numerical_PDE_solns.rds")

### PDE bar plotting function
func_plot_PDE_bars = function(df, id="", plot.out=TRUE){
  vec_names_outcomes = c("CC", "CN", "NC", "NN")
  vec_labels_outcomes = c("P(Chase - Crash)", "P(Chase - No_crash)", "P(No_chase - Crash)", "P(No_chase - No_crash)")
  vec_colours = c("#a6cee3", "#b2df8a", "#ca0020", "#f4a582", "#0571b0", "#92c5de")

  vec_names_variables = c("D", "r", "c", "dType")
  vec_labels_variables = c("Dispersion parameter\n(D)",
                          "Intrinsic growth rate\n(r)",
                          "Gene drive conversion rate\n(c)",
                          "Suppression gene drive system")
  for (i in 1:length(vec_names_variables)){
    ### test: i=1
    variable = vec_names_variables[i]
    variable_label = vec_labels_variables[i]
    variable_levels = sort(unique(eval(parse(text=paste0("df$", variable)))))
    if (plot.out) {if (id != "") {id_dashed = paste0(id, "-")} else {id_dashed = id}; svg(paste0("PDE_barplots-", id_dashed, variable, ".svg"), width=15, height=7)}
    mat_MU = matrix(NA, nrow=length(vec_names_outcomes), ncol=length(variable_levels))
    mat_SD = mat_MU
    for (j in 1:length(vec_names_outcomes)){
      ### test: j=3
      outcome = vec_names_outcomes[j]
      outcome_label = vec_labels_outcomes[j]
      MU = eval(parse(text=paste0("aggregate(", outcome, " ~ ", variable, ", data=df, FUN=mean)")))
      MU = MU[order(MU[,1]),]
      SD = eval(parse(text=paste0("aggregate(", outcome, " ~ ", variable, ", data=df, FUN=function(x){sd(x)/sqrt(length(x)-1)})")))
      SD = SD[order(SD[,1]),]
      mat_MU[j,] = MU[,2]
      mat_SD[j,] = SD[,2]
    }
    vec_ylim = c(0, max(c(1.0, max(mat_MU)+0.1)))
    bp = barplot(mat_MU, beside=TRUE,
                names.arg=MU[,1],
                xlab=variable_label,
                ylab="Probability",
                ylim=vec_ylim,
                col=vec_colours[3:6],
                border=NA)
    grid()
    for (k in 1:nrow(mat_MU)){
      for (l in 1:ncol(mat_MU)){
        text(x=bp[k,l], y=mat_MU[k,l], lab=round(mat_MU[k,l],3), pos=3, col=vec_colours[k+2])
        # arrows(x0=bp[i,j], y0=mat_OUTCOMES_mu[i,j]-mat_OUTCOMES_se[i,j], y1=mat_OUTCOMES_mu[i,j]+mat_OUTCOMES_se[i,j], code=3, angle=90)
      }
    }
    legend("top", legend=vec_labels_outcomes, horiz=FALSE, fill=vec_colours[3:6], bty="n")
    if (plot.out) {dev.off()}
  }
  return(bp)
}

### barplots
func_plot_PDE_bars(df=df, id="", plot.out=TRUE)
func_plot_PDE_bars(df=df[df$dType=="W_shredder", ], id="W_shredder", plot.out=TRUE)
func_plot_PDE_bars(df=df[df$dType=="X_shredder", ], id="X_shredder", plot.out=TRUE)
# func_plot_PDE_bars(df=df[df$dType=="TADS", ], id="TADS", plot.out=TRUE)

#############################################
### Analysing numerical simulation output ###
#############################################
df = readRDS("numerical_PDE_solns.rds")
vec_names_outcomes = c("CC", "CN", "NC", "NN")
vec_labels_outcomes = c("P(Chase - Crash)", "P(Chase - No_crash)", "P(No_chase - Crash)", "P(No_chase - No_crash)")

svg("PDE_numerical_solutions_D1_c1.svg", width=10, height=7)
idx_W = (df$D==10) & (df$c==1) & (df$dType=="W_shredder")
sub_df_W = df[idx_W, ]
idx_X = (df$D==10) & (df$c==1) & (df$dType=="X_shredder")
sub_df_X = df[idx_X, ]
par(mfrow=c(2,2))
for (i in 1:length(vec_names_outcomes)){
  outcome = vec_names_outcomes[i]
  outcome_label = vec_labels_outcomes[i]
  y_W = eval(parse(text=paste0("sub_df_W$", outcome)))
  y_X = eval(parse(text=paste0("sub_df_X$", outcome)))
  plot(x=sub_df_W$r, y=y_W, type="n", xlab="r", ylab="Boolean", main=outcome_label, yaxt="n")
  grid()
  axis(side=2, at=c(0,1), lab=c("FALSE", "TRUE"))
  lines(x=sub_df_W$r, y=y_W, col="black")
  lines(x=sub_df_X$r, y=y_X, col="red")
  legend("right", legend=c("W-shredder", "X-shredder"), lty=1, lwd=2, col=c("black", "red"))
}
dev.off()

