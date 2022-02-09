### Analysis of stochastioc simulation output

one.dimensional = FALSE

##########################
### Concatenate output ###
##########################
### prepare header
system("echo 'rep,Rmax,Nstar,sigma,bw,drive_type,conversion,introduced,landscape,release,re_release_gen,n_mu,n_sd,n_min,n_max,n_init,n_fin,q_mu,q_sd,q_min,q_max,q_init,q_fin,wWT_mu,wWT_sd,wWT_min,wWT_max,wWT_init,wWT_fin,wDrive_mu,wDrive_sd,wDrive_min,wDrive_max,wDrive_init,wDrive_fin,v_x_Drive_mu,v_x_Drive_sd,v_x_Drive_min,v_x_Drive_max,v_x_Drive_init,v_x_Drive_fin,v_y_Drive_mu,v_y_Drive_sd,v_y_Drive_min,v_y_Drive_max,v_y_Drive_init,v_y_Drive_fin,l_x_Drive_mu,l_x_Drive_sd,l_x_Drive_min,l_x_Drive_max,l_x_Drive_init,l_x_Drive_fin,l_y_Drive_mu,l_y_Drive_sd,l_y_Drive_min,l_y_Drive_max,l_y_Drive_init,l_y_Drive_fin,h_Drive_mu,h_Drive_sd,h_Drive_min,h_Drive_max,h_Drive_init,h_Drive_fin,v_x_WT_mu,v_x_WT_sd,v_x_WT_min,v_x_WT_max,v_x_WT_init,v_x_WT_fin,v_y_WT_mu,v_y_WT_sd,v_y_WT_min,v_y_WT_max,v_y_WT_init,v_y_WT_fin,l_x_WT_mu,l_x_WT_sd,l_x_WT_min,l_x_WT_max,l_x_WT_init,l_x_WT_fin,l_y_WT_mu,l_y_WT_sd,l_y_WT_min,l_y_WT_max,l_y_WT_init,l_y_WT_fin,h_WT_mu,h_WT_sd,h_WT_min,h_WT_max,h_WT_init,h_WT_fin,dn_mu,dn_sd,dn_min,dn_max,dn_init,dn_fin,dq_mu,dq_sd,dq_min,dq_max,dq_init,dq_fin,final_gen,crashed,chasing' > merged_sensitivity_analysis_input.csv")
### concatenate together with the header
nReps = 100
for (i in 1:nReps){
  system(paste0("cat drivechasingecol-iter.*-rep.", i, "-*.csv >> merged_sensitivity_analysis_input.csv"))
}
dat = read.csv("merged_sensitivity_analysis_input.csv")


### extract wild-type velocities, and width per Rmax-sigma combination
if (one.dimensional) {
  v_mu = dat$v_x_WT_mu
  l_mu = dat$l_x_WT_mu
} else {
  v_mu = apply(cbind(dat$v_x_WT_mu, dat$v_y_WT_mu), MAR=1, FUN=mean)
  l_mu = apply(cbind(dat$l_x_WT_mu, dat$l_y_WT_mu), MAR=1, FUN=mean)
}
df_output = data.frame(Rmax=dat$Rmax, sigma=dat$sigma, v_WT_mu=v_mu, l_WT_mu=l_mu)
if (one.dimensional) {
  saveRDS(df_output, file="1D-WT-velocity-and-width-Rmax_vs_sigma.rds")
} else {
  saveRDS(df_output, file="2D-WT-velocity-and-width-Rmax_vs_sigma.rds")
}
