### Run simulations

####################
### Script input ###
####################
args = commandArgs(trailing=TRUE)
# args = c("/home/jeff/Documents/drivechasingecol/src", "/home/jeff/Documents/drivechasingecol/src", "5661")
src_dir = args[1]
out_dir = args[2]
iter = as.numeric(args[3])
# time \
# Rscript /home/jeff/Documents/drivechasingecol/src/driveChaseEcolRun.R \
#             /home/jeff/Documents/drivechasingecol/src \
#             /home/jeff/Documents/drivechasingecol/ms/plots/STOCHASTIC \
#             5401 # 15401 # 25401

######################
### Load functions ###
######################
setwd(src_dir)
source(paste0(src_dir, "/driveChaseEcolFunctions.R"))

########################
### Output directory ###
########################
output_dir = paste0(out_dir, "/driveChaseEcol-OUTPUT/")
dir.create(output_dir, showWarnings=FALSE)
setwd(output_dir)

##################
### Parameters ###
##################
nGens = 1000
Nstar = 5
bw = 1
torus = FALSE ### TRUE  ### for unbouded landscape
absorb = TRUE ### FALSE ### for unbouded landscape
x_min = 0
y_min = 0
x_max = 2500
y_max = 0
nInit = NULL
burnIn = 0
init_WT_x_range = c((x_max*0.10), (x_max*0.50))
init_WT_y_range = c(NULL, NULL)
init_drive_x_range = c(NULL, NULL)
init_drive_y_range = c(NULL, NULL)
conversion = 1.00
fIntr = 0.01
verbose = TRUE
plot.show = FALSE
vid.out = FALSE
pop.out = FALSE
### Miscellaneous parameters:
nReps = 100
landscape="Absorbing_boundaries"

#################
### Variables ###
#################
Rmax_max = 10
vec_Rmax = seq(2, Rmax_max, length=10)
vec_sigma = seq(2, 20, length=10)
vec_drive_type = c("W_shredder", "X_shredder", "TADS")

##################################
### Generate simulation inputs ###
##################################
### i.e. a dataframe which points to the index of each variable element 
### which we will be iterating across all combinations of
df_iterator = data.frame(expand.grid(rep = 1:nReps,
                                     Rmax = 1:length(vec_Rmax),
                                     sigma = 1:length(vec_sigma),
                                     drive_type = 1:length(vec_drive_type)))

################
### Simulate ###
################
print("####################################################")
print(paste0("Iteration: ", iter))
set.seed(iter)
print("####################################################")
### identify the indices corresponding to the values of the input variable combinations at each iteration
idx_rep        = df_iterator$rep[iter]
idx_Rmax       = df_iterator$Rmax[iter]
idx_sigma      = df_iterator$sigma[iter]
idx_drive_type = df_iterator$drive_type[iter]
### identify the values of the input variables
rep             = c(1:nReps)[idx_rep]
Rmax             = vec_Rmax[idx_Rmax]
sigma            = vec_sigma[idx_sigma]
drive_type       = vec_drive_type[idx_drive_type]
### simulate
out = evolve(nGens=nGens,
             Rmax=Rmax,
             Nstar=Nstar,
             sigma=sigma,
             bw=bw,
             torus=torus,
             absorb=absorb,
             x_min=x_min,
             y_min=y_min,
             x_max=x_max,
             y_max=y_max,
             nInit=nInit,
             burnIn=burnIn,
             init_WT_x_range=init_WT_x_range,
             init_WT_y_range=init_WT_y_range,
             init_drive_x_range=init_drive_x_range,
             init_drive_y_range=init_drive_y_range,
             drive_type=drive_type,
             conversion=conversion,
             fIntr=fIntr,
             verbose=verbose,
             plot.show=plot.show,
             vid.out=vid.out,
             vid.id=iter,
             pop.out=pop.out)
### summarise and write output
extract.metrics.and.plot(out=out,
                         output_dir=output_dir,
                         id=iter,
                         rep=rep,
                         plot.time.series=plot.show,
                         plot.velocities=plot.show,
                         vid.out=vid.out,
                         Rmax=Rmax,
                         Nstar=Nstar,
                         sigma=sigma,
                         bw=bw,
                         drive_type=drive_type,
                         conversion=conversion,
                         landscape=landscape)

############
### NOTE ###
############
### Compile `PointMetricsDriveChaseEcol2D.c` first before running in parallel:
### sinteractive -p interactive --time=01:00:00
### module load r/4.0.0
### `R CMD SHLIB -L/usr/lib64/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4/ PointMetricsDriveChaseEcol2D.c`