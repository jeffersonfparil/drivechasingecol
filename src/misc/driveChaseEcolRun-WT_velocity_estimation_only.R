### Run simulations

####################################################
### Set working directory and load the functions ###
####################################################
args = commandArgs(trailing=TRUE)
# args = c("/data/cephfs/punim1173/driveChase/", "123") # args = c("/homevol/speel/", "123")
homedir = args[1]
iter = as.numeric(args[2])
setwd(paste0(homedir, "drivechasingecol/src"))
source("driveChaseEcolFunctions.R")

########################
### Output directory ###
########################
output_dir = paste0(homedir, "driveChaseEcol-OUTPUT/")
system(paste0("mkdir ", output_dir))

##################
### Parameters ###
##################
nGens = 50
nReps = 100
nInit = 1000
burnIn = 0
collect = TRUE
dType = "Gaussian"
x_min = 0
y_min = 0
x_max = 2500 #50
y_max = 0 #50
fracN = 0.1 ### fraction of the equilibrium population size at which the gene drive is re-introduced if re_release_gen=="auto"
plot.show = FALSE
vid.out = FALSE
verbose = TRUE
K.out = FALSE
pop.out = FALSE
loss.break = TRUE
one.dimensional = TRUE #FALSE

#################
### Variables ###
#################
vec_Nstar = 5
vec_Rmax = seq(2, 5, length=13)
vec_sigma = seq(2, 20, length=10) #seq(1, 10, length=10)
vec_bw = 1
vec_drive_type = c("W_shredder", "X_shredder", "TADS")
vec_conversion = 1
vec_introduced = 0.01
mat_torus_x_absorb = matrix(c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE), ncol=2, byrow=TRUE)
mat_torus_x_absorb = mat_torus_x_absorb[2, , drop=FALSE]
    vec_landscape = c("Flat_torus", "Absorbing_boundaries", "Reflecting_boundaries")
    vec_landscape = vec_landscape[2]
mat_init_release = matrix(c(0, 1, x_max/2, 6), ncol=2, byrow=TRUE)
mat_init_release = mat_init_release[1, , drop=FALSE]
    vec_release = c("One_release_point", "Six_release_points")
    vec_release = vec_release[1]
vec_re_release = c(Inf, "auto")
vec_re_release = vec_re_release[1]

##################################
### Generate simulation inputs ###
##################################
### i.e. a dataframe which points to the index of each variable element 
### which we will be iterating across all combinations of
df_iterator = data.frame(expand.grid(rep = 1:nReps,
                                     Rmax = 1:length(vec_Rmax),
                                     Nstar = 1:length(vec_Nstar),
                                     sigma = 1:length(vec_sigma),
                                     bw = 1:length(vec_bw),
                                     drive_type = 1:length(vec_drive_type),
                                     conversion = 1:length(vec_conversion),
                                     introduced = 1:length(vec_introduced),
                                     landscape = 1:length(vec_landscape),
                                     release = 1:length(vec_release),
                                     re_release_gen = 1:length(vec_re_release)))

################
### Simulate ###
################
print("####################################################")
print(paste0("Iteration: ", iter))
set.seed(iter) ### add 1 to prevent completely random population collapse leading to simulation failure
print("####################################################")
### identify the indices corresponding to the values of the input variable combinations at each iteration
idx_rep        = df_iterator$rep[iter]
idx_Rmax       = df_iterator$Rmax[iter]
idx_Nstar      = df_iterator$Nstar[iter]
idx_sigma      = df_iterator$sigma[iter]
idx_bw         = df_iterator$bw[iter]
idx_drive_type = df_iterator$drive_type[iter]
idx_conversion = df_iterator$conversion[iter]
idx_introduced = df_iterator$introduced[iter]
idx_landscape  = df_iterator$landscape[iter]
idx_release    = df_iterator$release[iter]
idx_reRelease  = df_iterator$re_release_gen[iter]
### identify the values of the input variables
introduced       = vec_introduced[idx_introduced]
Rmax             = vec_Rmax[idx_Rmax]
Nstar            = vec_Nstar[idx_Nstar]
sigma            = vec_sigma[idx_sigma]
bw               = vec_bw[idx_bw]
drive_type       = vec_drive_type[idx_drive_type]
conversion       = vec_conversion[idx_conversion]
torus            = mat_torus_x_absorb[idx_landscape,1]
absorb           = mat_torus_x_absorb[idx_landscape,2]
release_radius   = mat_init_release[idx_release,1]
n_release_points = mat_init_release[idx_release,2]
re_release_gen   = vec_re_release[idx_reRelease]
### labels for landscape and initial release
landscape        = vec_landscape[idx_landscape]
release          = vec_release[idx_release]
### simulate
out = evolve(nGens=nGens,
            fIntr=introduced,
            Rmax=Rmax,
            Nstar=Nstar,
            nInit=nInit,
            burnIn=burnIn,
            collect=collect,
            dType=dType,
            sigma=sigma,
            bw=bw,
            drive_type=drive_type,
            conversion=conversion,
            torus=torus,
            absorb=absorb,
            x_min=x_min,
            y_min=y_min,
            x_max=x_max,
            y_max=y_max,
            one.dimensional=one.dimensional,
            release_radius=release_radius,
            n_release_points=n_release_points,
            re_release_gen=re_release_gen,
            fracN=fracN,
            plot.show=plot.show,
            vid.out=vid.out,
            vid.id=iter,
            verbose=verbose,
            K.out=K.out,
            pop.out=pop.out,
            loss.break=loss.break,
            WT_wave_stats=TRUE,
            init_WT_x_min=x_max/2,
            init_WT_y_min=y_max/2,
            init_WT_x_max=x_max/2,
            init_WT_y_max=y_max/2)
### summarise and write output
extract.metrics.and.plot(out=out, 
            output_dir=output_dir,
            df_iterator=df_iterator,
            iter=iter,
            idx_rep=idx_rep,
            plot.time.series=plot.show,
            nGens=nGens,
            fIntr=introduced,
            Rmax=Rmax,
            Nstar=Nstar,
            nInit=nInit,
            burnIn=burnIn,
            collect=collect,
            dType=dType,
            sigma=sigma,
            bw=bw,
            drive_type=drive_type,
            conversion=conversion,
            torus=torus,
            absorb=absorb,
            x_max=x_max,
            y_max=y_max,
            release_radius=release_radius,
            n_release_points=n_release_points,
            re_release_gen=re_release_gen,
            fracN=fracN,
            plot.show=plot.show,
            vid.out=vid.out,
            vid.id=iter,
            verbose=verbose,
            K.out=K.out,
            pop.out=pop.out,
            loss.break=loss.break,
            WT_wave_stats=TRUE)
