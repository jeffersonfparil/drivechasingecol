######################################################
### STOCHASTIC INDIVIDUAL-BASED SPATIALLY-EXPLICIT ###
### SUPPRESSION GENE DRIVE SIMULATIONS (1D and 2D) ###
######################################################
### Authors: Ben Phillips, and Jeff Paril

### Compile and load C functions
where<-Sys.info()["sysname"]
if (where=="Darwin"){
  system(paste("R CMD SHLIB src/PointMetricsDriveChaseEcol2D.c"))
}
if (where=="Linux"){
  system(paste("R CMD SHLIB -L/usr/lib64/ -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4/ PointMetricsDriveChaseEcol2D.c"))
}
dyn.load("PointMetricsDriveChaseEcol2D.so")

############################
### INTIALISE POPULATION ###
############################
### Initialise a matrix with n individuals + an array for dispersal genotypes
### NOTE: for 1D simulations set y_min == y_max
##################################################################
### TEST:
# nInit=12000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test = initInds(nInit, x_min, y_min, x_max, y_max)
# summary(test)
##################################################################
initInds <- function(nInit, x_min=0, y_min=0, x_max=50, y_max=50){
  ### Define x-axis coordinates
  X = runif(nInit, x_min, x_max)
  ### Define y-axis coordinates
  Y = runif(nInit, y_min, y_max)
  ### Define wild-type genotype sexes (50:50 m:f), drive state (all WT), and TADS target genotype (all WT)
  Z0 = rbinom(nInit, 1, 0.5) # male=1; female=0
  zState = rep(0, nInit) # drive state, so either wildtype = 0; heterozygous = 1; homozygous drive = 2
  TADS_target = rep(0, nInit) ### TADS spermatogenesis target gene: 0 = homozygous wild type functional; 2 = drive-cleaved non-functional without drive allele; 1 = heterozygote
  ### output the initial population as a matrix
  return(cbind(X, Y, Z0, zState, TADS_target))
}

######################
### IDENTIFY MALES ###
######################
### Identify the males in ZW + drive allele system. 
###  - Z0 == 0 are females
###  - Z0 == 1 are males
#############################################################################
### TEST:
# nInit=12000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# for (i in 1:4){
#   print(mean(sexDet(pop=test.pop, drive_type=c("W_shredder", "X_shredder", "TADS")[i])))
# }
#############################################################################
sexDet <- function(pop, drive_type=c("W_shredder", "X_shredder", "TADS")[1]){
  switch(drive_type,
    W_shredder = {
      (pop[, "Z0"]==1)
    },
    X_shredder = {
      (pop[, "Z0"]==1)
    },
    TADS = {
      (pop[, "Z0"]==1)
    }
  )
}

############################
### INDIVIDUAL DENSITIES ###
############################
### Calculate the density of individuals around each individual using their pairwise distances
###   - Calculate the probability density on a bivariate normal distribution given the
###     pairwise toroidal or Euclidean distances between individuals
###     (i.e. distances on a toroidal surface or on a flat space, respectively).
###   - Specifically, calculate the density around eah individual ii
###     by calculating the pairwise toroidal distance
###     between individual ii and all the other individuals,
###     and using it to calculate the density on a bivariate normal distribution.
###   - Hence, as the distance increases (z = sqrt(|x|^2+|y|^2)),
###   	the distance from the mean of the distribution also increases, and
###   	the density decreases.
###   - The density around individual ii is the sum of densities 
###     across all pairwise distances with the other individuals.
### NOTE: for 1D simulations set y_min == y_max
###############################################################################
### TEST:
# nInit=12000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# bw=1; torus=TRUE ### torus=FALSE
# test = metrics(pop=test.pop, bw, torus, x_min, y_min, x_max, y_max)
# par(mfrow=c(1,2))
# plot(test.pop)
# plot(test)
###############################################################################
metrics <- function(pop, bw, torus=TRUE, x_min=0, y_min=0, x_max=50, y_max=50){
  if (is.null(pop)) return(NULL)
  out = .Call("metrics", R_X=pop[,"X"], 
             R_Y=pop[,"Y"],
             R_n=nrow(pop),
             R_bw=bw,
             R_torus=torus,
             R_x_min=x_min,
             R_y_min=y_min,
             R_x_max=x_max,
             R_y_max=y_max)
  return(out)
}

#####################################
### BEVERTON-HOLT GROWTH FUNCTION ###
#####################################
### Beverton-Holt growth function (calculates how many offisprings each female can generate)
###   - N = density around an individual (i.e. output of metrics())
###   - Rmax = max number of offspring an individual can have
###   - Nstar = carrying capasity given Rmax
####################################
### TEST:
# nInit=12000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# bw=1; torus=TRUE ### torus=FALSE
# test.density = metrics(pop=test.pop, bw, torus, x_min, y_min, x_max, y_max)
# Rmax=3; Nstar=5
# test = bevHolt(N=test.density, Rmax, Nstar)
# hist(test)
####################################
bevHolt <- function(N, Rmax, Nstar){
  a = (Rmax-2) / (2*Nstar) #calculates alpha for respective Rmax and Nstar values
  return(Rmax / (1 + (a*N)))
}

##################################
### FIND MALE MATES PER FEMALE ###
##################################
### For each female find a male mate (each male can mate with at least 1 female)
###   - Given a popmatrix, finds a mate for each individual, 2D case
###   - Identify the male mate of each female, 
###       - if nM > nF then output is of length nF with NAs <= nF
###       - if nM < nF then output is of length nF with NAs <= nM
###       - also if nM < nF, allows for the same male to mate with multiple females
### NOTE: for 1D simulations set y_min == y_max
######################################################################################
### TEST:
# nInit=12000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# bw=1; torus=TRUE ### torus=FALSE
# test.density = metrics(pop=test.pop, bw, torus, x_min, y_min, x_max, y_max)
# fPop = test.pop[test.pop[,"Z0"]==0, ]; nF = nrow(fPop)
# mPop = test.pop[test.pop[,"Z0"]==1, ]; nM = nrow(mPop)
# test = mate(fPop, mPop, nF, nM, bw, torus, x_min, y_min, x_max, y_max)
# print(length(test)==nF)
# print(table(test)[max(table(test))==table(test)])
######################################################################################
mate <- function(fPop, mPop, nF, nM, bw, torus, x_min=0, y_min=0, x_max=50, y_max=50){
  if (is.null(fPop) | is.null(mPop)) return(NULL)
  out<-.Call("mate", R_fX=fPop[,"X"], R_fY=fPop[,"Y"], 
             R_mX=mPop[,"X"], R_mY=mPop[,"Y"], 
             R_nF=nF, R_nM=nM,
             R_bw=bw,
             R_torus=torus,
             R_x_min=x_min,
             R_y_min=y_min,
             R_x_max=x_max,
             R_y_max=y_max)
  out[out==-99]<-NA # catch lonely ones
  return(out)
}

##################################
### TOROIDAL SPACE COORDINATES ###
##################################
### Flat toroidal space coordinates
###   - this space wraps around on itself, i.e. edgeless
### NOTE: for 1D simulations set y_min == y_max
##################################################################
### TEST:
# nInit=50; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# test = torus_coor(pop=test.pop, x_min, y_min, x_max/2, y_max/2)
# summary(test)
# par(mfrow=c(1,2))
# plot(test.pop[,1:2], type="n"); text(test.pop[,1:2], lab=1:nInit); grid()
# plot(test[,1:2], type="n"); text(test[,1:2], lab=1:nInit); grid()
##################################################################
torus_coor <- function(pop, x_min=0, y_min=0, x_max=50, y_max=50){
  pop[,"X"] = (pop[,"X"] %% (x_max-x_min)) + x_min
  pop[,"Y"] = (pop[,"Y"] %% max(c(1, (y_max-y_min)))) + y_min ### set one-dimelsional landscapes to be mod 1
  return(pop)
}

#########################################################
### NON-TOROIDAL COORDINATES WITH BOUNDARY CONDITIONS ###
#########################################################
### Non-toroidal coordinates where boundaries either absorb or reflect
###   - absorb: individuals going beyond the boundaries are lost
###   - reflect: individuals cannot move beyond the boundaries
### NOTE: for 1D simulations set y_min == y_max
##################################################################################
### TEST:
# nInit=50; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# absorb=TRUE # absorb=FALSE
# test = nonTorus_coor(pop=test.pop, x_min, y_min, x_max/2, y_max/2, absorb)
# summary(test)
# par(mfrow=c(1,2))
# plot(test.pop[,1:2], type="n"); text(test.pop[,1:2], lab=1:nInit); grid()
# plot(test[,1:2], type="n"); text(test[,1:2], lab=1:nInit); grid()
##################################################################################
nonTorus_coor <- function(pop, x_min=0, y_min=0, x_max=50, y_max=50, absorb=TRUE){
  idx_X_beyond_min = pop[,"X"] < x_min
  idx_Y_beyond_min = pop[,"Y"] < y_min
  idx_X_beyond_max = pop[,"X"] > x_max
  idx_Y_beyond_max = pop[,"Y"] > y_max
  if (absorb){
    pop[idx_X_beyond_min, "X"] = NA
    pop[idx_Y_beyond_min, "Y"] = NA
    pop[idx_X_beyond_max, "X"] = NA
    pop[idx_Y_beyond_max, "Y"] = NA
  } else {
    pop[idx_X_beyond_min, "X"] = x_min
    pop[idx_Y_beyond_min, "Y"] = y_min
    pop[idx_X_beyond_max, "X"] = x_max
    pop[idx_Y_beyond_max, "Y"] = y_max
  }
  return(pop)
}

###########################
### DISPERSAL FUNCTIONS ###
###########################
### Dispersal distance functions
### Randomly sample the distance travelled by an individual from its original location as a:
### - normal distribution
### - uniform distribution
##############################
### TEST:
# nInit=1000; x_min=0; y_min=0; x_max=10; y_max=10
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# sigma = 1
# test = gDisp(pop=test.pop, sigma)
# # test = uDisp(pop=test.pop, sigma)
# par(mfrow=c(1,2))
# plot(test.pop[, 1:2])
# hist(test)
##############################
gDisp <- function(pop, sigma){
  rnorm(nrow(pop), mean=0, sd=sigma) #dispersal distances
}
uDisp <- function(pop, min=-25, max=+25){
  runif(nrow(pop), min=min, max=max) #dispersal distances
}

############################
### INDIVIDUAL DISPESION ###
############################
### Random individual dispersion
###   - Calculate new x-y coordinates
###   - Individuals move within a circular area with radius defined by the
###     displacement distance 'disp~N(0,sigma)',
###   - and direction of movement defined by
###     direction 'dirn = (0, 2*pi)'
###     such that cos(dirn) and sin(dirn) ranges from -1 to 1;
###   - And given the trigonometric functions
###     'cos(dirn) = dx / disp' and 'cos(dirn) = dy / disp',
###   - the steps in the left-right and up-down directions are caluculated.
### NOTE: for 1D simulations set y_min == y_max
##############################################################################################################
### TEST:
# nInit=50; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# torus=TRUE # torus=FALSE
# absorb=TRUE # absorb=FALSE
# # dType="Gaussian"; sigma=5
# # test = disperse(pop=test.pop, dType, torus, absorb, x_min, y_min, x_max, y_max, sigma=sigma)
# dType="Uniform"; min=-10; max=+10
# test = disperse(pop=test.pop, dType, torus, absorb, x_min, y_min, x_max, y_max, min=min, max=max)
# summary(test)
# par(mfrow=c(1,2))
# plot(test.pop[,1:2], type="n"); text(test.pop[,1:2], lab=1:nInit); grid()
# plot(test[,1:2], type="n"); text(test[,1:2], lab=1:nInit); grid()
##############################################################################################################
disperse<-function(pop, dType="Gaussian", torus=TRUE, absorb=TRUE, x_min=0, y_min=0, x_max=50, y_max=50, ...){
  ### detect if we are simulating 1D
  if (y_min == y_max){
    one.dimensional=TRUE
  } else {
    one.dimensional=FALSE
  }
  ### choose dispersal type
  switch (dType,
          Gaussian = {disp<-gDisp(pop, ...)},
          Uniform = {disp<-uDisp(pop, ...)}
  )
  ### if we are simulating a one-dimensional space then do not disperse along the y-axis
  if (one.dimensional==TRUE){
    pop[,"X"] = pop[,"X"] + disp
  } else {
    dirn<-runif(nrow(pop), 0, 2*pi)
    pop[,"X"]<-pop[,"X"]+disp*cos(dirn)
    pop[,"Y"]<-pop[,"Y"]+disp*sin(dirn)
  }
  if (torus){
    out = torus_coor(pop, x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
  } else {
    out = nonTorus_coor(pop, absorb=absorb, x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
  }
  return(out)
}

###########################################
### SUPPRESSION GENE DRIVE INTRODUCTION ###
###########################################
### Introduce the drive heterozygote individuals with Gaussian dispersion
### NOTE: for 1D simulations set y_min == y_max
##############################################################################################################################################################
### TEST:
# nInit=100; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# test.pop = initInds(nInit, x_min, y_min, x_max, y_max)
# nIntr=10
# drive_type=c("W_shredder", "X_shredder", "TADS")[3]
# sigma = 10
# torus=TRUE # torus=FALSE
# absorb=TRUE # absorb=FALSE
# test = introDrive(nIntr, drive_type, sigma, torus, absorb, x_min, y_min, x_max/2, y_max/2)
# dim(test) ### can be less than nIntr if torus==FALSE & absorb==TRUE
# plot(test.pop[,1:2], pch=20, col=rgb(0.5,0.5,0.5,alpha=0.5)); grid()
# points(test[,1:2], pch=20, col="red")
##############################################################################################################################################################
introDrive <- function(nIntr, drive_type=c("W_shredder", "X_shredder", "TADS")[1], sigma=1, torus=TRUE, absorb=FALSE, x_min=24, y_min=24, x_max=26, y_max=26){
  ### initialise the drive heterozygotes in the middles of the spatial range
  dInds = initInds(nIntr, x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max)
  dInds[,"zState"] = 1
  ### If the drive system is X-shredder then only males can be carriers
  if (drive_type=="X_shredder"){
    dInds[,"Z0"] = 1
  }
  ### Gaussian dispersion of drive-carriers 
  out = disperse(dInds, dType="Gaussian", torus=torus, absorb=absorb, x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max, sigma=sigma)
  ### remove individuals moving beyond the landscape boundaries, i.e. NAs when torus==FALSE & absorb==TRUE
  out = out[complete.cases(out), , drop=FALSE]
  return(out)
}

#####################################
### SUPPRESSION GENE DRIVE ACTION ###
#####################################
### Gene drive conversion
###   - convert gene drive heterozygotes into gene drive homozygotes at conversion rate in ZW-to-ZZ, W-shredder and X-shredder systems
###   - convert the target spermatogenesis wild type alleles independently in each gamete given drive allele (in heterozygotes) in TADS
##########################################################################################
### TEST:
# nInit=100; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# nIntr=10
# drive_type=c("W_shredder", "X_shredder", "TADS")[3]
# sigma = 10
# torus=TRUE # torus=FALSE
# absorb=TRUE # absorb=FALSE
# conversion=0.6
# test.WT = initInds(nInit, x_min, y_min, x_max, y_max)
# test.pop = rbind(test.WT, introDrive(nIntr, drive_type, sigma, torus, absorb, x_min, y_min, x_max/2, y_max/2))
# test = drive(pop=test.pop, conversion, drive_type)
# par(mfrow=c(1,2))
# plot(test.pop[,1:2], pch=20, col=rgb(0.5,0.5,0.5,alpha=0.5)); grid()
# points(test.pop[test.pop[,"zState"]==1, 1:2], pch=20, col="red")
# plot(test[,1:2], pch=20, col=rgb(0.5,0.5,0.5,alpha=0.5)); grid()
# points(test[test[,"zState"]==2, 1:2], pch=20, col="red")
##########################################################################################
drive <- function(pop, conversion=1, drive_type=c("W_shredder", "X_shredder", "TADS")[1]){
  switch(drive_type,
    W_shredder = {
      ### ZW become Z alone
      hets = pop[, "zState"]==1 # drive heterozygotes
      pop[hets, "zState"] = 1 + rbinom(n = sum(hets), size = 1, prob = conversion)
    },
    X_shredder = {
      ### XY become Y alone
      ### Males (XY) carry the drive in the Y chromosome if:
      ###   - pop[, "zState"]==1
      ### zState becomes 2 at a rate of conversion or shredding efficiency
      drive_carrier_males = pop[, "zState"]==1
      pop[drive_carrier_males, "zState"] = 1 + rbinom(n = sum(drive_carrier_males), size = 1, prob = conversion)
    },
    TADS = {
      ### In drive carrying individuals, the target allele/s is/are redendered non-functional at the conversion rate
      drive_carriers = pop[,"zState"] > 0
      pop[drive_carriers, "TADS_target"] = pop[drive_carriers, "TADS_target"] +
                                           rbinom(n=sum(drive_carriers), size=(2-pop[drive_carriers, "TADS_target"]), prob=conversion)
      # ### equivalent to below:
      # ### homozygous for the wild type target spermatogenesis gene
      # TADS_homo = pop[, "TADS_target"]==0
      # pop[(drive_carriers)&(TADS_homo), "TADS_target"] = 0 + rbinom(n=sum((drive_carriers)&(TADS_homo)), size=2, prob=conversion)
      # ### heterozygous for the target spermatogenesis gene
      # TADS_hets = pop[, "TADS_target"]==1
      # pop[(drive_carriers)&(TADS_hets), "TADS_target"] = 1 + rbinom(n=sum((drive_carriers)&(TADS_hets)), size=1, prob=conversion)
    })
  return(pop)
}

###################################
### REPRODUCTION AND DISPERSION ###
###################################
### Reproduction
### Density-dependent reproduction, gene drive expression, and  offspring dispersal.
### (1) calculate density around each indivudal with metrics()
### (2) identify the males and females
### (3) find the mates of each female with mate()
### (4) extract the densitiy around each female with a mate
### (5) calculate the expected number of offspring per female with bevHolt()
### (6) randomly sample the number of offspring per female given the expected fecundity calculated above
### (7) gene drive expression, either pre-zygotic or post-zygotic
### (8) random offspring dispersal
### NOTE: for 1D simulations set y_min == y_max
#############################################################################################################################################################################################################
### TEST:
# nInit=1000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# torus=TRUE # torus=FALSE
# absorb=TRUE # absorb=FALSE
# drive_type=c("W_shredder", "X_shredder", "TADS")[1]
# nIntr=10; conversion=0.95; Rmax=3; Nstar=5; bw=1; sigma=2
# test.pop = rbind(initInds(nInit, x_min, y_min, x_max, y_max),
#                  introDrive(nIntr, drive_type, sigma, torus, absorb, x_min, y_min, x_max/2, y_max/2))
# test = repro(pop=test.pop, NULL, Rmax, Nstar, sigma, bw, drive_type, torus, absorb, conversion, x_min, y_min, x_max, y_max)
# dim(test) ### can be less than nIntr if torus==FALSE & absorb==TRUE
# par(mfrow=c(1,2))
# plot(test.pop[,1:2], pch=20, col=rgb(0.5,0.5,0.5,alpha=0.5)); grid()
# points(test.pop[test.pop[,"zState"]>0, 1:2], pch=20, col="red")
# plot(test[,1:2], pch=20, col=rgb(0.5,0.5,0.5,alpha=0.5)); grid()
# points(test[test[,"zState"]>0, 1:2], pch=20, col="red")
#############################################################################################################################################################################################################
repro <- function(pop, spatMetrics=NULL, Rmax=3, Nstar=5, sigma=1, bw=1, drive_type=c("W_shredder", "X_shredder", "TADS")[1], torus=TRUE, absorb=TRUE, conversion=1.0, x_min=0, y_min=0, x_max=50, y_max=50){
  if (is.null(pop)) return(NULL)
  ### calculate density around each individual
  if (is.null(spatMetrics)){
    spatMetrics = metrics(pop, bw, torus, x_min, y_min, x_max, y_max)
  }
  ### identify the males (Z0==1)
  mss = sexDet(pop, drive_type)
  ### identify the females (Z0==0)
  fss = !mss
  ### number of males
  nM = sum(mss)
  ### the number of females
  nF = sum(fss)
  ### Has the population gone extinct?
  if (nF==0 | nM==0) return(NULL)
  ### extract the males
  mPop = pop[mss, , drop=FALSE]
  ### extract the females
  fPop = pop[fss, , drop=FALSE]
  ### sample neighbouring sires to mate with each female
  ###   - where each male can mate with at least 1 female
  sireInd = mate(fPop, mPop, nF, nM, bw, torus, x_min, y_min, x_max, y_max)
  ### identify only the mating pairs, i.e. singles or NAs not allowed
  matingPairsInd = !is.na(sireInd)
  ### filter the index of males keeping only those with mates
  sireInd = sireInd[matingPairsInd]
  ### filter the index of females keeping only those with mates
  damInd = c(1:nF)[matingPairsInd]
  ### extract the densitiy around each mated female
  Nx = spatMetrics[fss][matingPairsInd]
  ### number of mated females
  nF = length(Nx)
  ### reproductive output of females:
  ####### first, what is the expected number of offspring per female?
  #######         Given the density of individuals around each female (calculated by metrics),
  #######         calculate the mean number of offspring of each female.
  EW = bevHolt(Nx, Rmax, Nstar)
  ####### second, how many offsprings in total?
  #######         Given the mean number of offspring of each female,
  #######         randomly sample the number of offspring of each female from a Poisson distribution.
  W = rpois(nF, lambda=EW)
  offN = sum(W)
  ###### Finally, Has the population collapsed because the number of males or females or the number of offsprings dropped to zero?
  if (nF==0 | offN==0) return(NULL)
  ### generate the population matrices for the dams, sires and offsprings
  ### where each individual dam and sire is replicated to align with the number of offsprings per mating pair, i.e. W
  ###### dams
  damInd = rep(damInd, times=W)
  damPop = fPop[damInd, , drop=FALSE]
  ###### sires
  sireInd = rep(sireInd, times=W)
  sirePop = mPop[sireInd, , drop=FALSE]
  ###### offsprings
  offPop = damPop
  ### offspting sex determination
  ### initially set 50:50 male:female prior to gene drive expression (pre- or pos-zygotic dependeing on the gene drive system)
  offPop[,"Z0"] = rbinom(n = offN, prob = 0.5, size = 1)
  ########################
  ### W-shredder drive ###
  ########################
  if (drive_type=="W_shredder"){
    ### gene drive expression at meiosis : ZW-->Z females
    damPop = drive(damPop, conversion, drive_type) ### zState=1 becomes zState=2 which is essentially the same as shredding the W chromosome
    ### identify W-shredded female (ZW-->Z0)
    w_shred = (damPop[,"zState"]==2)
    ### generate only male offspring from W-shredded dams
    offPop[w_shred, "Z0"] = 1
    ### inherit only Z gametes from W-shredded dams
    z_allele_from_dam = as.numeric(w_shred)
    ### randomly select the wild-type or drive Z chromosome from the sire
    z_allele_from_sire = rbinom(n=nrow(sirePop), size=1, prob=sirePop[, "zState"]/2)
    ### determine the genotype of the offspring at the Z chromosome's drive locus
    offPop[,"zState"] = z_allele_from_dam + z_allele_from_sire
  }
  ########################
  ### X-shredder drive ###
  ########################
  if (drive_type=="X_shredder"){
    ### gene drive expression at meiosis : XY--> Y
    sirePop = drive(sirePop, conversion, drive_type) ### converts X-shredded individuals into zState=2 from zState=1
    ### identify X-shredded sires (XY-->Y)
    x_shred = (sirePop[,"zState"]==2)
    ### generate only male offspring which are also drive-carriers, from X-shredded sires
    offPop[x_shred, "Z0"] = 1
    offPop[x_shred, "zState"] = 1
    ### male offspring from non-X-shredded carrier sires still cary the Y-linked drive allele
    idx_males = offPop[, "Z0"] == 1
    idx_male_carriers_no_shred = (sirePop[,"zState"]==1)
    offPop[idx_males & idx_male_carriers_no_shred, "zState"] = 1
  }
  ###########################################
  ### TADS: Toxin-antidote dominant sperm ###
  ###########################################
  if (drive_type=="TADS"){
    ### no homing, instead we have conversion of the wild type target spermatogenesis gene into non-functional form at meiosis
    sirePop = drive(sirePop, conversion=conversion, drive_type="TADS") ### at zState=1 TADS_target heterozygotes become 2 or TADS_target homozygote wild types become 1 or 2 based on conversion efficiency
    ### additionally, the wild-type target alleles in the dams are converted into non-functional alleles by the drive allele (homozygous and heterozygous drive females) with a constant embry cut rate or coversion rate of 95% efficiency
    damPop = drive(damPop, conversion=0.95*conversion, drive_type="TADS")
    ### identify the drive homozygote sires which are sterile and cull-off the offspring from these sterile sires (as well as the corrensponding sire-dam combinations)
    sterile_males = (sirePop[,"zState"]==2)
    offPop = offPop[!sterile_males, , drop=FALSE]
    sirePop = sirePop[!sterile_males, , drop=FALSE]
    damPop = damPop[!sterile_males, , drop=FALSE]
    ### proceed with meiosis if we still have fertile males
    if (sum(!sterile_males)>0){
      ### identify sires homozygous for the non-functional target spermatogenesis gene and missing the drive and cull them off since they will not produce any offspring
      sterile_males = (sirePop[,"TADS_target"]==2) & (sirePop[,"zState"]==0) ### the resulting (1,0) gametes are non-viable
      offPop = offPop[!sterile_males, , drop=FALSE]
      sirePop = sirePop[!sterile_males, , drop=FALSE]
      damPop = damPop[!sterile_males, , drop=FALSE]
    }
    ### continue meiosis if we still have fertile males
    if (sum(!sterile_males)>0){
      ### randomly sample dam alleles
      target_allele_from_dam = rbinom(n=nrow(offPop), size=1, prob=(damPop[,"TADS_target"]/2))
      drive_allele_from_dam = rbinom(n=nrow(offPop), size=1, prob=(damPop[,"zState"]/2))
      ### randomly sample sire alleles, and later correct for allele combinations in the resulting viable gametes
      target_allele_from_sire = rbinom(n=nrow(offPop), size=1, prob=(sirePop[,"TADS_target"]/2))
      drive_allele_from_sire = rbinom(n=nrow(offPop), size=1, prob=(sirePop[,"zState"]/2))
      ### allele segregation in the gametes, where non-functional target allele without the drive do not mature,
      ### and only the gametes with viable allele combinations are produced by the sires
      # (sirePop[,"TADS_target"]==2) & (sirePop[,"zState"]==0) ### (1) sterile: non-functional target and no drive (x1,0)
      # (sirePop[,"TADS_target"]==2) & (sirePop[,"zState"]==1) ### (2) fertile: non-functional target and drive (1,1)
      # (sirePop[,"TADS_target"]==2) & (sirePop[,"zState"]==2) ### (3) sterile: zState==2 are sterile males
      # (sirePop[,"TADS_target"]==1) & (sirePop[,"zState"]==0) ### (4) fertile: functional target and no drive (0,0)
      # (sirePop[,"TADS_target"]==1) & (sirePop[,"zState"]==1) ### (5) fertile: functional target with/without drive, as well as non-functional target with drive (0,1), (0,0), and (1,1)
      # (sirePop[,"TADS_target"]==1) & (sirePop[,"zState"]==2) ### (6) sterile: zState==2 are sterile males
      # (sirePop[,"TADS_target"]==0) & (sirePop[,"zState"]==0) ### (7) fertile: functional target and no drive (0,0)
      # (sirePop[,"TADS_target"]==0) & (sirePop[,"zState"]==1) ### (8) fertile: functional target with/without drive (0,1), and (0,0)
      # (sirePop[,"TADS_target"]==0) & (sirePop[,"zState"]==2) ### (9) sterile: zState==2 are sterile males
      ### (2) fertile: non-functional target and drive (1,1)
      idx = (sirePop[,"TADS_target"]==2) & (sirePop[,"zState"]==1)
      target_allele_from_sire[idx] = drive_allele_from_sire[idx] = 1
      ### (4) fertile: functional target and no drive (0,0)
      idx = (sirePop[,"TADS_target"]==1) & (sirePop[,"zState"]==0)
      target_allele_from_sire[idx] = drive_allele_from_sire[idx] = 0
      ### (5) fertile: functional target with/without drive, as well as non-functional target with drive (0,1), (0,0), and (1,1)
      idx = (sirePop[,"TADS_target"]==1) & (sirePop[,"zState"]==1)
      viable_gametes = matrix(c(0,1,0,0,1,1), ncol=2, byrow=TRUE) ### arrange the viable gametes in a matrix (column 1 for the target gene and column 2 for the drive)
      viable_gametes_sample = viable_gametes[sample(x=1:3, size=sum(idx), replace=TRUE), , drop=FALSE] ### randomly sample from these viable gametes
      target_allele_from_sire[idx] = viable_gametes_sample[,1]
      drive_allele_from_sire[idx] = viable_gametes_sample[,2]
      ### inherit target spermatogenesis and drive alleles from dams and sires
      offPop[,"TADS_target"] = target_allele_from_dam + target_allele_from_sire
      offPop[,"zState"] = drive_allele_from_dam + drive_allele_from_sire
    }
  }
  ### Disperse offspring and output the offspring population
  offPop = disperse(pop=offPop, dType="Gaussian", torus, absorb, x_min, y_min, x_max, y_max, sigma=sigma)
  ### remove individuals which went beyond the landscape boundaries if we have absorbing boundaries
  idx = is.na(offPop[,"X"]) | is.na(offPop[,"Y"])
  offPop  =  offPop[!idx, , drop=FALSE]
  damPop  =  damPop[!idx, , drop=FALSE]
  sirePop = sirePop[!idx, , drop=FALSE]
  if (nrow(offPop)==0) {
    offPop = NULL
  }
  ### return offspring population
  return(offPop)
}

##################################
### DETECT FLUCTUATION SIGNALS ###
##################################
### This is Jean-Paul van Brakel's algorithm he posted on stackoverflow on 2014-03-25 [https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362]
### inputs:
###     - y [R] is the vector of time-series data
###     - lag [N] is the number of data points or the window size to compute the mean and standard deviation from
###     - threshold [+R] is the number of number of standard deviations above which a signal is detected
###     - influence [0, 1] is the influence of the data in the next window on the current mean and standard deviation
#####################################################################
### TEST:
# y = sin(seq(0, 4*pi, length=100))
# test = ThresholdingAlgo(y)
# plot(y, type="l")
# lines(test$signals, col="red")
#####################################################################
ThresholdingAlgo <- function(y, lag=10, threshold=3, influence=0.05){
    signals <- rep(0,length(y))
    filteredY <- y[0:lag]
    avgFilter <- NULL
    stdFilter <- NULL
    avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
    stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
    for (i in (lag+1):length(y)){
        if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
            if (y[i] > avgFilter[i-1]) {
                signals[i] <- 1;
            } else {
                signals[i] <- -1;
            }
            filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
        } else {
            signals[i] <- 0
            filteredY[i] <- y[i]
        }
        avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
        stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
    }
    return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

################################
### INVASION WAVE PROPERTIES ###
################################
### Invasion wave properties

### NOTE: for 1D simulations set y_min == y_max
#############################################################################################
### TEST:
# nInit=5000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# torus=FALSE
# absorb=TRUE # absorb=FALSE
# drive_type=c("W_shredder", "X_shredder", "TADS")[1]
# nIntr=nInit*0.01; conversion=0.95; Rmax=3; Nstar=5; bw=1; sigma=1
# column="zState"; nclass=20; pch=20; alpha=0.5; cex=1; gg=123; nGens=500; Rmax=3; Nstar=5
# fIntr=nIntr/nInit
# test.WT = initInds(nInit, x_max*0.1, y_max*0.1, x_max*0.5, y_max*0.5)
# test.DC = introDrive(nIntr, drive_type, sigma, torus=TRUE, absorb, x_max*0.1, y_max*0.1, x_max*0.15, y_max*0.15)
# test.pop = rbind(test.WT, test.DC)
# test.pop = repro(pop=test.pop,
#                  NULL, Rmax, Nstar, sigma, bw, drive_type, torus, absorb, conversion, x_min, y_min, x_max, y_max)
# plot(test.pop)
# edge_reached = FALSE
# Nstar = Nstar
# test = wave_properties(test.pop, sigma, bw, torus, x_min, y_min, x_max, y_max)
#############################################################################################
wave_properties <- function(pop, sigma, bw, torus, x_min, y_min, x_max, y_max, edge_reached=FALSE, Nstar=5){
  ### initialise output
  drive_wave_x_min   = drive_wave_y_min   = NA
  drive_wave_x_max   = drive_wave_y_max   = NA
  drive_wave_x_width = drive_wave_y_width = NA
  drive_wave_x_width_trailing = drive_wave_y_width_trailing = NA
  drive_wave_height  = NA
  WT_wave_x_min   = WT_wave_y_min   = NA
  WT_wave_x_max   = WT_wave_y_max   = NA
  WT_wave_x_width = WT_wave_y_width = NA
  WT_wave_height  = NA
  penetration = FALSE
  ### subset the drive-carriers
  drive_carriers = pop[pop[,"zState"]>0, , drop=FALSE]
  ### subset only the wild-types
  wild_types = pop[pop[,"zState"]==0, , drop=FALSE]
  ### calculate wave statistics if none of the edges has been reached (this is to prevent measuring the chaotic movements when the WT pentrated the drive wave and also as a result of toroidal landscapes)
  if (!edge_reached) {
    #####################################
    ### DRIVE-CARRIER WAVE PROPERTIES ###
    #####################################
    ### calculate density among drive-carriers
    y = metrics(drive_carriers, bw, torus, x_min, y_min, x_max, y_max)
    ### estimate the range, width, and height of the invading gene drive wave,
    ### if just a single wave is present,
    ### and there are not individuals at the left-most edge of the landscape
    for (spatial_axis in c("X", "Y")) {
      # spatial_axis = "X"
      ### break off the loop is we are at the y-axis and the landscape is 1D
      ### or if the drive-carriers are not yet dispersed
      ### of if there are still drive-carriers
      if (nrow(drive_carriers)<2) {break}
      if (((spatial_axis=="Y") & (y_min==y_max)) | (sd(drive_carriers[, spatial_axis], na.rm=TRUE)==0)) {break}
      ### linear interpolation of these densities across the whole spatial domain
      x = drive_carriers[, spatial_axis]
      approx_density_func = approxfun(x=x, y=y)
      x_new = c(0:x_max)
      y_pred = approx_density_func(x_new)
      y_pred[(is.na(y_pred)) | (is.infinite(y_pred))] = 0 ### set NAs (i.e. not within the range of the drive coordinates) to zero
      ### detect peaks and throughs along the spatial domain (set lag to be at least 2)
      vec_signal = ThresholdingAlgo(y=y_pred, lag=max(c(2, ceiling(x_max*0.01))), threshold=4, influence=0.05)$signal
      ### find the transition points (we expect peaks to be represented by +1 followed by -1 after none or multiple zeros)
      vec_diff = diff(vec_signal)
      ### calculate drive wave properties only if we have a single peak to avoid messiness,
      ### and there are no individuals in the lef-most edge of the landscape (since we're assuming/measuring a unidirectional wave movement)
      individual_on_right_edge = (eval(parse(text=paste0(tolower(spatial_axis), "_max"))) - max(pop[,spatial_axis])) < sigma
      if ((sum(vec_diff==+1)==1) & (sum(vec_diff==-1)==1) & !individual_on_right_edge){
        idx = y>0
        vec_peak_loc = drive_carriers[idx, spatial_axis]
        vec_peak_density = y[idx]
        wave_range_min = min(vec_peak_loc, na.rm=TRUE)
        wave_range_max = max(vec_peak_loc, na.rm=TRUE)
        wave_width = wave_range_max - wave_range_min
        wave_width_trailing = vec_peak_loc[which(vec_peak_density==max(vec_peak_density, na.rm=TRUE))[1]] - wave_range_min
        wave_height = max(vec_peak_density, na.rm=TRUE)
      } else {
        wave_range_min = NA
        wave_range_max = NA
        wave_width = NA
        wave_width_trailing = NA
        wave_height = NA
      }
      assign(paste0("drive_wave_", tolower(spatial_axis), "_min"), wave_range_min)
      assign(paste0("drive_wave_", tolower(spatial_axis), "_max"), wave_range_max)
      assign(paste0("drive_wave_", tolower(spatial_axis), "_width"), wave_width)
      assign(paste0("drive_wave_", tolower(spatial_axis), "_width_trailing"), wave_width_trailing)
      if (is.na(drive_wave_height)){
        drive_wave_height = wave_height
      }
    }
    #################################
    ### WILD-TYPE WAVE PROPERTIES ###
    #################################
    ### calculate density among drive-carriers
    y = metrics(wild_types, bw, torus, x_min, y_min, x_max, y_max)
    ### estimate the range, width, and height of the
    ### wild-type wave invading free space to the right,
    ### if and only if the WT has not yet reached the boundaries
    for (spatial_axis in c("X", "Y")) {
      # spatial_axis = "X"
      if ((spatial_axis=="Y") & (y_min==y_max)) {break}
      ### calculate WT wave proporties only if we have not reached the boundaries,
      ### and the wild-types have not penetrated the wave of drive-carriers
      WT_on_left_edge = (min(wild_types[,spatial_axis], na.rm=TRUE) - eval(parse(text=paste0(tolower(spatial_axis), "_min")))) < sigma ### It's fine if we measure this even if there is wave penetration since we're interested in WT_on_right_edge to measure WT velocity moving to the right
      WT_on_right_edge = (eval(parse(text=paste0(tolower(spatial_axis), "_max"))) - max(wild_types[,spatial_axis], na.rm=TRUE)) < sigma
      if ( !WT_on_left_edge & !WT_on_right_edge ) {
        ### find the farthest location at which density is at least at carrying capacity
        idx = (y >= Nstar)
        idx[is.na(idx)] = FALSE
        if (sum(idx)==0){
          idx = which(y==max(y, na.rm=TRUE))
          idx[is.na(idx)] = FALSE
        }
        wave_right_side_peak_loc = max(wild_types[idx, spatial_axis], na.rm=TRUE)
        idx = wild_types[, spatial_axis] > wave_right_side_peak_loc
        wave_right_side_tail_loc = max(wild_types[idx, spatial_axis], na.rm=TRUE)
        ### WT wave properties to output
        wave_range_min = min(wild_types[, spatial_axis], na.rm=TRUE) ### defined as the left-most WT
        wave_range_max = wave_right_side_tail_loc ### defined as the right-most WT
        wave_width = 2*(wave_right_side_tail_loc - wave_right_side_peak_loc) ### defined as twice the distance between the right-most peak and the right-most WT, since we have a flat top as the WT proliferates unlike the suppression gene drive-carriers
        wave_height = y[wild_types[, spatial_axis]==wave_right_side_peak_loc]
      } else {
        wave_range_min = NA
        wave_range_max = NA
        wave_width = NA
        wave_height = NA
      }
      assign(paste0("WT_wave_", tolower(spatial_axis), "_min"), wave_range_min)
      assign(paste0("WT_wave_", tolower(spatial_axis), "_max"), wave_range_max)
      assign(paste0("WT_wave_", tolower(spatial_axis), "_width"), wave_width)
      if (is.na(WT_wave_height) & (length(wave_height)==1)){
        WT_wave_height = wave_height
      }
      ############################################################
      ### Detect if at least one of the edges has been reached ###
      ############################################################
      edge_reached = edge_reached | WT_on_left_edge | WT_on_right_edge
    }
  }
  ############################################
  ### Drive wave penetration by wild-types ###
  ############################################
  if (exists("y")==FALSE){
    y = metrics(wild_types, bw, torus, x_min, y_min, x_max, y_max)
  }
  for (spatial_axis in c("X", "Y")) {
    if ((spatial_axis=="Y") & (y_min==y_max)) {break}
    WT_on_right_edge = (eval(parse(text=paste0(tolower(spatial_axis), "_max"))) - max(wild_types[,spatial_axis], na.rm=TRUE)) < sigma
    WT_to_the_left_of_drive_wave = (min(drive_carriers[,spatial_axis]) - wild_types[,spatial_axis]) > sigma
    WT_penetrated_with_desity_greater_than_2 = y[WT_to_the_left_of_drive_wave] > 2
    if (length(WT_penetrated_with_desity_greater_than_2)==0){WT_penetrated_with_desity_greater_than_2 = FALSE}
    if (torus){
      penetration = penetration | (WT_penetrated_with_desity_greater_than_2 & !WT_on_right_edge)
    } else {
      penetration = penetration | (WT_penetrated_with_desity_greater_than_2)
    } 
  }
  ##############
  ### Output ###
  ##############
  out = list(drive_wave_x_min   = drive_wave_x_min,
             drive_wave_x_max   = drive_wave_x_max,
             drive_wave_y_min   = drive_wave_y_min,
             drive_wave_y_max   = drive_wave_y_max,
             drive_wave_x_width = drive_wave_x_width,
             drive_wave_x_width_trailing = drive_wave_x_width_trailing,
             drive_wave_y_width = drive_wave_y_width,
             drive_wave_y_width_trailing = drive_wave_y_width_trailing,
             drive_wave_height  = drive_wave_height,
             WT_wave_x_min      = WT_wave_x_min,
             WT_wave_x_max      = WT_wave_x_max,
             WT_wave_y_min      = WT_wave_y_min,
             WT_wave_y_max      = WT_wave_y_max,
             WT_wave_x_width    = WT_wave_x_width,
             WT_wave_y_width    = WT_wave_y_width,
             WT_wave_height     = WT_wave_height,
             penetration        = penetration,
             edge_reached       = edge_reached)
  return(out)
}

##########################################
### PLOT POINTS AND DENSITIES IN SPACE ###
##########################################
### Plots: scatterplots, and invasion wave plots
### NOTE: for 1D simulations set y_min == y_max
############################################################################################################################################################################################################################################################################################
### TEST:
# nInit=1000; x_min=0; y_min=0; x_max=50; y_max=50
# # y_max=0 ### for 1D
# torus=TRUE # torus=FALSE
# absorb=TRUE # absorb=FALSE
# drive_type=c("W_shredder", "X_shredder", "TADS")[1]
# nIntr=10; conversion=0.95; Rmax=3; Nstar=5; bw=1; sigma=2
# column="zState"; nclass=20; pch=20; alpha=0.5; cex=1; gg=123; nGens=500; Rmax=3; Nstar=5
# fIntr=nIntr/nInit
# test.pop = repro(pop=rbind(initInds(nInit, x_min, y_min, x_max/2, y_max/2), introDrive(nIntr, drive_type, sigma, torus, absorb, x_min, y_min, x_max/4, y_max/4)),
#                  Rmax, Nstar, sigma, bw, drive_type, torus, absorb, conversion, x_min, y_min, x_max, y_max)
# pop.plot2D(pop=test.pop, NULL, column, nclass, pch, alpha, cex, gg, nGens, Rmax, Nstar, x_min, y_min, x_max, y_max, sigma, bw, torus, absorb, fIntr, nIntr, drive_type, conversion)
############################################################################################################################################################################################################################################################################################
pop.plot2D = function(pop, vec_density=NULL, column="zState", nclass=20, pch=20, alpha=0.5, cex=1, gg=1, nGens=1, Rmax=3, Nstar=5, x_min=0, y_min=0, x_max=50, y_max=50, sigma=1, bw=1, torus=TRUE, absorb=FALSE, fIntr=0.01, nIntr=NULL, drive_type=NULL, conversion=NULL, wave_prop=NULL){  
  ### detect if we are simulating 1D
  if (y_min == y_max){
    one.dimensional=TRUE
  } else {
    one.dimensional=FALSE
  }
  ### whole population densities
  if (is.null(vec_density)){
    vec_density = metrics(pop, bw, torus, x_min, y_min, x_max, y_max)
  }
  ### drive genotype colours
  col_DD = rgb(0.39,0.39,0.39, alpha=alpha) # #636363
  col_Dd = rgb(0.99,0.60,0.16, alpha=alpha) # #fe9929
  col_dd = rgb(0.87,0.18,0.15, alpha=alpha) # #de2d26
  colours = rep(col_DD, times=nrow(pop))
  colours[pop[,column]==1] = col_Dd
  colours[pop[,column]==2] = col_dd
  ########################
  ### plot area layout ###
  ########################
  if (one.dimensional==FALSE){
    layout(matrix(c(rep(1,6), rep(4,2), rep(2,3), rep(3,3), rep(4,2)), byrow=TRUE, nrow=2))
  } else {
    layout(matrix(c(rep(1,6), rep(3,2), rep(2,6), rep(3,2)),nrow=2, byrow=TRUE))
  }
  ########################################
  ### scatterplot across the landscape ###
  ########################################
  par(mar=c(5,5,1,5))
  plot(x=pop[,"X"], y=pop[,"Y"], xlim=c(0,x_max), ylim=c(0,y_max), pch=pch, col=colours, xlab="X", ylab="Y", cex=cex, cex.lab=cex/2, cex.axis=cex/2)
  grid()
  ### densities and drive allele frequencies across the spatial axis/axes
  if (one.dimensional==FALSE){
    vec_axes = c("X", "Y")
    K = Nstar*(x_max-x_min)*(y_max-y_min)/nclass
  } else {
    vec_axes = c("X")
    K = Nstar*(x_max-x_min)/nclass
  }
  ###################
  ### histogram/s ###
  ###################
  for (spatial_axis in vec_axes){
    df_1D = data.frame(x=pop[,spatial_axis], d=vec_density, q=pop[,"zState"]/2)
    df_1D = df_1D[order(df_1D$x, decreasing=FALSE), ]
    par(mar=c(5,5,1,5))
    vec_breaks = seq(x_min, x_max, length=nclass) ### define histogram breaks
    ### population density
    hist(df_1D$x, breaks=vec_breaks, xlim=c(x_min,x_max), ylim=c(0, 2*K), xlab=spatial_axis, ylab="Population density", main="", bord=NA)
    lines(x=c(x_min, x_max), y=c(K, K), col="black")
    if (length(df_1D$x)>=2){
      par(new=TRUE)
      plot(density(df_1D$x), xlim=c(x_min,x_max), col="gray", xaxt="n", yaxt="n", xlab="", ylab="", main="")
    }
    ### drive allele-carriers
    df_q = df_1D[df_1D$q > 0, ]
    vec_q = unlist(apply(df_q, MAR=1, FUN=function(x){rep(x[1], times=2*x[3])})) ### count drive-carriers proportional to the number of drive alleles they hold
    if (length(vec_q)>0){
      par(new=TRUE)
      hist(vec_q, breaks=vec_breaks, col="red", xlim=c(x_min,x_max), ylim=c(0, 2*K), xlab="", ylab="", main="", bord=NA)
    }
    if (length(vec_q)>=2){
      par(new=TRUE)
      plot(density(vec_q), xlim=c(x_min,x_max), col="red", xaxt="n", yaxt="n", xlab="", ylab="", main="")
    }
    axis(side=4, at=0, , lab=""); mtext(side=4, padj=4, text="Density")
  }
  ###############
  ### legends ###
  ###############
  par(mar=c(0,0,0,0))
  plot(0,0,type="n",xaxt="n",yaxt="n", bty="n")
  legend_1_text = c(paste0("Generation: ", gg, " / ", nGens),
                    paste0("Population size: ", nrow(pop), " (", round(sum(pop[,"Z0"]==0)*100/nrow(pop),2), "% females)"),
                    paste0("Mean density: ", round(mean(vec_density, na.rm=TRUE),4)),
                    paste0("Standard deviation of density: ", round(sd(vec_density, na.rm=TRUE),4)),
                    paste0("Drive allele frequency: ", round(sum(pop[,"zState"])/(2*nrow(pop)),4)),
                    paste0("Rmax (max. fecundity): ", Rmax),
                    paste0("Nstar (max. density): ", Nstar),
                    paste0("X range: ", x_min, "-", x_max), 
                    paste0("Y range: ", y_min, "-", y_max), 
                    paste0("Landscape: Flat Torus"), 
                    paste0("Density sigma: ", bw),
                    paste0("Dispersal sigma: ", sigma),
                    paste0("Drive type: ", drive_type),
                    paste0("Conversion rate: ", conversion),
                    paste0("Drive-carriers introduced: ", nIntr, " (", round(fIntr*100), "%)"))
                  if ((torus==FALSE) & (absorb==TRUE)){
                    legend_1_text[grepl("Landscape", legend_1_text)] = paste0("Landscape: Absorbing boundaries")
                  } else if ((torus==FALSE) & (absorb==FALSE)){
                    legend_1_text[grepl("Landscape", legend_1_text)] = paste0("Landscape: Reflecting boundaries")
                  }
  if (!is.null(wave_prop)){
    legend_1_text = c(legend_1_text,
                      paste0("Drive wave left tail position (x-axis): ", round(wave_prop$drive_wave_x_min,2)),
                      paste0("Drive wave right tail position (x-axis): ", round(wave_prop$drive_wave_x_max,2)),
                      paste0("Drive wave width (x-axis): ", round(wave_prop$drive_wave_x_width,2)),
                      paste0("Drive wave left tail position (y-axis): ", round(wave_prop$drive_wave_y_min,2)),
                      paste0("Drive wave right tail position (y-axis): ", round(wave_prop$drive_wave_y_max,2)),
                      paste0("Drive wave width (y-axis): ", round(wave_prop$drive_wave_y_width,2)),
                      paste0("Drive wave height: ", round(wave_prop$drive_wave_height,2)),
                      paste0("WT trailing wave left tail pos. (x-axis): ", round(wave_prop$WT_wave_x_min,2)),
                      paste0("WT leading wave right tail pos. (x-axis): ", round(wave_prop$WT_wave_x_max,2)),
                      paste0("WT leading wave width (x-axis): ", round(wave_prop$WT_wave_x_width,2)),
                      paste0("WT trailing wave left tail pos. (y-axis): ", round(wave_prop$WT_wave_y_min,2)),
                      paste0("WT leading wave right tail pos. (y-axis): ", round(wave_prop$WT_wave_y_max,2)),
                      paste0("WT leading wave width (y-axis): ", round(wave_prop$WT_wave_y_width,2)),
                      paste0("WT leading wave height: ", round(wave_prop$WT_wave_height,2)))
  }
  legend_2_text = c("Wild type", "Drive heterozygote", "Drive homozygote")
  legend("topleft", legend=legend_1_text, cex=cex/2, bty="n")
  legend("center", legend=legend_2_text, pch=19, col=c(col_DD,col_Dd,col_dd), cex=cex/1.5, bty="n")
}

################################
### MAIN SIMULATION FUNCTION ###
################################
### (1) Initialise a population
### (2) Perform burnin simulations without gene drive introduction until the population size equillibriate
### (3) Introduce gene drive heterozygote individuals
### (4) Re-introduce gene drive heterozygote individuals if specified
### (5) Continue reproduction until population crashes or the gene drive allele is lost
### NOTE: for 1D simulations set y_min == y_max
############################################################################################################################################################################################################################################################################
### TEST:
# ### mandatory input
# nGens=78; Rmax=3; Nstar=5; sigma=5 # dispersal_mean = sigma * sqrt(2/pi)
# ### landscape parameters
# # bw=1; torus=TRUE; absorb=TRUE; x_min=0; y_min=0; x_max=50; y_max=50
# bw=1; torus=TRUE; absorb=FALSE; x_min=0; y_min=0; x_max=1000; y_max=0
# ### initilisation parameters
# nInit=NULL; burnIn=0
# # init_WT_x_range=c(NULL, NULL)
# # init_WT_y_range=c(NULL, NULL)
# init_WT_x_range=c((x_max*0.10), (x_max*0.70))
# init_WT_y_range=c((y_max*0.10), (y_max*0.70))
# init_drive_x_range=c(NULL, NULL)
# init_drive_y_range=c(NULL, NULL)
# ### suppression gene drive parameters
# drive_type=c("W_shredder", "X_shredder", "TADS")[1]
# conversion=1.0; fIntr=0.01
# ### output paramters
# verbose=TRUE
# plot.show=TRUE
# vid.out=FALSE
# vid.id="test"
# pop.out=TRUE
# set.seed(123)
# test = evolve(nGens, Rmax, Nstar, sigma, bw, torus, absorb, x_min, y_min, x_max, y_max, nInit, burnIn, init_WT_x_range, init_WT_y_range, init_drive_x_range, init_drive_y_range, drive_type, conversion, fIntr, verbose, plot.show, vid.out, vid.id, pop.out)
############################################################################################################################################################################################################################################################################
evolve <- function(nGens, Rmax, Nstar, sigma,
                   bw=1, torus=TRUE, absorb=FALSE, x_min=0, y_min=0, x_max=50, y_max=50,
                   nInit=NULL, burnIn=10, init_WT_x_range=c(NULL, NULL), init_WT_y_range=c(NULL, NULL), init_drive_x_range=c(NULL, NULL), init_drive_y_range=c(NULL, NULL),
                   drive_type=c("W_shredder", "X_shredder", "TADS")[1], conversion=1.0, fIntr=0.01, 
                   verbose=TRUE, plot.show=FALSE, vid.out=FALSE, vid.id=NULL, pop.out=FALSE){
  #################################################
  ### Define the landscape and intial positions ###
  #################################################
  if (y_min==y_max){
    one.dimensional = TRUE
  } else {
    one.dimensional = FALSE
  }
  ### initial positions of wild-types (WT)
  init_WT_x_min = init_WT_x_range[1]
  init_WT_x_max = init_WT_x_range[2]
  init_WT_y_min = init_WT_y_range[1]
  init_WT_y_max = init_WT_y_range[2]
  if (is.null(init_WT_x_min)){ init_WT_x_min = x_min }
  if (is.null(init_WT_y_min)){ init_WT_y_min = y_min }
  if (is.null(init_WT_x_max)){ init_WT_x_max = x_max }
  if (is.null(init_WT_y_max)){ init_WT_y_max = y_max }
  #################################
  ### Initialise the population ###
  #################################
  if (is.null(nInit)==TRUE){
    ### if nInit is not specified, then calculate the carrying capacity and use it as nInit
    if (one.dimensional==TRUE){
      nInit = ceiling(Nstar * (init_WT_x_max - init_WT_x_min))
    } else {
      nInit = ceiling(Nstar * (init_WT_x_max - init_WT_x_min) * (init_WT_y_max - init_WT_y_min))
    }
  }
  if (verbose) {
    cat("Initialising population of size", nInit, "\n")
  }
  mid_x_init_domain = (init_WT_x_max-init_WT_x_min)/2
  mid_y_init_domain = (init_WT_y_max-init_WT_y_min)/2
  pop = initInds(nInit, x_min=init_WT_x_min, y_min=init_WT_y_min, x_max=init_WT_x_max, y_max=init_WT_y_max)
  # pop.plot2D(pop, NULL, column="zState", nclass=20, pch=20, alpha=0.5, cex=1, gg, nGens, Rmax, Nstar, x_min, y_min, x_max, y_max, sigma, bw, torus, absorb, fIntr, nIntr=0, drive_type, conversion)
  #################################
  ### Initialise output vectors ###
  #################################
  n = rep(NA, nGens) ### population sizes
  m = rep(NA, nGens) ### number of males (ZZ)
  q = rep(NA, nGens) ### drive allele frequencies
  u = rep(NA, nGens) ### mean of the individual densities
  s = rep(NA, nGens) ### standard deviation of the individual densities
  k = rep(NA, nGens) ### number of re-introduction locations (re-introduction in the middle of optimally detected clusters)
  K = rep(NA, nGens) ### number of detected clusters
  drive_wave_x_min = rep(NA, nGens) ### lower tail of the drive wave along the x-axis
  drive_wave_x_max = rep(NA, nGens) ### upper tail of the drive wave along the x-axis
  drive_wave_y_min = rep(NA, nGens) ### lower tail of the drive wave along the y-axis
  drive_wave_y_max = rep(NA, nGens) ### upper tail of the drive wave along the y-axis
  drive_wave_x_width = rep(NA, nGens) ### width of the drive wave along the x-axis
  drive_wave_x_width_trailing = rep(NA, nGens) ### width of the drive wave along the x-axis
  drive_wave_y_width = rep(NA, nGens) ### width of the drive wave along the y-axis
  drive_wave_y_width_trailing = rep(NA, nGens) ### width of the drive wave along the y-axis
  drive_wave_height = rep(NA, nGens) ### height of the drive wave, i.e. the maximum drive-carrier density
  WT_wave_x_min = rep(NA, nGens) ### lower tail of the WT wave along the x-axis
  WT_wave_x_max = rep(NA, nGens) ### upper tail of the WT wave along the x-axis
  WT_wave_y_min = rep(NA, nGens) ### lower tail of the WT wave along the y-axis
  WT_wave_y_max = rep(NA, nGens) ### upper tail of the WT wave along the y-axis
  WT_wave_x_width = rep(NA, nGens) ### width of the WT wave along the x-axis
  WT_wave_y_width = rep(NA, nGens) ### width of the WT wave along the y-axis
  WT_wave_height = rep(NA, nGens) ### height of the WT wave, i.e. the maximum drive-carrier density
  penetration = rep(NA, nGens) ### logical: Did the drive wave got penetrated by the wild-types?
  edge_reached = FALSE ### initialise marker indicating when at least one of the landscape borders have been reached at which point wave property estimation stops to prevent miscalculations when movements become chaotic.
  ####################
  ### Burn-in step ###
  #################### Burn-in within the defined intial spatial domain assuming it is bounded with absorbing boundaries to generate a WT wave.
  if (verbose) {
    cat("Burning in for", burnIn, "generations... \n")
  }
  for (bb in rep(1,burnIn)){
    pop = repro(pop, NULL, Rmax, Nstar, sigma, bw, drive_type, torus=FALSE, absorb=TRUE, conversion=conversion, x_min=init_WT_x_min, y_min=init_WT_y_min, x_max=init_WT_x_max, y_max=init_WT_y_max)
    if (verbose) print(paste0("Gen=", bb, "; N=", nrow(pop)))
  }
  # pop.plot2D(pop, NULL, column="zState", nclass=20, pch=20, alpha=0.5, cex=1, gg, nGens, Rmax, Nstar, x_min, y_min, x_max, y_max, sigma, bw, torus, absorb, fIntr, nIntr=0, drive_type, conversion)
  ################################################################
  ### Intial population size and drive loss generation counter ###
  ################################################################
  N0 = nrow(pop)
  drive_loss_counter = 0
  ###########################
  ### Introduce the drive ###
  ########################### Introduce the drive within the defined initial spatial domain.
  ### calculate the number of gene drive heterozygote individuals we need to introduce
  nIntr = ceiling(fIntr * N0) ### equivalent to ceiling(fIntr * nrow(pop))
  if (verbose) {
    print(paste0("Introducing ", nIntr, " (", round(fIntr*100), "%)", " drive individuals."))
  }
  ### initial positions of suppression gene drive-carriers
  init_drive_x_min = init_drive_x_range[1]
  init_drive_x_max = init_drive_x_range[2]
  init_drive_y_min = init_drive_y_range[1]
  init_drive_y_max = init_drive_y_range[2]
  if (is.null(init_drive_x_min)){ init_drive_x_min = min(pop[,"X"]) }
  if (is.null(init_drive_y_min)){ init_drive_y_min = min(pop[,"Y"]) }
  if (is.null(init_drive_x_max)){ init_drive_x_max = min(pop[,"X"]) }
  if (is.null(init_drive_y_max)){ init_drive_y_max = min(pop[,"Y"]) }
  pop =  rbind(pop, introDrive(nIntr, drive_type, sigma, torus=FALSE, absorb=FALSE,
                               x_min=init_drive_x_min,
                               y_min=init_drive_y_min,
                               x_max=init_drive_x_max,
                               y_max=init_drive_y_max))
  # pop.plot2D(pop, NULL, column="zState", nclass=20, pch=20, alpha=0.5, cex=1, gg, nGens, Rmax, Nstar, x_min, y_min, x_max, y_max, sigma, bw, torus, absorb, fIntr, nIntr, drive_type, conversion)
  ##################################
  ###                            ###
  ### ITERATE ACROSS GENERATIONS ###
  ###                            ###
  ##################################
  for (gg in 1:nGens){
    #####################################
    ### Has the population collapsed? ###
    #####################################
    not_extinct = !is.null(pop) & !is.null(nrow(pop)) & length(nrow(pop))!=0
    if (not_extinct){
      ###################################################################
      ### Extract summary statistics if the population is not extinct ###
      ###################################################################
      vec_density = metrics(pop, bw, torus, x_min, y_min, x_max, y_max)
      n[gg] = nrow(pop)
      m[gg] = sum(pop[,"Z0"]==1)
      q[gg] = sum(pop[,"zState"])/(2*n[gg])
      u[gg] = mean(vec_density, na.rm=TRUE)
      s[gg] = sd(vec_density, na.rm=TRUE)
      ################################
      ### Invasion wave properties ###
      ################################
      wave_prop = wave_properties(pop, sigma, bw, torus, x_min, y_min, x_max, y_max, edge_reached, Nstar)
      drive_wave_x_min[gg]   = wave_prop$drive_wave_x_min
      drive_wave_x_max[gg]   = wave_prop$drive_wave_x_max
      drive_wave_y_min[gg]   = wave_prop$drive_wave_y_min
      drive_wave_y_max[gg]   = wave_prop$drive_wave_y_max
      if (gg > 1){
        ### start saving the widths and heights of the drive wave after 1 generation of reproduction not after just dispersion
        drive_wave_x_width[gg] = wave_prop$drive_wave_x_width
        drive_wave_x_width_trailing[gg] = wave_prop$drive_wave_x_width_trailing
        drive_wave_y_width[gg] = wave_prop$drive_wave_y_width
        drive_wave_y_width_trailing[gg] = wave_prop$drive_wave_y_width_trailing
        drive_wave_height[gg]  = wave_prop$drive_wave_height
      }
      WT_wave_x_min[gg]   = wave_prop$WT_wave_x_min
      WT_wave_x_max[gg]   = wave_prop$WT_wave_x_max
      WT_wave_y_min[gg]   = wave_prop$WT_wave_y_min
      WT_wave_y_max[gg]   = wave_prop$WT_wave_y_max
      WT_wave_x_width[gg] = wave_prop$WT_wave_x_width
      WT_wave_y_width[gg] = wave_prop$WT_wave_y_width
      WT_wave_height[gg]  = wave_prop$WT_wave_height
      penetration[gg]  = wave_prop$penetration
      edge_reached = wave_prop$edge_reached
      ######################
      ### PRINT PROGRESS ###
      ######################
      if (verbose) {
        ### print population size
        cat("Gen = ", gg, ": Pop size =", nrow(pop), "\n")
      }
      ############
      ### Plot ###
      ############
      if (vid.out | plot.show) {
        if (is.null(vid.id)){
          vid.id=round(abs(runif(1)*10000))
        }
        ### save the plot if we want to output a video of the chase
        leading_zeros = nchar(as.character(nGens)) - nchar(as.character(gg))
        gen = paste0(paste(rep(0, times=leading_zeros),collapse=""), gg)
        png(paste0("drive_chase-id_", vid.id, "-gen_", gen, ".png"), width=1200, height=1000)
      }
      if (vid.out | plot.show) {
        ### plot the distribution of the genotypes across the landscape
        pop.plot2D(pop, vec_density, nclass=50, cex=3, gg=gg, nGens=nGens, Rmax=Rmax, Nstar=Nstar, x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max, sigma=sigma, bw=bw, torus=torus, absorb=absorb, fIntr=fIntr, nIntr=nIntr, drive_type=drive_type, conversion=conversion, wave_prop=wave_prop)
      }
      if (vid.out | plot.show) {
        ### save the plot if we want to output a video of the chase
        dev.off()
      }
      ###############################################
      ### Find out if we've lost the drive allele ###
      ###############################################
      driveLoss = sum(pop[,"zState"])==0
      ### If the drive is lost, then terminate after the population has reached equillibrium size
      if (driveLoss) {
        if (verbose) {
          cat("Drive extinct at", gg, "gens\n")
        }
        ### 10 generations after losing the drive allele
        ### break from the for-loop iterating across generations
        if (drive_loss_counter < 10) {
          drive_loss_counter = drive_loss_counter + 1
        } else {
          break
        }
      }
    } else {
      ######################################################
      ### Stop this malarky if the population is extinct ###
      ######################################################
        if (verbose) {
          cat("Population extinct at", gg, "generations\n")
        }
        ### population size of zero at extinction
        n[gg] = 0
        break
    }
    ###################################################################################################
    ### Simulate reproduction: find mates, generate offsprings, gene drive, and disperse offspring  ###
    ### (gene drive may come before generating offspring, i.e. pre-zygotic, or after, post-zygotic) ###
    ###################################################################################################
    pop = repro(pop, vec_density, Rmax, Nstar, sigma, bw, drive_type, torus, absorb, conversion, x_min, y_min, x_max, y_max)
  }
  #########################################
  ###                                   ###
  ### END ITERATIONS ACROSS GENERATIONS ###
  ###                                   ###
  #########################################

  ####################################
  ### GENERATE OUTPUT AND CLEAN-UP ###
  ####################################
  ### generate video from the scatterplots and cleanup, i.e. remove the plot images
  if (vid.out){
    library(av)
    av_encode_video(system(paste0("ls drive_chase-id_", vid.id, "-gen_*.png"), intern=TRUE), framerate=10, output=paste0("drive_chase-id_", vid.id, ".mp4"))
    system(paste0("rm drive_chase-id_", vid.id, "-gen_*.png"))
    # av_encode_video(system(paste0("ls filled*-pop-gen-*.png"), intern=TRUE), framerate=10, output=paste0("filledcontour_pop_", vid.id, ".mp4"))
    # system(paste0("rm filled*-pop-gen-*.png"))
    # av_encode_video(system(paste0("ls filled*-drive-gen-*.png"), intern=TRUE), framerate=10, output=paste0("filledcontour_drive_", vid.id, ".mp4"))
    # system(paste0("rm filled*-drive-gen-*.png"))
  }
  if (verbose) {
    cat("Simulation complete.\n")
  }
  ### output list of population sizes, drive allele frequencies, means and standard deviations of the density distributions across generations, etc..., 
  OUT = list(n=n, m=m, q=q, u=u, s=s, k=k, 
             drive_wave_x_min=drive_wave_x_min,
             drive_wave_x_max=drive_wave_x_max,
             drive_wave_y_min=drive_wave_y_min,
             drive_wave_y_max=drive_wave_y_max,
             drive_wave_x_width=drive_wave_x_width,
             drive_wave_x_width_trailing=drive_wave_x_width_trailing,
             drive_wave_y_width=drive_wave_y_width,
             drive_wave_y_width_trailing=drive_wave_y_width_trailing,
             drive_wave_height=drive_wave_height,
             WT_wave_x_min = WT_wave_x_min,
             WT_wave_x_max = WT_wave_x_max,
             WT_wave_y_min = WT_wave_y_min,
             WT_wave_y_max = WT_wave_y_max,
             WT_wave_x_width = WT_wave_x_width,
             WT_wave_y_width = WT_wave_y_width,
             WT_wave_height = WT_wave_height,
             penetration=penetration)
  if (pop.out){
    OUT$pop = pop
  }
  return(OUT)
}

###################################
### SUMMARISE SIMULATION OUTPUT ###
###################################
### summarise the output dataframe of evolve(); plot some time-series graphs and write-out csv files
### NOTE: for 1D simulations set y_min == y_max
#################################################################################################################################################
### TEST:
# nGens=100; Rmax=4; Nstar=5; sigma=5 # dispersal_mean = sigma * sqrt(2/pi)
# # bw=1; torus=TRUE; absorb=TRUE; x_min=0; y_min=0; x_max=50; y_max=50
# bw=1; torus=TRUE; absorb=TRUE; x_min=0; y_min=0; x_max=1000; y_max=0
# nInit=NULL; burnIn=0
# init_WT_x_range=c((x_max*0.2), (x_max*0.6))
# init_WT_y_range=c((y_max*0.2), (y_max*0.6))
# init_drive_x_range=c(NULL, NULL)
# init_drive_y_range=c(NULL, NULL)
# ### suppression gene drive parameters
# drive_type=c("W_shredder", "X_shredder", "TADS")[1]
# conversion=1.0; fIntr=0.01
# ### output paramters
# verbose=TRUE
# plot.show=TRUE
# vid.out=TRUE
# vid.id="test"
# pop.out=FALSE
# set.seed(123)
# output_dir="./"; id=vid.id; rep=1
# landscape = "Torus"
# plot.time.series=TRUE
# plot.velocities=TRUE
# landscape="Torus"
# test.out = evolve(nGens, Rmax, Nstar, sigma, bw, torus, absorb, x_min, y_min, x_max, y_max, nInit, burnIn, init_WT_x_range, init_WT_y_range, init_drive_x_range, init_drive_y_range, drive_type, conversion, fIntr, verbose, plot.show, vid.out, vid.id, pop.out)
# test = extract.metrics.and.plot(test.out, output_dir, id, rep, plot.time.series, plot.velocities, vid.out, Rmax, Nstar, sigma, bw, drive_type, conversion, landscape)
#################################################################################################################################################
extract.metrics.and.plot = function(out, output_dir="./", id=1, rep=1, plot.time.series=TRUE, plot.velocities=TRUE, vid.out=FALSE, Rmax=3, Nstar=5, sigma=1, bw=1, drive_type=c("W_shredder", "X_shredder", "TADS")[1], conversion=1, landscape="Torus"){
  ########################
  ### Plot time-series ###
  ########################
  if (plot.time.series){
    png(paste0(output_dir, "drive_chase-timeseries-id-", id, ".png"), width=1200, height=600)
    vec_colours = c("black", "green", "red", "black", "orange", "purple")
    par(mfrow=c(2,1), mar=c(5,5,2,5))
    ### (1/2) plot population size, no. of males, and drive allele frquency fluctuations
    n_max = max(c(out$n,out$m),na.rm=TRUE)
    if (n_max>1000){
      nth = length(strsplit(as.character(n_max), "")[[1]])-1
      nth_factor = as.numeric(paste0(c("1", rep("0", times=nth)), collapse=""))
      ylab = paste0("Number of individuals (x", nth_factor, ")")
    } else {
      ylab = "Number of individuals"
    }
    plot(out$n, xlim=c(0, nGens), ylim=c(0, max(c(out$n,out$m),na.rm=TRUE)), type="l", col=vec_colours[1], xlab="Generation", ylab=ylab, yaxt="n")
    axis(side=2, las=2, at=seq(0, n_max, length=5), labels=formatC(seq(0, n_max/nth_factor, length=5), format="f", flag="0", digits=2))
    lines(out$m, xlim=c(0, nGens), col=vec_colours[2])
    par(new=TRUE)
    plot(out$q, xlim=c(0, nGens), ylim=c(0, max(out$q,na.rm=TRUE)), type="l", col=vec_colours[3], xaxt="n", yaxt="n", xlab="", ylab="")
    axis(side=4, las=2, at=seq(min(out$q,na.rm=TRUE),max(out$q,na.rm=TRUE),length=4), labels=formatC(seq(min(out$q,na.rm=TRUE),max(out$q,na.rm=TRUE),length=4), format="f", flag="0", digits=2))
    mtext(side=4, text="Gene drive allele frequency (q)", padj=5)
    grid()
    legend("bottomleft", legend=c("Population size", "Number of Males", "Drive allele frequency"), col=vec_colours[1:3], lwd=2, bty="n")
    ### (/2) plot density distribution paramters excluding the intial generation which are often too high which reduces plot resolution of smaller variations for the rest of the generations
    d_min = min(tail(out$u,-1),na.rm=TRUE)
    d_max = max(tail(out$u,-1),na.rm=TRUE)
    plot(out$u, xlim=c(0, nGens), ylim=c(d_min, d_max), type="l", col=vec_colours[4], xlab="Generation", ylab="Mean density", yaxt="n")
    axis(side=2, las=2, at=seq(d_min, d_max, length=5), labels=formatC(seq(d_min, d_max, length=5), format="f", flag="0", digits=2))
    par(new=TRUE)
    s_min = min(tail(out$s,-1),na.rm=TRUE)
    s_max = max(tail(out$s,-1),na.rm=TRUE)
    plot(out$s, xlim=c(0, nGens), ylim=c(s_min, s_max), type="l", col=vec_colours[5], xaxt="n", yaxt="n", xlab="", ylab="")
    axis(side=4, las=2, at=seq(s_min, s_max, length=4), labels=formatC(seq(s_min, s_max, length=4), format="f", flag="0", digits=2))
    mtext(side=4, text="Standard deviation of density", padj=5)
    grid()
    legend("bottomleft", legend=c("Mean of density", "Standard deviation of density"), col=vec_colours[4:5], lwd=2, bty="n")
    dev.off()
  }
  ############################
  ### Rename the video out ###
  ############################
  if (vid.out){
    system(paste0("mv drive_chase-id_", id, ".mp4 ", output_dir, "drive_chase-scatterplot_density-id-", id, ".mp4"))
  }
  ######################################
  ### Prepare 1-row output dataframe ###
  ######################################
  df_output = data.frame(id=id,
                         rep=rep,
                         Rmax=Rmax,
                         Nstar=Nstar,
                         sigma=sigma,
                         bw=bw,
                         drive_type=drive_type,
                         conversion=conversion,
                         introduced=fIntr*(out$n[1]/(fIntr+1)),
                         landscape=landscape)
  ###############################
  ### Compute wave velocities ###
  ###############################
  ### i.e. - drive_velocity_x
  ###      - drive_velocity_y
  ###      - WT_velocity_x
  ###      - WT_velocity_y
  if (plot.velocities) {
    png(paste0(output_dir, "drive_chase-wave_velocities-id-", id, ".png"), width=900, height=600)
    par(mfrow=c(2,2))
  }
  vec_genotypes = c("drive", "WT")
  for (genotype in vec_genotypes){
    for (spatial_axis in c("x", "y")){
      # genotype = "WT"
      # spatial_axis = "x"
      vec_gen = c(1:length(out$n))
      vec_pos = eval(parse(text=paste0("out$", genotype, "_wave_", spatial_axis, "_max")))
      ### restrict right-most position once it's reached
      epsilon = 10
      idx = (eval(parse(text=paste0(spatial_axis, "_max"))) - vec_pos) <= epsilon
      if (sum(idx, na.rm=TRUE)>0){
        gen_edge_reached = min(c(1:length(idx))[idx], na.rm=TRUE)
        vec_pos[gen_edge_reached:length(vec_pos)] = NA
      }
      if (sum(!is.na(vec_pos))>=2){
        idx = !is.infinite(vec_pos)
        mod = lm(vec_pos[idx] ~ vec_gen[idx]) ### or manually: idx = !is.infinite(vec_pos) & !is.na(vec_pos); X = cbind(rep(1,sum(idx)), vec_gen[idx]); y = vec_pos[idx]; solve(t(X)%*%X) %*% t(X) %*% y
        velocity = mod$coefficients[2]
        if (plot.velocities) {plot(vec_gen, vec_pos, xlim=c(0,length(vec_gen)), ylim=c(0,eval(parse(text=paste0(spatial_axis, "_max")))), main=paste0(genotype, " (", spatial_axis, "-axis)")); abline(mod, col="red"); grid(); legend("bottomright", legend=paste0("velocity=", round(velocity,2)))}
      } else {
        velocity = NA
        if (plot.velocities) {plot(0, xlim=c(0,length(vec_gen)), ylim=c(0,eval(parse(text=paste0(spatial_axis, "_max")))), main=paste0(genotype, " (", spatial_axis, "-axis)"), type="n"); grid(); legend("bottomright", legend=paste0("velocity=", round(velocity,2)))}
      }
      eval(parse(text=paste0("df_output$", genotype, "_velocity_", spatial_axis, " = velocity")))
    }
  }
  if (plot.velocities) {
    dev.off()
  }
  ##############################################################
  ### Compute relative drive wave velocity, height and width ###
  ##############################################################
  ### relative drive velocities (drive / WT)
  for (spatial_axis in c("x", "y")){
    # spatial_axis = "x"
    eval(parse(text=paste0("df_output$relative_drive_velocity_", spatial_axis, " = df_output$drive_velocity_", spatial_axis, " / df_output$WT_velocity_", spatial_axis)))
  }
  ### relative drive height (drive / WT) on the out list (summary stats will be appended into df_output below)
  out$relative_drive_wave_height = out$drive_wave_height / out$WT_wave_height
  ### relative drive width (drive / WT) on the out list (summary stats will be appended into df_output below)
  for (spatial_axis in c("x", "y")){
    # spatial_axis = "x"
    eval(parse(text=paste0("out$relative_drive_wave_", spatial_axis, "_width = out$drive_wave_", spatial_axis, "_width / out$WT_wave_", spatial_axis, "_width")))
  }
  ###############################################################################
  ### Summary statistics of the absolute and relative wave heights and widths ###
  ###############################################################################
  ### define the summary statistics calculator function
  sumstats <- function(x, no.neg=FALSE){
    x = x[!is.na(x) & !is.infinite(x)] ### remove NAs and infinities
    if (no.neg){
      x = x[!(x<0)] ### remove negatives from density calculations resulting from C data leak
    }
    if (length(x)>2){
      init = head(x, 1)
      fin = tail(x, 1)
      x = head(tail(x,-1),-1) ### remove possible outliers at the first and final generations
      mu  = mean(x)
      sd  =   sd(x)
      min =  min(x)
      max =  max(x)
      out = cbind(mu, sd, min, max, init, fin)
    } else {
      out = cbind(mu=NA, sd=NA, min=NA, max=NA, init=NA, fin=NA)
    }
    return(out)
  }
  ### compute summary statistics
  vec_metrics = c("n", "q", "m",
                  "drive_wave_x_max", "drive_wave_x_min", "drive_wave_y_max", "drive_wave_y_min",
                  "WT_wave_x_max", "WT_wave_x_min", "WT_wave_y_max", "WT_wave_y_min",
                  "drive_wave_height", "drive_wave_x_width", "drive_wave_y_width", "drive_wave_x_width_trailing", "drive_wave_y_width_trailing",
                  "WT_wave_height", "WT_wave_x_width", "WT_wave_y_width",
                  "relative_drive_wave_height", "relative_drive_wave_x_width", "relative_drive_wave_y_width")
  for (metric in vec_metrics){
    # metric = "n"
    eval(parse(text=paste0("df_output = cbind(df_output, sumstats(out$", metric, "))")))
    eval(parse(text=paste0("colnames(df_output)[(ncol(df_output)-5): ncol(df_output)] = paste0('", metric, "_', colnames(df_output)[(ncol(df_output)-5): ncol(df_output)])")))
  }
  #################################
  ### Did the population crash? ###
  #################################
  vec_n = out$n[!is.na(out$n)]
  if (tail(vec_n, 1) == 0 ){
    df_output$crash = TRUE
  } else {
    df_output$crash = FALSE
  }
  #############################################################################
  ### When did the wild-types able to penetrate the wave of drive-carriers? ###
  #############################################################################
  vec_penetration = out$penetration[!is.na(out$penetration)]
  if (sum(vec_penetration)>0){
    df_output$penetration_gen = min(c(1:length(vec_penetration))[vec_penetration])
  } else {
    df_output$penetration_time = NA
  }
  ##############
  ### Output ###
  ##############
  fname_out = paste0(output_dir, "drive_chase-output-id.", id, ".csv")
  write.table(df_output, file=fname_out, sep=",", col.names=TRUE, row.names=FALSE)
  return(df_output)
}
