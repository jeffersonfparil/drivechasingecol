##################################################################################################
###                                                                                            ###
### PDEs for W-shredder, X-shredder, and TADS suppression gene drives in one spatial dimension ###
###                                                                                            ###
##################################################################################################

##########################################
### Load packages for solving the PDEs ###
##########################################
library(deSolve)

#############################################################################
### Partial differential equation of population size with respect to time ###
#############################################################################
func_pdn_pdt = function(t, n, parameters){
  ### input
  # t = current time-point (not used in this function per se but required by deSolve::ode.1D() function)****
  # n  = population size along the landscape (for each interval across the spatial x-axis)
  # parameters = list(q, D, r, k, x_min, xmax)
  # where: q = suppression drive allele frequency
  #        D = diffusion coefficient
  #        r = intrinsic growth rate (unit: n per n)
  #        K = maximum population size in each interval, i.e. carrying capacity
  #        dx = step size along the **spatial** x-axis
  #        dType = drive type; select from "W_shredder" or "TADS"
  #        c = conversion efficiency for W-shredder, i.e. the efficiency of the gene drive allele to perform its task [NOTE: if dType is not among the ones listed above, this parameter is treated as the fitness, w]
  #        k = P(TTdd) for TADS
  #        f = P(Ttdd) for TADS
  #        h = P(ttdd) for TADS
  with(as.list(c(n, parameters)), {
    n = c(head(n,1), n, tail(n,1))
    dn_dx = ( tail(n, -1) - head(n, -1) ) / dx
    d2n_dx2 = ( tail(dn_dx, -1) - head(dn_dx, -1) ) / dx
    n = head(tail(n, -1), -1)
    dispersion_rate        = D*d2n_dx2      ### equilibriates between n_max and n_min
    growth_rate            = r*n*(1-(n/K)) ### logistic growth - the continuous analogue of Beverton-Holt growth rate
    if ((dType=="W_shredder") | (dType=="X_shredder")){
      P_female = (2-c-c*q)/(4-2*c)
      drive_suppression_rate = 2 * n * r * ((P_female/0.5) - 1)
      # drive_suppression_rate = -2 * r * n * (0.5 - P_female)
    } else if (dType=="TADS"){
      P_non_productive_matings = k + f + h
      drive_suppression_rate = 2 * n * r * (-1 * P_non_productive_matings)
    } else {
      ### if dType is not among the ones listed above we use this generic growth rate formula with c as the fitness, w
      w = c
      drive_suppression_rate = -r*q*w
    }
    pdn_pdt = dispersion_rate + growth_rate + drive_suppression_rate
    list(pdn_pdt)
  })
}

#####################################################################################################
### Partial differential equation of suppression gene drive allele frequency with respect to time ###
#####################################################################################################
func_pdq_pdt = function(t, q, parameters){
  ### input
  # t = current time-point (not used in this function per se but required by deSolve::ode.1D() function)
  # q = suppression drive allele frequency
  # parameters = list(D, dx, dType, c, n, ...)
  # where: D = diffusion coefficient
  #        dx = step size along the spatial x-axis
  #        dType = drive type; select from "W_shredder" or "TADS"
  #        c = conversion efficiency for W-shredder, i.e. the efficiency of the gene drive allele to perform its task [NOTE: if dType is not among the ones listed above, this parameter is treated as the fitness, w]
  #        n = population size along the landscape (for each interval across the spatial x-axis)
  #        b = P(TTDd) for TADS
  #        k = P(TTdd) for TADS
  #        e = P(TtDd) for TADS
  #        f = P(Ttdd) for TADS
  #        g = P(ttDd) for TADS
  #        h = P(ttdd) for TADS
  #        v = conversion efficiency in males for TADS
  with(as.list(c(q, parameters)), {
    q = c(head(q,1), q, tail(q,1))
    dq_dx = ( head(q, -1) - tail(q, -1) ) / dx
    d2q_dx2 = ( head(dq_dx, -1) - tail(dq_dx, -1) ) / dx
    q = head(tail(q, -1), -1)
    dispersion_rate = D*d2q_dx2 ### equilibriates between q_max and q_min
    if (dType=="W_shredder"){
      drive_allele_growth_rate = (c/(8-4*c)) * q*(1-q) ### dq_dt
      # q_t = ( q_0*exp(c*t/(8-4*c)) ) / ( 1+(q_0*((exp(c*t/(8-4*c)))-1)) ) ### solution to the ode above
    } else if (dType=="X_shredder"){
      drive_allele_growth_rate = (c/(2-c+c*q)) * q*(1-q) ### dq_dt
      # q_t = ( q_0*exp(c*t/(8-4*c)) ) / ( 1+(q_0*((exp(c*t/(8-4*c)))-1)) ) ### solution to the ode above
    } else if (dType=="TADS"){
      drive_allele_growth_rate = (1/4) * ((((2*b*(3-v))-(2*(2-v)*(-2*e-g*(3-v))))/((3-v)*(2-v)*(1-k-f-h))) - (b + 2*k + e + 2*f + g + 2*h))
    } else {
      ### if dType is not among the ones listed above we use this generic growth rate formula with c as the fitness, w
      w = c
      drive_allele_growth_rate = w*q*(1-q) ### caps at q=1 or q=0
    }
    ### advection term
    logn = log(n)
    logn = c(head(logn,1), logn, tail(logn,1))
    dlogn_dx = ( tail(logn, -1) - head(logn, -1) ) / dx
    d2logn_dx2 = ( tail(dlogn_dx, -1) - head(dlogn_dx, -1) ) / dx
    advection_term = 2 * d2logn_dx2 * d2q_dx2
    ### overall rate across space and time    
    pdq_pdt = dispersion_rate + advection_term + drive_allele_growth_rate
    # pdq_pdt = dispersion_rate + drive_allele_growth_rate
    list(pdq_pdt)
  })
}

func_pdu_pdt = function(t, u, parameters){
  with(as.list(c(q, parameters)), {
    u = c(head(u,1), u, tail(u,1))
    du_dx = ( head(u, -1) - tail(u, -1) ) / dx
    d2u_dx2 = ( head(du_dx, -1) - tail(du_dx, -1) ) / dx
    u = head(tail(u, -1), -1)
    dispersion_rate = D*d2u_dx2
    if (dType=="W_shredder"){
      drive_allele_growth_rate = (c/(8-4*c)) * q*(1-q) ### dq_dt
      # q_t = ( q_0*exp(c*t/(8-4*c)) ) / ( 1+(q_0*((exp(c*t/(8-4*c)))-1)) ) ### solution to the ode above
    } else if (dType=="X_shredder"){
      drive_allele_growth_rate = (c/(2-c+c*q)) * q*(1-q) ### dq_dt
      # q_t = ( q_0*exp(c*t/(8-4*c)) ) / ( 1+(q_0*((exp(c*t/(8-4*c)))-1)) ) ### solution to the ode above
    } else if (dType=="TADS"){
      drive_allele_growth_rate = (1/4) * ((((2*b*(3-v))-(2*(2-v)*(-2*e-g*(3-v))))/((3-v)*(2-v)*(1-k-f-h))) - (b + 2*k + e + 2*f + g + 2*h))
    } else {
      ### if dType is not among the ones listed above we use this generic growth rate formula with c as the fitness, w
      w = c
      drive_allele_growth_rate = w*q*(1-q) ### caps at q=1 or q=0
    }
    ### pdn_pdt
    n = c(head(n,1), n, tail(n,1))
    dn_dx = ( tail(n, -1) - head(n, -1) ) / dx
    d2n_dx2 = ( tail(dn_dx, -1) - head(dn_dx, -1) ) / dx
    n = head(tail(n, -1), -1)
    dispersion_rate_n        = D*d2n_dx2      ### equilibriates between n_max and n_min
    growth_rate            = r*n*(1-(n/K)) ### logistic growth - the continuous analogue of Beverton-Holt growth rate
    if ((dType=="W_shredder") | (dType=="X_shredder")){
      P_female = (2-c-c*q)/(4-2*c)
      drive_suppression_rate = 2 * n * r * ((P_female/0.5) - 1)
      # drive_suppression_rate = -2 * r * n * (0.5 - P_female)
    } else if (dType=="TADS"){
      P_non_productive_matings = k + f + h
      drive_suppression_rate = 2 * n * r * (-1 * P_non_productive_matings)
    } else {
      ### if dType is not among the ones listed above we use this generic growth rate formula with c as the fitness, w
      w = c
      drive_suppression_rate = -r*q*w
    }
    pdn_pdt = dispersion_rate_n + growth_rate + drive_suppression_rate
    pdn_pdt_TIMES_dq_dt = pdn_pdt * drive_allele_growth_rate
    CORRECTION_TERM = (n*drive_allele_growth_rate) + (q*pdn_pdt)
    pdu_pdt = dispersion_rate + pdn_pdt_TIMES_dq_dt - CORRECTION_TERM
    list(pdu_pdt)
  })
}

###########################################
### Mechanistic TADS recursive function ###
###########################################
func_next_gen_freqs = function(list_geno_freqs, list_conversions, xqz_freqs=FALSE){
  ### test:
  # list_geno_freqs=list(a=0.99, b=0.01, k=0, d=0, e=0, f=, g=0, h=0)
  # list_geno_freqs=list(a=0.000165674459515909, b=0.000194773913310173, c=5.23754551912644e-05, d=0.00284877854123094, e=0.0500106011884987, f=0.0213042280150589, g=0.0489227695209421, h=0.789614024576454, x=0.875618597969791, q=0.86053470035808, z=0.128867601555645)
  # list_conversions=list(c_f=0.90, c_m=0.95)
  #########
  ### initial genotype frequencies (vectors)
  vec_a = list_geno_freqs$a
  vec_b = list_geno_freqs$b
  vec_k = list_geno_freqs$k
  vec_d = list_geno_freqs$d
  vec_e = list_geno_freqs$e
  vec_f = list_geno_freqs$f
  vec_g = list_geno_freqs$g
  vec_h = list_geno_freqs$h
  ### conversion rates in females and males
  u = list_conversions$c_f
  v = list_conversions$c_m
  ### prepare output vectors
  vec_a1 = c()
  vec_b1 = c()
  vec_k1 = c()
  vec_d1 = c()
  vec_e1 = c()
  vec_f1 = c()
  vec_g1 = c()
  vec_h1 = c()
  vec_x1 = c()
  vec_q1 = c()
  vec_z1 = c()
  for (i in 1:length(vec_a)){
    a = vec_a[i]
    b = vec_b[i]
    k = vec_k[i]
    d = vec_d[i]
    e = vec_e[i]
    f = vec_f[i]
    g = vec_g[i]
    h = vec_h[i]
    ### vectors of genotype frequencies in females and males (P(female)==P(male)==0.5)
    vec_female_freqs = matrix(c(a, b, k, d, e, f, g, h), ncol=1)
    vec_male_freqs   = vec_female_freqs
    rownames(vec_female_freqs) = c("a_f", "b_f", "k_f", "d_f", "e_f", "f_f", "g_f", "h_f")
    rownames(vec_male_freqs) =   c("a_m", "b_m", "k_m", "d_m", "e_m", "f_m", "g_m", "h_m")
    ### matrices of gamete frequencies per genotype in females and males
    mat_female_gametes = matrix(c(      1,       0,       0,       0,
                                  (1-u)/2, (1-u)/2,     u/2,     u/2,
                                        0,     1-u,       0,       u,
                                      1/2,       0,     1/2,       0,
                                  (1-u)/4, (1-u)/4, (1+u)/4, (1+u)/4,
                                        0, (1-u)/2,       0, (1+u)/2,
                                        0,       0,     1/2,     1/2,
                                        0,       0,       0,       1), byrow=TRUE, nrow=8)
    mat_male_gametes =   matrix(c(          1,           0,       0,           0,
                                  (1-v)/(2-v), (1-v)/(2-v),       0,     v/(2-v),
                                            0,           0,       0,           0,
                                            1,           0,       0,           0,
                                  (1-v)/(3-v), (1-v)/(3-v),       0, (1+v)/(3-v),
                                            0,           0,       0,           0,
                                            0,           0,       0,           1,
                                            0,           0,       0,           0), byrow=TRUE, nrow=8)
    rownames(mat_female_gametes) = c("a_f", "b_f", "k_f", "d_f", "e_f", "f_f", "g_f", "h_f")
    rownames(mat_male_gametes) =   c("a_m", "b_m", "k_m", "d_m", "e_m", "f_m", "g_m", "h_m")
    colnames(mat_female_gametes) = c("A", "B", "C", "D")
    colnames(mat_male_gametes) =   c("A", "B", "C", "D")
    ### matrix of offspring frequencies of the current genotypes  
    mat_offspring_freqs = matrix(0, nrow=4, ncol=4)
    for (j in 1:nrow(mat_female_gametes)){
      for (l in 1:nrow(mat_male_gametes)){
        ### test: j=1; l=2
        vec_female_gametes = mat_female_gametes[j,]
        vec_male_gametes =   mat_male_gametes[l,]
        mat_offsprings = (vec_female_gametes*vec_female_freqs[j]) %*% t(vec_male_gametes*vec_male_freqs[l])
        mat_offspring_freqs = mat_offspring_freqs + mat_offsprings
      }
    }
    ### frequency of productive matings (i.e. not involving sterile males: TTdd, Ttdd, and ttdd)
    z1 = sum(mat_offspring_freqs)
    ### normalise the offspring frequency using the frequency of productive matings
    mat_offspring_freqs = mat_offspring_freqs/z1
    ### genotype frequencies in the next generation
    a1 = mat_offspring_freqs[1,1]
    b1 = mat_offspring_freqs[1,2] + mat_offspring_freqs[2,1]
    k1 = mat_offspring_freqs[2,2]
    d1 = mat_offspring_freqs[3,1]
    e1 = mat_offspring_freqs[1,4] + mat_offspring_freqs[3,2] + mat_offspring_freqs[4,1]
    f1 = mat_offspring_freqs[2,4] + mat_offspring_freqs[4,2]
    g1 = mat_offspring_freqs[3,4]
    h1 = mat_offspring_freqs[4,4]
    ### frequency of the non-functional allele in the target locus in the next generations
    x1 = ((d1+e1+f1)/2) + g1 + h1
    ### frequency of the suppression gene drive allele in the next generations
    q1 = ((b1+e1+g1)/2) + k1 + f1 + h1
    ### append genotype frequencies in the next generation into the output vectors
    vec_a1 = c(vec_a1, a1)
    vec_b1 = c(vec_b1, b1)
    vec_k1 = c(vec_k1, k1)
    vec_d1 = c(vec_d1, d1)
    vec_e1 = c(vec_e1, e1)
    vec_f1 = c(vec_f1, f1)
    vec_g1 = c(vec_g1, g1)
    vec_h1 = c(vec_h1, h1)
    vec_x1 = c(vec_x1, x1)
    vec_q1 = c(vec_q1, q1)
    vec_z1 = c(vec_z1, z1)
  }
  ### output
  if (xqz_freqs){
    out = list(a=vec_a1, b=vec_b1, k=vec_k1, d=vec_d1, e=vec_e1, f=vec_f1, g=vec_g1, h=vec_h1, x=vec_x1, q=vec_q1, z=vec_z1)
  } else {
    out = list(a=vec_a1, b=vec_b1, k=vec_k1, d=vec_d1, e=vec_e1, f=vec_f1, g=vec_g1, h=vec_h1)
  }
  return(out)
}

##########################################################################
### Partial differential equation for the dispersion of TADS genotypes ###
##########################################################################
func_pdF_pdt = function(t, F, parameters){
  ### input
  # t = current time-point (not used in this function per se but required by deSolve::ode.1D() function)
  # F = TADS genotype frequencies across space
  # where: D = diffusion coefficient
  #        dx = step size along the spatial x-axis
  #        w = change in F from the previous to the current generation (an approximation of the fitness of the genotype)
  with(as.list(c(F, parameters)), {
    F = c(head(F,1), F, tail(F,1))
    dF_dx = ( head(F, -1) - tail(F, -1) ) / dx
    d2F_dx2 = ( head(dF_dx, -1) - tail(dF_dx, -1) ) / dx
    F = head(tail(F, -1), -1)
    dispersion_rate = D*d2F_dx2 + (w*F*(1-F))
    list(dispersion_rate)
  })
}

###############################
### Disperse TADS genotypes ###
###############################
func_disperse_TADS_genotypes = function(F, w, dx, D=1.00){
  F_parameters = list(D=D, w=w, dx=dx)
  ode_pdF_pdt = ode.1D(y=F,
                  times=c(0, 1),
                  func=func_pdF_pdt,
                  parms=F_parameters,
                  nspec=1,
                  names=c("F"))
  F = ode_pdF_pdt[2, -1]
  return(F)
}

####################################
### Define the intial contidions ###
####################################
func_init = function(nx, K, coef_K=2, q_intro=0.50, intro_pop=FALSE, intro_drive=TRUE, dType="TADS"){
  ### input
  # nx = number of steps or discrete intervals along the spatial x-axis
  # K = maximum population size in each interval, i.e. carrying capacity
  # a = critical population size below which population growth is negative
  # q_intro = initial frequency of the drive allele in the points of introduction (max at 0.5, i.e. 100% drive heterozygotes)
  # intro_pop = introduce a population below carrying capacity [TRUE/FALSE]
  # intro_drive = introduce the gene drive at q_intro [TRUE/FALSE]
  # dType = drive type; select from "W_shredder" or "TADS"
  if (intro_pop){
    n = rep(0, nx)
    n[round(nx/2)] = coef_K*K ### dounbling the density so that the other individuals can disperse and the Allee effect won't initially crash the population
  } else {
    n = rep(K, nx)
  }
  if (intro_drive){
    q = rep(0, nx)
    q[round(nx/2)] = q_intro
  } else {
    q = rep(0, nx)
  }
  if (dType=="W_shredder"){
    out = list(n=n, q=q)
  } else if (dType=="TADS"){
    ### initialise genotype frequencies across the landscape
    a = rep(1, nx) - q
    b = k = d = e = f = g = h = rep(0, nx)
    ### introduce TTDd drive carriers (remember that q_max = 0.5, i.e. 100% drive heterozygotes)
    b = 2*q
    ### output
    out = list(n=n, q=q, a=a, b=b, k=k, d=d, e=e, f=f, g=g, h=h)
  } else {
    out = list(n=n, q=q)
  }
  return(out)
}

##########################
### Execution function ###
##########################
func_run = function(nt, nx, dx, K, D, r, dType, c, coef_K=2, q_intro=0.5, intro_pop=FALSE, intro_drive=TRUE, threshold=0.05, progress.bar=FALSE){
  ### initial conditions
  init = func_init(nx=nx, K=K, coef_K=coef_K, q_intro=q_intro, intro_pop=intro_pop, intro_drive=intro_drive, dType=dType)
  q = init$q
  n = init$n
  u = n*q
  array_n = array(NA, dim=c(nt+1, nx))
  array_q = array(NA, dim=c(nt+1, nx))
  array_u = array(NA, dim=c(nt+1, nx))
  array_n[1, ] = n
  array_q[1, ] = q
  array_u[1, ] = u
  if (dType=="TADS"){
    a=init$a ### P(TTDD)
    b=init$b ### P(TTDd)
    k=init$k ### P(Ttdd)
    d=init$d ### P(TtDd)
    e=init$e ### P(Ttdd)
    f=init$f ### P(ttDD)
    g=init$g ### P(ttDd)
    h=init$h ### P(ttdd)
    list_geno_freqs = list(a=a, b=b, k=k, d=d, e=e, f=f, g=g, h=h)
    list_conversions = list(c_f=(0.95*c), c_m=c)
  }
  ### iterate across generations (NOTE: cannot wrap these into a function because deSolve::ode.1D does not seem to function as expected with parms)
  if (progress.bar){pb = txtProgressBar(min=2, max=(nt+1), char="=", style=3)}
  for (t in 2:(nt+1)){
    ### prepare PDEs input
    if ((dType=="W_shredder") | (dType=="X_shredder")){
      q_parameters = list(D=D, dx=dx, dType=dType, c=c, n=n)
      # u_parameters = list(q=q, D=D, dx=dx, dType=dType, c=c, n=n, r=r, K=K)
      n_parameters = list(q=q, D=D, r=r, K=K, dx=dx, dType=dType, c=c)
    } else if (dType=="TADS"){
      q_parameters = list(D=D, dx=dx, dType=dType, b=b, k=k, e=e, f=f, g=g, h=h, v=c, n=n)
      # u_parameters = list(q=q, D=D, dx=dx, dType=dType, b=b, k=k, e=e, f=f, g=g, h=h, v=c, n=n, r=r, K=K)
      n_parameters = list(q=q, D=D, r=r, K=K, dx=dx, dType=dType, k=k, f=f, h=h)
    }
    ### numerically solve PDEs
    ode_pdq_pdt = ode.1D(y=q,
                    times=c(0, 1),
                    func=func_pdq_pdt,
                    parms=q_parameters,
                    nspec=1,
                    names=c("q"),
                    method="lsodes")
    # ode_pdu_pdt = ode.1D(y=u,
    #                 times=c(0, 1),
    #                 func=func_pdu_pdt,
    #                 parms=u_parameters,
    #                 nspec=1,
    #                 names=c("u"),
    #                 method="lsodes")
    ode_pdn_pdt = ode.1D(y=n,
                    times=c(0, 1),
                    func=func_pdn_pdt,
                    parms=n_parameters,
                    nspec=1,
                    names=c("n"),
                    method="lsodes")
    ### derive additional frequencies in TADS
    if (dType=="TADS"){
      ### Solve PDEs of each genotype with fitness equivalent to the change in frequencies to the previous to the next generation
      ### using the mechanistic TADS recursive function - func_next_gen_freqs()
      geno_freqs = func_next_gen_freqs(list_geno_freqs, list_conversions)
      list_geno_freqs = list(a=geno_freqs$a, b=geno_freqs$b, k=geno_freqs$k, d=geno_freqs$d, e=geno_freqs$e, f=geno_freqs$f, g=geno_freqs$g, h=geno_freqs$h)
      list_geno_w = list(a=(geno_freqs$a-a), b=(geno_freqs$b-b), k=(geno_freqs$k-k), d=(geno_freqs$d-d), e=(geno_freqs$e-e), f=(geno_freqs$f-f), g=(geno_freqs$g-g), h=(geno_freqs$h-h))
      for (i in 1:length(list_geno_freqs)){
        F = unlist(list_geno_freqs[i])
        w = unlist(list_geno_w[i])
        F_new = tryCatch(func_disperse_TADS_genotypes(F=F, w=w, dx=dx, D=D), error=function(e){rep(0,rep=length(F))})
        list_geno_freqs[i] = list(F_new)
      }
      ### update genotype frequencies
      a = round(list_geno_freqs$a, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      b = round(list_geno_freqs$b, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      k = round(list_geno_freqs$k, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      d = round(list_geno_freqs$d, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      e = round(list_geno_freqs$e, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      f = round(list_geno_freqs$f, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      g = round(list_geno_freqs$g, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
      h = round(list_geno_freqs$h, 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
    }
    ### remove the first elements which are incorrectly computed boundary position
    q = round(ode_pdq_pdt[2, -1], 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
    # u = round(ode_pdu_pdt[2, -1], 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
    n = round(ode_pdn_pdt[2, -1], 5) ### rounding to 5 decimal digits to avoid weird random fluctuations even when we start at a flat q=0
    u = n * q
    # n = n + u
    # q = u/n
    ### detect if the populaion density has dropped to essentially zero, i.e. if n < threshold
    n[n<threshold] = 0 ### the population density is essentially zero
    q[(n<threshold) | is.na(q)] = 0 ### when the population density is below the threshold then the drive is lost
    u[(n<threshold) | is.na(u)] = 0 ### when the population density is below the threshold then the drive is lost
    
    ### drop q to zero if q < 1/(2*n)
    # q[q < 1/(2*ceiling(n))] = 0
    
    if (dType=="TADS"){
      ### update genotype frequencies
      a[n<threshold] = 0
      b[n<threshold] = 0
      k[n<threshold] = 0
      d[n<threshold] = 0
      e[n<threshold] = 0
      f[n<threshold] = 0
      g[n<threshold] = 0
      h[n<threshold] = 0
    }
    if (progress.bar){setTxtProgressBar(pb, t)}
    ### allocate
    array_q[t, ] = q
    array_u[t, ] = u
    array_n[t, ] = n
  }
  if (progress.bar){close(pb)}
  # return(list(array_n=array_n, array_q=array_q))
  return(list(array_n=array_n, array_q=array_q, array_u=array_u))
}

#########################
### Plotting function ###
#########################
plot_spectrum_densities_and_frequencies = function(array_n, array_q, nt, nx, K, spec, period){
  par(mfrow=c(2,1))
  ### spectral graph
  plot(spec$freq, spec$spec, type="l", xlab="Frequency\n(multiply by the no. of generations to get generations)", ylab="Spectrum")
  lines(x=c(period,period)/nt, y=c(0,max(spec$spec)), lty=2, lwd=2, col="red")
  text(x=period/nt, y=max(spec$spec), pos=4, lab=paste0("period=", round(period,2)))
  par(new=TRUE)
  plot(vec_n, type="l", col="blue", xaxt="n", yaxt="n", xlab="", ylab="", main="")
  grid()
  legend("topright", legend=c("Spectral decomposition", "Population density across generations", paste0("Period = ", round(period,2))), lty=c(1,1,2), lwd=2, col=c("black", "blue", "red"))

  ### line plot across time (middle of the landscape)
  par(mar=c(5,5,2,5))
  plot(array_n[,round(nx/2)], ylim=c(0, K), type="l", las=2, xlab="Time", ylab="n")
  par(new=TRUE)
  plot(array_q[,round(nx/2)], ylim=c(0, 1), type="l", col="red", xaxt="n", yaxt="n", xlab="", ylab="")
  axis(side=4, at=seq(0,1,by=0.1), lab=seq(0,1,by=0.1), las=2)
  mtext(side=4, text="q", padj=4)
  grid()
  legend("topright", legend=c("Population density (n)", "Drive allele frequency (q)"), lty=1, col=c("black", "red"))
}

########################################
### Computed and Fisher's velocities ###
########################################
func_velocity = function(array, D, r, thresh=0.5, by=2, lab="", plot.out=FALSE){
  ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ### test
  # array = array_q
  # D = 10
  # r = 0.25
  # thresh = 0.5
  # by = 2
  # lab = "Drive allele frequency"
  # plot.out = TRUE
  ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ### detect the edge of the growing wave front
  array_idx = array > thresh
  ### loop across time
  vec_v = c()
  vec_t = seq(from=by, to=nt, by=by)
  for (t in vec_t){
    t1 = t-by
    t2 = t
    x1 = c(1:nx)[array_idx[t1, ] == 1][1]
    x2 = c(1:nx)[array_idx[t2, ] == 1][1]
    vec_v = c(vec_v, (x1 - x2)/(t2 - t1))
  }

  ### remove NAs
  idx = !is.na(vec_v)
  vec_v = vec_v[idx]
  vec_t = vec_t[idx]
  ### remove negative velocities since we're expeting only increase across space and time
  idx = vec_v > 0
  vec_v = vec_v[idx]
  vec_t = vec_t[idx]
  ### finally select the median velocity
  velocity = median(vec_v)
  ### predicted velocy
  Fishers_velocity = round(2*sqrt(D*r), 2)
  ### plot
  if(plot.out){
    par(ask=TRUE)
    image(array_idx)
    plot(x=vec_t, y=vec_v, xlab="Time", ylab=paste0("Velocity of\n", lab), type="l"); grid()
    legend("top", bg="gray", legend=c(paste0("Fisher's velocity = ", Fishers_velocity),
                          paste0("Median computed velocity = ", round(velocity, 2))))
  }
  return(list(array_idx=array_idx, Fishers_velocity=Fishers_velocity, velocity=velocity))
}

############################################
### Extract drive allele wave statistics ###
############################################
func_q_wave_stats = function(array_q, K, nx, epsilon=1e-2){
  ###@@@@@@@@@@@@@@@@@@@@@@@
  ### Test
  # K = K
  # nx = nx
  # epsilon = 1e-2
  # plot(array_n[t,]/K, type="b", ylim=c(0,1))
  # points(array_q[t,], col="red")
  # lines(array_q[t,], col="red")
  # grid()
  ###@@@@@@@@@@@@@@@@@@@@@@@
  ### Find the time at which we maximise q
  t = c(1:(nt+1))[apply(array_q, MAR=1, FUN=max) == max(array_q)][1]
  ### Find the height of the invading drive allele wave
  height = max(array_q[t,])
  coor_height = c(1:nx)[array_q[t,] == height][1]
  ### Find the leading edge of the drive allele wave (x0)
  for (x in coor_height:1){
    if (array_q[t,x] == 0){
      x0 = x
      break
    }
  }
  ### Find the trailing edge of the drive allele wave (x1)
  for (x in coor_height:(nx/2)){
    if (array_q[t,x] == 0){
      x1 = x
      break
    }
  }
  if ((!exists("x0")) | (!exists("x1"))) {
    x0 = NA
    x1 = NA
    width = NA
    coor_ixn = NA
    q_ixn = NA
    steepness = NA
  } else {
    ### Find the width of the invading drive allele wave
    width = x1 - x0
    ### Find the intersection between n and q
    vec_n = array_n[t, x0:x1]
    vec_q = array_q[t, x0:x1]
    idx = (vec_n != 0) & (vec_q != 0)
    vec_n = vec_n[idx]
    vec_q = vec_q[idx]
    vec_n_diff_q = abs((vec_n/K) - vec_q)
    coor_ixn = c(x0:x1)[idx][vec_n_diff_q==min(vec_n_diff_q)]
    q_ixn = array_q[t, coor_ixn]
    ### Compute the steepness of the wave
    steepness = height / width
  }
  ### Output
  return(list(height=height, width=width, steepness=steepness, coor_ixn=coor_ixn, q_ixn=q_ixn))
}

############
### TEST ###
############
### test input parameters
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
progress.bar = TRUE  
D = 10
r = 1.6
c = 1
dType = "W_shredder"

### test run
test = func_run(nt=nt, nx=nx, dx=dx, K=K, D=D, r=r, dType=dType, c=c, coef_K=coef_K, q_intro=q_intro, intro_pop=intro_pop, intro_drive=intro_drive, threshold=threshold, progress.bar=progress.bar)
array_q = test$array_q
array_n = test$array_n
array_u = test$array_u
q = array_q[(nt+1), ]
n = array_n[(nt+1), ]
u = array_u[(nt+1), ]

### test plots
par(ask=TRUE)

t_plot = round(seq(1, nt, length=10))
par(mar=c(5,5,2,5))
plot(0, ylim=c(0, 1), xlim=c(0,nx), type="n", yaxt="n", xlab="Time", ylab="Population desnity")
axis(side=2, at=seq(0,1,by=0.2), lab=round(seq(0,1,by=0.2)*K))
axis(side=4, at=seq(0,1,by=0.2), lab=seq(0,1,by=0.2))
mtext(side=4, text="Drive allele frequency", padj=4)
grid()
for (t in t_plot){
  lines(array_n[t, ]/K)
  lines(array_q[t, ], col="red")
}

### gradient plots
filled.contour(x=1:nrow(array_n), z=array_n, y=1:ncol(array_n), zlim=c(0,K), xlab="Time", ylab="Spatial axis", main="Population density (n)")
filled.contour(x=1:nrow(array_q), z=array_q, y=1:ncol(array_q), zlim=c(0,1), xlab="Time", ylab="Spatial axis", main="Drive allele frequency (q)")

### test find the velocities
if(intro_pop==TRUE){
  func_velocity(array=array_n, D=D, r=r, plot.out=TRUE, lab="population density")
}
if (intro_drive==TRUE) {
  print(paste0("q_max = ", max(array_q)))
  print(paste0("v_wt = ", round(2*sqrt(D*r),2)))
  if(dType=="W_shredder"){
    func_velocity(array=array_q, D=D, r=c/(8-4*c), plot.out=TRUE, lab="drive allele frequency")
    print(paste0("v_dr = ", round(2*sqrt(D*c/(8-4*c)),2)))
  } else if(dType=="X_shredder"){
    func_velocity(array=array_q, D=D, r=c/(2-c+c*0), plot.out=TRUE, lab="drive allele frequency")
    print(paste0("v_dr = ", round(2*sqrt(D*c/(2-c+c*0)),2)))
  }
}
### Hence the velocities are correct!!!
### Therefore chasing is more to do with the drive wave penetration and re-establishment capacity
### of the wild-types into the free space!!!
### Now I'm thinking we should look more into the steepness of the wave rather than its velocity!

### Nope, the drive allele frequencies do not seem to be a big factor
### And so let's try to plot n and q per generation
plot.out = TRUE
if (plot.out) {
  par(ask=FALSE)
  for (t in 1:(nt+1)){
    png(paste0("drive_PDE-id_", t, ".png"), width=700, height=500)
    t_plot = round(seq(1, nt, length=10))
    par(mar=c(5,5,2,5))
    plot(0, ylim=c(0, 1), xlim=c(0,nx), type="n", yaxt="n", xlab="Spatial axis", ylab="Population desnity", main=paste0("Generation ", t-1))
    axis(side=2, at=seq(0,1,by=0.2), lab=round(seq(0,1,by=0.2)*K))
    axis(side=4, at=seq(0,1,by=0.2), lab=seq(0,1,by=0.2))
    mtext(side=4, text="Drive allele frequency", padj=4)
    grid()
    lines(array_n[t, ]/K)
    lines(array_q[t, ], col="red")
    v_wt = func_velocity(array=array_n, thresh=threshold, D=D, r=r, plot.out=FALSE)
    if (dType=="W_shredder"){
      stats = func_q_wave_stats(array_q=array_q[t, , drop=FALSE], K=K, nx=nx)
      v_q = func_velocity(array=array_q, thresh=threshold, D=D, r=c/(8-4*c), plot.out=FALSE)
      legend("left", legend=c("Population density", "Drive allele frequency", paste0("q_max=",round(stats$height,2)), paste0("width=",round(stats$width,2)), paste0("WT velocity=", round(v_wt$velocity,2)), paste0("Drive velocity=", round(v_q$velocity,2))), lty=c(1,1,NA,NA,NA,NA), lwd=2, col=c("black","red",NA,NA,NA,NA), bty="n")
    } else if (dType=="X_shredder"){
      stats = func_q_wave_stats(array_q=array_q[t, , drop=FALSE], K=K, nx=nx)
      v_q = func_velocity(array=array_q, thresh=threshold, D=D, r=c/(2-c+c*0), plot.out=FALSE)
      legend("left", legend=c("Population density", "Drive allele frequency", paste0("q_max=",round(stats$height,2)), paste0("width=",round(stats$width,2)), paste0("WT velocity=", round(v_wt$velocity,2)), paste0("Drive velocity=", round(v_q$velocity,2))), lty=c(1,1,NA,NA,NA,NA), lwd=2, col=c("black","red",NA,NA,NA,NA), bty="n")
    }
    dev.off()
  }
  library(av)
  av_encode_video(system(paste0("ls -1v drive_PDE-id_*.png"), intern=TRUE), framerate=10, output=paste0("drive_PDE-D", D, "-r", r, "-c", c, "-", dType, ".mp4"))
  system("rm drive_PDE-id_*.png")
  ### plot at a single timepoint for illustration
  if (dType=="W_shredder"){
    png(paste0("drive_PDE-D", D, "-r", r, "-c", c, "-", dType, ".png"), width=900, height=600)
    t = 96
    par(mar=c(5,5,2,5))
    plot(0, ylim=c(0, 1), xlim=c(0,nx), type="n", yaxt="n", xlab="Time", ylab="Population desnity", main=paste0("Generation ", t-1))
    axis(side=2, at=seq(0,1,by=0.2), lab=round(seq(0,1,by=0.2)*K))
    axis(side=4, at=seq(0,1,by=0.2), lab=seq(0,1,by=0.2))
    mtext(side=4, text="Drive allele frequency", padj=4)
    grid()
    lines(array_n[t, ]/K)
    lines(array_q[t, ], col="red")
    stats = func_q_wave_stats(array_q=array_q[t, , drop=FALSE], K=K, nx=nx)
    legend("left", legend=c("Population density", "Drive allele frequency", paste0("q_max=",round(stats$height,2)), paste0("width=",round(stats$width,2))), lty=c(1,1,NA,NA), lwd=2, col=c("black", "red", NA, NA), bty="n")
    dev.off()
  }
  if (dType=="X_shredder"){
    png(paste0("drive_PDE-D", D, "-r", r, "-c", c, "-", dType, ".png"), width=900, height=600)
    t = 36
    par(mar=c(5,5,2,5))
    plot(0, ylim=c(0, 1), xlim=c(0,nx), type="n", yaxt="n", xlab="Time", ylab="Population desnity", main=paste0("Generation ", t-1))
    axis(side=2, at=seq(0,1,by=0.2), lab=round(seq(0,1,by=0.2)*K))
    axis(side=4, at=seq(0,1,by=0.2), lab=seq(0,1,by=0.2))
    mtext(side=4, text="Drive allele frequency", padj=4)
    grid()
    lines(array_n[t, ]/K)
    lines(array_q[t, ], col="red")
    stats = func_q_wave_stats(array_q=array_q[t, , drop=FALSE], K=K, nx=nx)
    legend("left", legend=c("Population density", "Drive allele frequency", paste0("q_max=",round(stats$height,2)), paste0("width=",round(stats$width,2))), lty=c(1,1,NA,NA), lwd=2, col=c("black", "red", NA, NA), bty="n")
    dev.off()
  }
  ### plot velocities as a function of q
  png(paste0("drive_PDE-velocities-D", D, "-r", r, "-c", c, "-", dType, ".png"), width=900, height=600)
  vec_q = seq(0,1,length=1000)
  growth_rate = 2*sqrt(D*r)
  suppression_rate = 2*sqrt(D*2*r*(1-((2-c-c*vec_q)/(2-c))))
  growth_suppression_rates = 2*sqrt(as.complex( D * r * (1 - 2*vec_q)))
  xlim = c(0,1)
  ylim = c( -max(c(suppression_rate, Im(growth_suppression_rates))), max(Re(growth_suppression_rates)))
  plot(x=xlim, y=ylim, type="n", xlab="Drive allele frequency (q)", ylab="Population velocity")
  grid()
  lines(x=vec_q, y=Re(growth_suppression_rates), col="black") ### both growth and suppression (Real part)
  lines(x=vec_q, y=-Im(growth_suppression_rates), col="red") ### both growth and suppression (Imaginary part as negative velocities)
  abline(h=growth_rate, lty=2, col="green") ### Growth alone (positive velocity)
  lines(x=vec_q, y=-suppression_rate, lty=2, col="orange") ### Suppression alone (negative velocity)
  abline(v=0.5, lty=4, col="grey") ### at q=0.5
  abline(h=-2*sqrt(D*c/(8-4*c)), lty=2, col="blue") ### W-shredder (negative velocity)
  lines(x=vec_q, y=-2*sqrt(D*c/(2-c+c*vec_q)), lty=2, col="purple") ### X-shredder (negative velocity)
  legend("topright", legend=c("+growth", "-growth", "growth alone", "suppression alone", "W-shredder", "X-shredder", "q=0.5"),
                     lty=c(1,1,2,2,2,2,4), lwd=2,
                     col=c("black", "red", "green", "orange", "blue", "purple", "grey"))
  dev.off()
}

### It looks like the growth just needs to overcome supression to re-establish the population

### Finding the r threshold at which W-shredder and X-shredder transition from crash to chase
vec_n = seq(0, K, length=1000)
vec_q = seq(0, 1, length=1000)
vec_dGrowth_dt  = r * vec_n * (1 - (vec_n/K))
mat_dSupression_dt = matrix(NA, nrow=length(vec_q), ncol=length(vec_n))
for (i in 1:length(vec_q)){
  q = vec_q[i]
  vec_dSupression_dt  = 2 * r * vec_n * (1 - ((2-c-c*q)/(2-c)))
  mat_dSupression_dt[i, ] = vec_dSupression_dt
}

xlim = c(0, K)
ylim = c(0, max(mat_dSupression_dt))
plot(x=xlim, y=ylim, type="n", xlab="n", ylab="d./dt")
lines(x=vec_n, y=vec_dGrowth_dt)
for (i in seq(1, length(vec_q), length=10)){
  q = mat_dSupression_dt[i,]
  lines(x=vec_n, y=q, col="red")
  text(x=4, y=max(q), pos=4, lab=paste0("q=", round(vec_q[i],2)))
}
