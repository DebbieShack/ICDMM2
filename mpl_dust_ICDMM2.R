#------------------------------------------------
#' Model Parameter List Creation
#'
#' \code{model_param_list_create} creates list of model parameters to be used
#' within \code{equilibrium_init_create}
#'
#' @param n_days Number of days for which the simulation will run. Default = 5*365
#' @param age_vector Lower bound of each age category. See function for default. 
#' @param het_brackets Number of biting heterogeneity groups. Default = 5
#' @param lag_rates #Number of sub-compartments within FOI and FOIv which approximate delay-differential equation. Higher values are a closer approximation, but computationally more expensive. Default = 10. 
#' @param num_int Number of intervention parameters.  Default = 1 (No intervention)
#' @param tsd Number of time-steps per day. Fewer is faster, more better approximates the continuous solution. Default = 4.
#' @param admin_unit Specifying an admin2 region will automatically add average-rainfall-based seasonality into the model
#' @param country If admin2 region name is not unique, the country also needs to be stated. 
#' @param eta Death rate for exponential population distribtuion, i.e. 1/Mean Population Age. Default = 0.0001305
#' @param rho Age-dependent biting parameter. Default = 0.85
#' @param a0 Age-dependent biting parameter. Default = 2920
#' @param sigma2 Variance of the log heterogeneity in biting rates. Default = 1.67
#' @param max_age Maximum age in days. Default = 100*365
#' @param rA Rate of leaving asymptomatic infection. Default = 0.00512821
#' @param rT Rate of leaving treatment. Default = 0.2
#' @param rD Rate of leaving clinical disease. Default = 0.2
#' @param rU Rate of recovering from subpatent infection. Default = 0.00906627
#' @param rP Rate of leaving prophylaxis. Default = 0.06666667
#' @param dE Latent period of human infection. Default = 12
#' @param delayGam Lag from parasites to infectious gametocytes. Default = 12.5
#' @param cD Untreated disease contribution to infectiousness. Default = 0.0676909
#' @param cT Treated disease contribution to infectiousness. Default =   0.322 * cD
#' @param cU Subpatent disease contribution to infectiousness. Default = 0.006203
#' @param gamma1 Parameter for infectiousness of state A. Default = 1.82425
#' @param d1 Minimum probability due to maximum immunity. Default = 0.160527
#' @param dID Inverse of decay rate. Default = 3650
#' @param ID0 Scale parameter. Default = 1.577533
#' @param kD Shape parameter. Default = 0.476614
#' @param uD Duration in which immunity is not boosted. Default = 9.44512
#' @param aD Scale parameter relating age to immunity. Default = 8001.99
#' @param fD0 Time-scale at which immunity changes with age. Default = 0.007055
#' @param gammaD Shape parameter relating age to immunity. Default = 4.8183
#' @param alphaA PCR detection probability parameters state A. Default = 0.757
#' @param alphaU PCR detection probability parameters state U. Default = 0.186
#' @param b0 Maximum probability due to no immunity. Default = 0.590076
#' @param b1 Maximum relative reduction due to immunity. Default = 0.5
#' @param dB Inverse of decay rate. Default = 3650
#' @param IB0 Scale parameter. Default = 43.8787
#' @param kB Shape parameter. Default = 2.15506
#' @param uB Duration in which immunity is not boosted. Default = 7.19919
#' @param phi0 Maximum probability due to no immunity. Default = 0.791666
#' @param phi1 Maximum relative reduction due to immunity. Default = 0.000737
#' @param dCA Inverse of decay rate. Default = 10950
#' @param IC0 Scale parameter. Default = 18.02366
#' @param kC Shape parameter. Default = 2.36949
#' @param uCA Duration in which immunity is not boosted. Default = 6.06349
#' @param PM New-born immunity relative to mother’s. Default = 0.774368
#' @param dCM Inverse of decay rate of maternal immunity. Default = 67.6952
#' @param delayMos Extrinsic incubation period. Default = 10
#' @param tau1 Duration of host seeking, assumed to be constant between species. Default = 0.69
#' @param tau2 Duration of mosquito resting after feed. Default = 2.31
#' @param mu0 Daily mortality of adult mosquitos. Default = 0.132
#' @param Q0 Anthrophagy probability. Default = 0.92
#' @param chi Endophily probability. Default = 0.86
#' @param bites_Bed Percentage of bites indoors and in bed. Default = 0.89
#' @param bites_Indoors Percentage of bites indoors . Default = 0.97
#' @param muEL Per capita daily mortality rate of early stage larvae (low density). Default = 0.0338
#' @param muLL Per capita daily mortality rate of late stage larvae (low density). Default = 0.0348
#' @param muPL Per capita daily mortality rate of pupae. Default = 0.249
#' @param dEL Development time of early stage larvae. Default = 6.64
#' @param dLL Development time of late stage larvae. Default = 3.72
#' @param dPL Development time of pupae. Default = 0.643
#' @param gammaL Relative effect of density dependence on late instars relative to early instars. Default = 13.25
#' @param km Seasonal carrying capacity. Default = 11
#' @param cm Seasonal birth rate. Default = 0.05
#' @param betaL Number of eggs laid per day per mosquito. Default = 21.2
#' @param itn_cov The proportion of people that use an ITN. Default = 0
#' @param DY Duration of year (days). Default = 365
#' @param daily_itn_cov Vector of daily ITN coverage (population proportion). If single value is given, it is assumed constant over simulation period. 
#' @param d_itn0 Probability of dying with an encounter with ITN (max). Default = 0.41
#' @param r_itn0 Probability of repeating behaviour with ITN (max). Default = 0.56
#' @param r_itn1 Probability of repeating behaviour with ITN (min). Default = 0.24
#' @param itn_half_life ITN half life. Default =   2.64 * DY
#' @param itn_interval How long before ITN is repeated
#' @param daily_smc_cov Vector of daily SMC coverage (population proportion). If single value given, model assumes constant over simulation period. 
#' @param smc_days Days (number of days from simulation onset) on which SMC treatments are administered. If not stated, model assumes 'optimal' dates 
#' @param n_smc_doses If SMC_days not stated - how many doses per SMC season? Default = 4.
#' @param days_between_smc_doses If smc_days not stated - how many days between doses? Default = 30
#' @param smc_upper_age Maximum age to receive SMC doses
#' @param shape_smc Shape parameter of Weibull distribution describing smc prophylactic decay. Default = 4.3
#' @param scale_smc Scale parameter of Weibull distribution describing smc prophylactic decay. Default = 38.1
#' @param alpha_smc Proportion of existing infections cleared by smc. Default = 0.95
#' @param smc_clearance_lag Number of days between smc administration and clearance of any existing infection (days). Default = 5. 
#' @param smc_rel_c Relative infectiousness of existing infections (after receipt of smc dose) during the period between dose and clearance.  
#' @param ... Any other parameters needed for non-standard model. If they share the same name
#' as any of the defined parameters \code{model_param_list_create} will stop. You can either write
#' any extra parameters you like individually, e.g. model_param_list_create(extra1 = 1, extra2 = 2)
#' and these parameteres will appear appended to the returned list, or you can pass explicitly
#' the ellipsis argument as a list created before, e.g. model_param_list_create(...=list(extra1 = 1, extra2 = 2))
#'
#' @export

if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

get_parameters <- function(
  n_days = 5*365,
  tsd = 4, #Time steps per day
  admin_unit = NULL,
  country = NULL,
  age_vector = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
  het_brackets = 5,
  lag_rates = 10, #Number of sub-compartments within FOI and FOIv which approximate delay-differential equation. Higher values are a closer approximation, but computationally more expensive
    # age, heterogeneity in exposure,
  eta = 1/(21*365),
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100*365,
  #  rate of leaving infection states
  rA = 1/195,
  rT = 0.2,
  rD = 0.2,
  rU = 1/110.299,
  rP = 1/15,
  #  human latent period and time lag from asexual parasites to
  dE  = 12,
  delayGam = 12.5,
  # human infectiousness to mosquitoes
  cD  = 0.0676909,
  cT  =  0.322 * cD,
  cU  = 0.006203,
  gamma1  = 1.82425,
  #  Immunity reducing probability of detection
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  # Immunity reducing probability of infection
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  # Immunity reducing probability of clinical disease
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  # entomological parameters
  delayMos = 10,
  tau1 = 0.69,
  tau2 = 2.31,
  mu0 = 0.132,
  Q0 = 0.92,
  chi = 0.86,
  bites_Bed = 0.89,
  bites_Indoors = 0.97,
  # larval parameters daily density dependent mortality rate of egg
  muEL = 0.0338,
  muLL = 0.0348,
  muPL = 0.249,
  dEL = 6.64,
  dLL = 3.72,
  dPL = 0.643,
  gammaL = 13.25,
  km = 11,
  cm = 0.05,
  betaL = 21.2,
  # intervention parameters
  DY = 365,
  ...
  
){
  # set up param list
  params <- list()
  
  # catach extra params and place in list
  extra_param_list <- list(...)
  if(length(extra_param_list)>0){
    if(is.list(extra_param_list[[1]])){
      extra_param_list <- extra_param_list[[1]]
    }
  }
  
  ## DEFAULT PARAMS
  params$n_days <- n_days
  params$tsd <- tsd
  params$n_ts <- n_days*tsd
  # duration of year
  params$DY <- DY
    # age, heterogeneity in exposure
  params$age_vector <- age_vector
  params$het_brackets <- het_brackets
  params$eta <- eta
  params$rho <- rho
  params$a0 <- a0
  params$sigma2 <- sigma2
  params$max_age <- max_age
  na <- as.integer(length(age_vector))  # number of age groups
  nh <- as.integer(het_brackets) 
  params$na <- na
  params$nh <- nh
  
  
  # rate of leaving infection states
  params$rA <- rA
  params$rT <- rT
  params$rD <- rD
  params$rU <- rU
  params$rP <- rP
  
  # human latent period and time lag from asexual parasites to
  # infectiousness
  params$dE <- dE
  params$delayGam <- delayGam
  
  # infectiousness to mosquitoes
  params$cD <- cD
  params$cT <- cT
  params$cU <- cU
  params$gamma1 <- gamma1
  
  # Immunity reducing probability of detection
  params$d1 <- d1
  params$dID <- dID
  params$ID0 <- ID0
  params$kD <- kD
  params$uD <- uD
  params$aD <- aD
  params$fD0 <- fD0
  params$gammaD <- gammaD
  
  # PCR prevalence parameters
  params$alphaA <- alphaA
  params$alphaU <- alphaU
  
  # anti-infection immunity
  params$b0 <- b0
  params$b1 <- b1
  params$dB <- dB
  params$IB0 <- IB0
  params$kB <- kB
  params$uB <- uB
  
  # clinical immunity
  params$phi0 <- phi0
  params$phi1 <- phi1
  params$dCA <- dCA
  params$IC0 <- IC0
  params$kC <- kC
  params$uCA <- uCA
  params$PM <- PM
  params$dCM <- dCM
  
  # entomological parameters
  params$delayMos <- delayMos
  params$tau1 <- tau1
  params$tau2 <- tau2
  params$mu0 <- mu0
  params$Q0 <- Q0
  params$bites_Bed <- bites_Bed
  params$bites_Indoors <- bites_Indoors
  params$fv0 <- 1 / (tau1 + tau2)
  params$av0 <- Q0 * params$fv0 # daily feeeding rate on humans
  params$Surv0 <- exp(-mu0 * delayMos) # probability of surviving incubation period
  params$p10 <- exp(-mu0 * tau1)  # probability of surviving one feeding cycle
  params$p2 <- exp(-mu0 * tau2)  # probability of surviving one resting cycle
  
  # larval parameters
  params$muEL <- muEL
  params$muLL <- muLL
  params$muPL <- muPL 
  params$dEL <- dEL
  params$dLL <- dLL
  params$dPL <- dPL
  params$gammaL <- gammaL
  params$km <- km
  params$cm <- cm
  params$betaL <- betaL
  # {White et al. 2011 Parasites and Vectors}
  params$eov <- betaL/mu0 * (exp(mu0/params$fv0) - 1)
  params$b_lambda <- (gammaL * muLL/muEL - dEL/dLL + (gammaL - 1) * muLL * dEL)
  params$lambda <- -0.5 * params$b_lambda +
    sqrt(0.25 * params$b_lambda^2 + gammaL * betaL * muLL * dEL/(2 * muEL * mu0 * dLL * (1 + dPL * muPL)))
  

  #Additional parameters for dust model
  params$lag_rates <- lag_rates
  params$dt <- 1/tsd
  
  #Track which parameters have been set so far.
  params$smc_set <- 0
  params$itn_set <- 0
  params$equilibrium_set <- 0
  params$seasonality_set <- 0
  
  ##Default parameters to be potentially overriden
  params$EIP_tempsens <- 0 #odin.dust likes numerical variable types
  params$mu_tempsens <- 0

  #-------------------------------------------------------------

  #check that none of the spare parameters in the extra
  if(sum(!is.na(match(names(extra_param_list),names(params))))!=0){

    stop (message(cat("Extra params in ... share names with default param names. Please check:\n",
                      names(extra_param_list)[!is.na(match(names(extra_param_list),names(params)))]
    )
    ))
  }

  return(append(params,extra_param_list))
}
cD <- 0.0676909
DY <- 365

set_equilibrium <- function(params, country = NULL, admin_unit = NULL, ft = 0,
                                         init_EIR, seasonality_loc = NULL)
{
  #Set intervention group coverages
  if(params$smc_set == 0 & params$itn_set == 0){
    num_int <- 1
    cov <- 1
  } else if(params$itn_set == 1 & params$smc_set == 0){
    num_int <- 2
    cov <- c((1 - params$max_itn_cov), params$max_itn_cov)
  } else if(params$itn_set == 0 & params$smc_set == 1){
    num_int <- 3
    cov <- c((1 - params$max_smc_cov),0,params$max_smc_cov)
  } else if(params$itn_set == 1 & params$itn_set == 1){
    num_int <- 4
    cov <- c((1 - params$max_itn_cov) * (1 - params$max_smc_cov), 
                    params$max_itn_cov * (1 - params$max_smc_cov), 
                    (1 - params$max_itn_cov) * params$max_smc_cov, 
                    params$max_itn_cov * params$max_smc_cov)
  }
  params$cov <- cov
  params$num_int <- num_int
  age_vector <- params$age_vector
  het_brackets <- params$het_brackets
  
  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(het_brackets)) stop("het_brackets provided is not numeric")
  if(!(is.null(country) | is.character(country))) stop("country specified is not character string")
  if(!(is.null(admin_unit) | is.character(admin_unit))) stop("admin_unit specified is not character string")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(init_EIR)) stop("init_EIR provided is not numeric")
  
  ## Handle parameters
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")
  
  ## population demographics
  age <- age_vector * params$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group
    
  }
  age_rate[na] = 0
  
  
  den <- 1/(1 + age_rate[1]/params$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + params$eta)  # proportion in each age_vector group
  }
  
  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years
  
  #Used to estimate PfPR_2_10 (as in Malariasimulation)
  upto_age02 <- max(which(age_vector < 2))  # index of age vector from which age is 2 years
  upto_age10 <- max(which(age_vector < 10))  # index of age vector before age is 10 years
  ## force of infection
  foi_age <- c()
  for (i in 1:na)
  {
    foi_age[i] <- 1 - (params$rho * exp(-age[i]/params$a0))  #force of infection for each age group
  }
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant
  
  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-params$sigma2/2 + sqrt(params$sigma2) * het_x)/sum(het_wt * exp(-params$sigma2/2 + sqrt(params$sigma2) * het_x))
  
  ## EIR
  EIRY_eq <- init_EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/params$DY
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq
  
  ## Immunity and FOI
  x_I <- den[1]/params$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - params$fD0)/(1 + (age/params$aD)^params$gammaD)
  
  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * params$DY) & age[i - 1] < (20 * params$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * params$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])
  
  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  ICM_eq <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * params$uB + 1) * x_I[i])/(1 + x_I[i]/params$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, params$b0,
                                            params$b0 * ((1 - params$b1)/(1 + (IB_eq[i, j]/params$IB0)^params$kB) + params$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * params$uD + 1) * x_I[i])/(1 + x_I[i]/params$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * params$uCA + 1) * x_I[i])/(1 + x_I[i]/params$dCA)
      p_det_eq[i, j] <- params$d1 + (1 - params$d1)/(1 + fd[i] * (ID_eq[i, j]/params$ID0)^params$kD)
      cA_eq[i, j] <- params$cU + (params$cD - params$cU) * p_det_eq[i, j]^params$gamma1
    }
  }
  # needs to be calculated after because it references ICA
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      ICM_init_eq[j] <- params$PM * (ICA_eq[age20l, j] + age_20_factor *
                                    (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
      ICM_eq[i, j] <- ifelse(i == 1,
                             ICM_init_eq[j], ICM_eq[i - 1,j])/(1 + x_I[i]/params$dCM)
    }
  }
  
  IC_eq <- ICM_eq + ICA_eq
  phi_eq <- params$phi0 * ((1 - params$phi1)/(1 + (IC_eq/params$IC0)^params$kC) + params$phi1)
  
  # human states
  gamma <- params$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(params$eta, age_rate[1:(na - 1)])
  
  betaT <- matrix(rep(params$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(params$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(params$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  
  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- params$rT * aT/betaP
  
  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]
  
  for (j in 1:nh)
  {
    for (i in 2:na)
    {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (params$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (params$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]
      
    }
  }
  
  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)
  
  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + params$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + params$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  
  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)
  
  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       params$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (params$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      FOIvij_eq[i, j] <- foi_age[i] * params$av0 * (params$cT * T_eq[i, j] + params$cD *
                                                   D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + params$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }
  
  # mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  Iv_eq <- FOIv_eq * params$Surv0/(FOIv_eq + params$mu0)
  Sv_eq <- params$mu0 * Iv_eq/(FOIv_eq * params$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq
  
  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * params$av0)
  
  # larval states
  K0 <- 2 * mv0 * params$dLL * params$mu0 * (1 + params$dPL * params$muPL) * params$gammaL * (params$lambda + 1)/(params$lambda/(params$muLL *
                                                                                                              params$dEL) - 1/(params$muLL * params$dLL) - 1)
  PL_eq <- 2 * params$dPL * params$mu0 * mv0
  LL_eq <- params$dLL * (params$muPL + 1/params$dPL) * PL_eq
  EL_eq <- (LL_eq/params$dLL + params$muLL* LL_eq * (1 + params$gammaL * LL_eq/K0))/(1/params$dEL - params$muLL * params$gammaL * LL_eq/K0)
  
  # add in final dimension - interventions
  cov <- params$cov
  
  mat <- matrix(0, na, nh)
  
  S_eq <- vapply(cov, FUN = function(x)
  {
    x * S_eq
  }, mat)
  T_eq <- vapply(cov, FUN = function(x)
  {
    x * T_eq
  }, mat)
  D_eq <- vapply(cov, FUN = function(x)
  {
    x * D_eq
  }, mat)
  A_eq <- vapply(cov, FUN = function(x)
  {
    x * A_eq
  }, mat)
  U_eq <- vapply(cov, FUN = function(x)
  {
    x * U_eq
  }, mat)
  P_eq <- vapply(cov, FUN = function(x)
  {
    x * P_eq
  }, mat)
  
  IB_eq = array(IB_eq, c(na, nh, num_int))
  ID_eq = array(ID_eq, c(na, nh, num_int))
  ICA_eq = array(ICA_eq, c(na, nh, num_int))
  ICM_eq = array(ICM_eq, c(na, nh, num_int))
  
  # TODO: Remove this part and put it as an edit to the equilibrium solution
  if(!is.null(params$ncc)){
    IB_eq = array(IB_eq, c(na, nh, num_int, params$ncc))
    ID_eq = array(ID_eq, c(na, nh, num_int, params$ncc))
    ICA_eq = array(ICA_eq, c(na, nh, num_int, params$ncc))
    ICM_eq = array(ICM_eq, c(na, nh, num_int, params$ncc))
    
    # add in final dimension - interventions
    all_rounds = params$MDA_grp_prop*params$MDA_cov
    ccov = c(all_rounds, 1-all_rounds)
    
    mat2 <- array(0, c(na,nh, num_int))
    S_eq <- vapply(ccov,FUN = function(x){x * S_eq},mat2)
    T_eq <- vapply(ccov,FUN = function(x){x * T_eq},mat2)
    D_eq <- vapply(ccov,FUN = function(x){x * D_eq},mat2)
    A_eq <- vapply(ccov,FUN = function(x){x * A_eq},mat2)
    U_eq <- vapply(ccov,FUN = function(x){x * U_eq},mat2)
    P_eq <- vapply(ccov,FUN = function(x){x * P_eq},mat2)
  }
  
  
  admin_units_seasonal <- load_file("admin_units_seasonal.rds")
  
  admin_matches <- admin_match(admin_unit = admin_unit, country = country,
                               admin_units_seasonal = admin_units_seasonal)
  
  if(admin_matches == 0){
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  } else {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]
  }
  
  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -params$sigma2/2, sdlog = sqrt(params$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -params$sigma2/2, sdlog = sqrt(params$sigma2))
  }
  
  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (params$max_age/365)+1
  
  #--------------------------------------------------------------
  ## Unused Default Parameters
  #--------------------------------------------------------------
  #Interventions
  if(params$itn_set == 0){
    params$max_itn_cov <- 0
    params$itn_decay_daily <- rep(0,params$n_days)
    params$itn_eff_cov_daily <- rep(0,params$n_days)
    params$mean_itn_decay <- 0
    params$d_itn0 <- 0
    params$r_itn0 <- 0
    params$r_itn1 <- 0
    
  }

  
  params$smc_mask <- array(0, dim = c(na,nh,num_int))
  if(params$smc_set == 0){
    params$max_smc_cov <- 0
    params$eff_smc_prop <- rep(0,params$n_days)
    params$P_smc <- rep(0,params$n_days)
    params$alpha_smc <- rep(0,params$n_ts)
    params$rel_c_days <- rep(1,params$n_days)
  }
  ##SMC mask must be defined in set_equilibrium, because num_int must be defined
  if(params$smc_set == 1){
    params$smc_mask[which(params$age_vector %in% params$smc_age),1:nh,3:num_int] <- 1  #Produce an array which is 1 for compartments receiving SMC, else 0. 
  }

  
  params$equilibrium_set <- 1
  
  #Climate
  if(is.null(params$daily_rain_input)) params$daily_rain_input <- rep(1,params$n_days)
  if(is.null(params$daily_temp)) params$daily_temp <- rep(1,params$n_days)
  


  
  ## collate init
  res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_Y = Y_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq, init_Iv = Iv_eq, init_Sv = Sv_eq,
              init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              omega = omega, foi_age = foi_age, rel_foi = rel_foi,
              K0 = K0, mv0 = mv0, na = na, nh = nh, ni = num_int, x_I = x_I,
              FOI_eq = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05, upto_age02 = upto_age02, 
              upto_age10 = upto_age10, ssa0 = ssa0, ssa1 = ssa1,
              ssa2 = ssa2, ssa3 = ssa3, ssb1 = ssb1, ssb2 = ssb2, ssb3 = ssb3,
              theta_c = theta_c, age = age_vector*params$DY, ft = ft, FOIv_eq = FOIv_eq,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)
  
  res <- append(res,params)
  
  return(res)
}


get_seasonality_ff_params <- function(admin_unit = NULL, country = NULL){
  if(!(is.null(country) | is.character(country))) stop("country specified is not character string")
  if(!(is.null(admin_unit) | is.character(admin_unit))) stop("admin_unit specified is not character string")
  
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")
  
  
  admin_units_seasonal <- load_file("admin_units_seasonal.rds")
  
  admin_matches <- admin_match(admin_unit = admin_unit, country = country,
                               admin_units_seasonal = admin_units_seasonal)
  
  if(admin_matches == 0){
    ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  } else {
    ssa0 <- admin_units_seasonal$a0[admin_matches]
    ssa1 <- admin_units_seasonal$a1[admin_matches]
    ssa2 <- admin_units_seasonal$a2[admin_matches]
    ssa3 <- admin_units_seasonal$a3[admin_matches]
    ssb1 <- admin_units_seasonal$b1[admin_matches]
    ssb2 <- admin_units_seasonal$b2[admin_matches]
    ssb3 <- admin_units_seasonal$b3[admin_matches]
    theta_c <- admin_units_seasonal$theta_c[admin_matches]
  }
  output <- list(g0 = ssa0,
                 g = c(ssa1, ssa2, ssa3),
                 h = c(ssb1, ssb2, ssb3),
                 theta_c = theta_c)
  return(output)
}

get_theta2 <- function(season_params, n_days){
  ssa0 <- season_params$g0
  ssa1 <- season_params$g[1]
  ssa2 <- season_params$g[2]
  ssa3 <- season_params$g[3]
  ssb1 <- season_params$h[1]
  ssb2 <- season_params$h[2]
  ssb3 <- season_params$h[3]
  theta_c <- season_params$theta_c
  time <- 1:n_days
  theta2 <- (ssa0+
               ssa1*cos(2*pi*time/365)+
               ssa2*cos(2*2*pi*time/365)+
               ssa3*cos(3*2*pi*time/365)+
               ssb1*sin(2*pi*time/365)+
               ssb2*sin(2*2*pi*time/365)+ 
               ssb3*sin(3*2*pi*time/365) ) /theta_c
  theta2[theta2 < 0.001] <- 0.001 #Set a floor of 0.001
  return(theta2)
}


get_specific_outputs <- function(ICDMM_output,
                                 ICDMM_params = NULL,
                                 ages_tot,
                                 het_brackets_tot,
                                 int_cats_tot,
                                 ages = 1:ages_tot,
                                 het_brackets_run = 1:het_brackets_tot,
                                 int_cats = 1:int_cats_tot,
                                 compartment = "all",
                                 combine_by = "sum") {
  if(!is.null(ICDMM_params)){
    ages_tot <- length(ICDMM_params$age_vector)
    het_brackets_tot <- ICDMM_params$het_brackets
    int_cats_tot <- ICDMM_params$num_int
  }
  empty_array <- array(data = NA,
                       dim = c(ages_tot, het_brackets_tot, int_cats_tot)) %>% melt() %>% dplyr::select(-value)
  colnames(empty_array) <- c("age", "het_bracket", "int_cat")
  age_indicies <- which(empty_array$age %in% ages)
  het_indicies <- which(empty_array$het_bracket %in% het_brackets_run)
  int_indicies <- which(empty_array$int_cat %in% int_cats)
  inter <- Reduce(intersect, list(age_indicies, het_indicies, int_indicies))
  ref_indicies <- paste0(".", inter-1)
  ref_indicies[ref_indicies == ".0"] <- ""
  if (compartment == "all") {
    comp = c("S", "T", "D", "A", "U", "P")
  } else if(compartment == "clin_inc") {
    comp = c("clin_inc_out")
  } else {
    comp = compartment
  }
  colnames_of_interest <- apply(expand.grid(comp, ref_indicies), 1, paste, collapse ="")
  if(combine_by == "sum"){
    spec_output <- rowSums(ICDMM_output[,colnames_of_interest])
  } else if (combine_by == "mean"){
    spec_output <- rowMeans(ICDMM_output[,colnames_of_interest])
  }
  return(spec_output)
}

get_peak_season <- function(ICDMM_output){
  yr1 <- ICDMM_output[1:365,]
  peak <- which(yr1$inc == max(yr1$inc))
  return(peak)
}

#----------------  MALARIASIMULATION-LIKE WRAPPER  -------------------------------
run_simulation <- function(gen, params){
  n_ts <- params$n_days * params$tsd #Total number of time steps to run simulation
  mod <- gen$new(pars = params,
                           time = 1,
                           n_particles = 1, #Deterministic simulation, so only one particle required  
                           n_threads = 4L,
                           seed = 1L) 
  out <- array(NA, dim = c(mod$info()$len, 1, n_ts)) #Initialize output array. Second dimension is set to 1 because model is deterministic
  for (t in seq_len(n_ts)) {
    out[ , , t] <- mod$run(t)
  }
  out <- t(drop(out)) #Transpose output matrix and make it 2D
  
  #Generate column names for output matrix
  index <- vector()
  for(i in 1:length(mod$info()$index)){
    index <- c(index, rep(names(mod$info()$index[i]),length(mod$info()$index[[i]])))
  }
  colnames(out) <- index
  
  ##Convert to dataframe and plot
  out_df <- data.frame(out[-1,])
  out_df <- out_df[!duplicated(out_df$day),] #Shrink output size by only taking the first timestep per day
  return(out_df)
}



#' @title Preset parameters for the DHA-PQP drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the DHA-PQP drug (dihydroartemisinin-piperaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.09434, drug_prophylaxis_shape: 4.4, drug_prophylaxis_scale: 28.1
#' @export
DHA_PQP_params <- c(.95, 0.09434, 4.4, 28.1)

#' @title Preset parameters for the AL drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the AL drug (artemether-lumefantrine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.05094, drug_prophylaxis_shape: 11.3, drug_prophylaxis_scale: 10.6
#' @export
AL_params <- c(.95, 0.05094, 11.3, 10.6)

#' @title Preset parameters for the SP-AQ drug
#' @details Use a vector of preset parameters for the SP-AQ drug (sulphadoxine-pyrimethamine and amodiaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.9, drug_rel_c: 0.32, drug_prophylaxis_shape: 4.3, drug_prophylaxis_scale: 38.1
#' @export
SP_AQ_params <- c(0.9, 0.32, 4.3, 38.1)


set_smc <- function(params,
                    drug_params = SP_AQ_params,
                    smc_cov = 0,
                    smc_clearance_lag = 5, #No. of days between receipt of smc and infectious clearance (default = 1/rT)
                    smc_days = NULL,
                    smc_lower_age = 0.25,
                    smc_upper_age = 5,
                    drug_efficacy = NULL,
                    drug_rel_c = NULL,
                    shape_smc = NULL,
                    scale_smc = NULL){
  
  if(is.null(smc_days)){
    stop(message("smc_days must provided"))
  }
  if(params$equilibrium_set == 1) stop(message("Equilbrium must be set last"))
  #Define specific smc parameters
  if(is.null(drug_efficacy)){drug_efficacy <- drug_params[1]}
  if(is.null(drug_rel_c)){drug_rel_c <- drug_params[2]}
  if(is.null(shape_smc)) {shape_smc <- drug_params[3]}
  if(is.null(scale_smc)) {scale_smc <- drug_params[4]}
  params$drug_efficacy <- drug_efficacy
  #Check inputs are sensible
  if(max(smc_cov) > 1) stop(message("daily_smc_cov must be < 1"))
  if(!(length(smc_cov) == 1 | length(smc_cov) == length(smc_days))) stop(message("length(smc_cov) must be either 1 or length(smc_days)"))
  if(max(smc_days) > params$n_days){
    warning(message("SMC_days which exceed n_days have been removed."))
    smc_days <- smc_days[smc_days <= params$n_days]
  }
  ##Calculate 'effective' proportion of SMC coverage
  max_smc_cov <- max(smc_cov)
  if(length(smc_cov) == 1) smc_cov <- rep(smc_cov,length(smc_days))
  params$max_smc_cov <- max_smc_cov
  
  #Calculate daily probability of prophylactic protection
  smc.dates.df <- data.frame(smc_days_join = smc_days, smc_days = smc_days, smc_cov = smc_cov)
  t <- 1:params$n_days
  smc.df <- data.frame(t = t, t_eff = 0) 
  smc.df <- merge(smc.df,smc.dates.df,by.x="t",by.y="smc_days_join",all.x=TRUE)
  
  smc.df[c("smc_days", "smc_cov")] <- lapply(smc.df[c("smc_days", "smc_cov")], function(x) {
    last <- NA
    sapply(x, function(y) ifelse(is.na(y), last, last <<- y))
  }) #Replace NAs with closest previous non-NA value
  smc.df$smc_cov[is.na(smc.df$smc_cov)] <- 0
  
  smc.df$t_eff = t - smc.df$smc_days #Days since last SMC administration
  P_smc <- exp(-((smc.df$t_eff / scale_smc) ^ shape_smc)) #Weibull distribution
  P_smc[is.na(P_smc)] <- 0
  params$P_smc <- P_smc
  
  params$eff_smc_prop <- smc.df$smc_cov / max_smc_cov
  ##The full SMC mask is defined in set_equilibrium (num_int has not been defined yet). 
  params$smc_age <- params$age_vector[params$age_vector >= smc_lower_age & params$age_vector < smc_upper_age]

  #Schedule infection clearance
  params$alpha_smc <- rep(0,params$n_days*params$tsd)
  params$alpha_smc[(smc_days + smc_clearance_lag)*params$tsd] <- drug_efficacy
  
  #The in the time between smc treatment and infection clearance - existing infections are less infectious by a factor of #smc_rel_c
  rel_c_days <- rep(1,params$n_days)
  rel_c_days[as.vector(outer(smc_days, 0:(smc_clearance_lag-1), "+"))] <- drug_rel_c
  params$rel_c_days <- rel_c_days #Time vector. Days where SMC influences infectivity are set to SMC_rel_c. Else 1.
  params$smc_set <- 1
  return(params)
}

set_itn <- function(params,
                    distribution_type = "discrete", #Assume discrete distributions or continuous? 
                    itn_days = NULL,
                    itn_cov = NULL,
                    gamman =   2.64 * 365, #ITN half life
                    itn_retention_time =   3 * 365, #Average number of days a net is kept for
                    d_itn0 = 0.41,
                    r_itn0 = 0.56,
                    r_itn1 = 0.24,
                    ...){
  n_days <- params$n_days
  
  params$d_itn0 = d_itn0
  params$r_itn0 = r_itn0
  params$r_itn1 = r_itn1
  
  if(params$equilibrium_set == 1) stop(message("Equilbrium must be set last"))
  #Check if inputs are sensible
  if(distribution_type != "discrete" & distribution_type != "continuous") stop (message("distribution_type must be either 'continuous' or 'discrete'"))
  if(max(itn_cov) > 1) stop(message("daily_itn_cov must be < 1"))
  
  if(distribution_type == "discrete"){
    usage_mat <- get_itn_usage_mat(itn_days,itn_cov,retention,n_days)
    daily_itn_cov <- get_daily_itn_cov(usage_mat) 
    
    decay_mat <- get_itn_decay_mat(itn_days,itn_cov,gamman,n_days)
    daily_itn_decay <- get_daily_itn_decay(usage_mat, decay_mat)
    params$itn_decay_daily <- daily_itn_decay
    params$max_itn_cov <- max(daily_itn_cov)
    params$itn_eff_cov_daily <- daily_itn_cov / max(daily_itn_cov)
  }
  

  # if(distribution_type == "continuous"){
  #   mean_itn_decay <- mean(exp(-((1:itn_interval)%%itn_interval) * itn_loss))  
  #   params$itn_decay <- rep(mean_itn_decay, n_days)  
  # }
  # 
  #                   
  #                   
  # ##Calculate 'effective' proportion of ITN coverage
  # if(length(daily_itn_cov == 1)) daily_itn_cov <- rep(daily_itn_cov,n_days)
  # params$max_itn_cov <- max_itn_cov
  # params$eff_itn_prop <- daily_itn_cov / max_itn_cov
  
  #Allow overwriting of ITN default parameters 
  extra_args <- list(...)
  for(name in names(extra_args)){
    params[[name]] <- extra_args[[name]]
  }
  

  params$itn_set <- 1
  return(params)
}



set_climate_parameters <- function(params,
                                   mu_tempsens = TRUE,
                                   EIP_tempsens = TRUE,
                                   daily_temp = NULL,
                                   daily_rain = NULL){
  
  #State preferred temperature sensitivity
  if(!(is.null(daily_rain) | length(daily_rain) >= n_days)) stop(message("n_days must be <= length(daily_rain)"))
  if(!(is.null(daily_temp) | length(daily_temp) >= n_days)) stop(message("n_days must be <= length(daily_temp)"))
  
  if(is.null(daily_temp)){
    mu_tempsens <- FALSE
    EIP_tempsens <- FALSE
  }
  params$mu_tempsens <- as.numeric(mu_tempsens)
  params$EIP_tempsens <- as.numeric(EIP_tempsens)
  
  if(!is.null(daily_rain)){
    rain <- daily_rain[1:n_days] 
    for(t in tau:(n_days-1)){
      exp <- vector(length = t+1)
      for(t_dash in 1:(t+1)){
        exp[t_dash] <- exp((-t+t_dash)/tau) * daily_rain[t_dash]
      }
      rain[t] <- (1/(tau*(1-exp(-(t/tau)))))*sum(exp)
      params$daily_rain_input <- rain
    }
    
    if(params$seasonality_set == 1) warning(message("daily_rain has overriden seasonality setting"))
  }
  params$set_climate <- 1
  return(params)
}

set_seasonality <- function(params, admin_unit = NULL, country = NULL){
  if(!is.null(params$daily_rain_input)) warning("Seasonality profile has replaced daily_rain input")
  if(params$equilibrium_set == 1) stop(message("Equilbrium must be set last"))
  season_params <- get_seasonality_ff_params(admin_unit = admin_unit, country = country) #Get Fourier Function parameters
  theta2 <- get_theta2(season_params, params$n_days)
  params$daily_rain_input <- theta2
  params$seasonality_set <- 1
  return(params)
}

get_itn_usage_mat <- function(itn_days,itn_cov,retention,n_days){
  n_dists <- length(itn_cov) #Number of ITN distribution events
  
  ##Matrix of the population prop. using a net from distribution event i (col number) on day t (row number)
  #Final column is those with no bednet
  itn_mat <- matrix(nrow = n_days, ncol = n_dists+1) 
  itn_mat[1,] <- c(rep(0,n_dists),1) #Initialise no bednets
  for(i in 2:n_days){
    #People discard nets
    itn_mat[i,1:n_dists] <- itn_mat[(i-1),(1:n_dists)]*exp(-1/retention)
    itn_mat[i,(n_dists+1)] <- 1 - sum(itn_mat[i,(1:n_dists)])
    
    #New nets are distributed on itn_days
    if(i %in% itn_days){
      dist_index <- which(itn_days == i) #Which distribution event
      new_nets <- itn_cov[dist_index]
      itn_mat[i,dist_index] <- new_nets
      itn_mat[i,-dist_index] <- itn_mat[i,-dist_index] - new_nets*itn_mat[i,-dist_index]
    }
  }
  return(itn_mat)
}

get_itn_decay_mat <- function(itn_days,itn_cov,gamman,n_days){
  n_dists <- length(itn_cov) #Number of ITN distribution events
  decay_mat <- matrix(0,nrow = n_days, ncol = n_dists+1)
  
  for(i in 1:n_dists){
    dist_day <- itn_days[i]
    t <- 1:(n_days - dist_day + 1)
    decay_mat[dist_day:n_days,i] <- exp(-t/gamman)
  }
  return(decay_mat)
}


get_daily_itn_cov <- function(itn_mat){
  n_dists <- ncol(itn_mat) - 1
  return(rowSums(itn_mat[,1:n_dists]))
}

get_daily_itn_decay <- function(itn_mat, decay_mat){
  n_dists <- ncol(itn_mat) - 1
  prop_net_cat <- itn_mat[,1:n_dists] / rowSums(itn_mat[,1:n_dists]) #Proportion of pop. in each net category
  prop_net_cat[is.nan(prop_net_cat)] <- 0 #Dividing by zero causes NaNs
  mean_decay <- rowSums(prop_net_cat * decay_mat[,1:n_dists])
  return(mean_decay)
}

