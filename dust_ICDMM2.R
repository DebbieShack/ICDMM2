## ICDMM2 DISCRETE DETERMINISTIC MODEL  ##

## DECLARE TIME STEPS
dt <- user()
# initial(time) <- 0
# update(time) <- (step + 1)*dt
initial(day) <- 1
update(day) <- round((step + 1)*dt) + 1
initial(iteration) <- 0
update(iteration) <- step + 1
n_days <- user() #Number of days
n_ts <- user() #Number of timesteps
###DELARE DESIRED MODEL CLIMATE SENSITIVITY
EIP_tempsens <- user() #Is EIP temperature-sensitive? Default = 0.
mu_tempsens <- user() #Is mosquito mortality temperature-sensitive? Default = 0. 

dim(daily_temp) <- n_days
daily_temp[] <- user() 
initial(temp) <- 25 
update(temp) <- daily_temp[day] 

dim(daily_rain_input) <- n_days
daily_rain_input[] <- user() 
#rainfall_input_type <- user()
# initial(rain) <- 0
# update(rain) <- daily_rain[day]

## MODEL VARIABLES
##------------------------------------------------------------------------------

na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------

# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

# fitted parameters for human compartments:
eta <- user() # death rate for exponential population distribution
age_rate[] <- user() # rate at which humans move through age categories
dim(age_rate) <- na
het_wt[] <- user() # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
rA <- user() # rate of movement from A -> U
rT <- user() # rate of treatment working: T -> P
rD <- user() #  rate from D -> A
rU <- user() # rate of clearance of subpatent infection U -> S
rP <- user() # rate at which prophylaxis wears off P -> S

################################## INITIALISE HUMAN COMPARTMENTS ######################################
# S - SUSCEPTIBLE
init_S[,,] <- user()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

# T- SUCCESSFULLY TREATED
init_T[,,] <- user()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

# D - CLEAR DISEASE
init_D[,,] <- user()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

# A - ASYMPTOMATIC DISEASE
init_A[,,] <- user()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

# U - SUBPATENT DISEASE
init_U[,,] <- user()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

# P - PROPHYLAXIS
init_P[,,] <- user()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

################################## DEFINE HUMAN EVENTS ######################################
#Events are organised by the compartment from which an individual leaves
## Births ##
dim(births) <- c(1,nh,num_int)
births[1,,] <- H*dt*cov[k]*eta*het_wt[j]

# S ##

##Individual transition rates
dim(ST_rate) <- c(na,nh,num_int)
dim(SD_rate) <- c(na,nh,num_int)
dim(SA_rate) <- c(na,nh,num_int)

ST_rate[,,] <- ft*phi[i,j,k]*FOI_smc[i,j,k]
SD_rate[,,] <- (1-ft)*phi[i,j,k]*FOI_smc[i,j,k]
SA_rate[,,] <- (1-phi[i,j,k])*FOI_smc[i,j,k]


##Individual transitions numbers
dim(ST_trans) <- c(na,nh,num_int)
dim(SD_trans) <- c(na,nh,num_int)
dim(SA_trans) <- c(na,nh,num_int)
dim(S_death) <- c(na,nh,num_int)
dim(S_age) <- c(na,nh,num_int)

ST_trans[,,] <- if(S[i,j,k]*dt*ST_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*ST_rate[i,j,k]
SD_trans[,,] <- if(S[i,j,k]*dt*SD_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SD_rate[i,j,k]
SA_trans[,,] <- if(S[i,j,k]*dt*SA_rate[i,j,k] < 0) 0 else S[i,j,k]*dt*SA_rate[i,j,k]
S_death[,,] <- if(S[i,j,k]*dt*eta < 0) 0 else S[i,j,k]*dt*eta
S_age[,,] <- if(S[i,j,k]*dt*age_rate[i] < 0) 0 else S[i,j,k]*dt*age_rate[i]

## T ##
dim(TP_trans) <- c(na,nh,num_int)
dim(T_death) <- c(na,nh,num_int)
dim(T_age) <- c(na,nh,num_int)

TP_trans[,,] <- if(T[i,j,k]*dt*rT < 0) 0 else T[i,j,k]*dt*rT
T_death[,,] <- if(T[i,j,k]*dt*eta < 0) 0 else T[i,j,k]*dt*eta
T_age[,,] <- if(T[i,j,k]*dt*age_rate[i] < 0) 0 else T[i,j,k]*dt*age_rate[i]

# D ##
dim(DA_trans) <- c(na,nh,num_int)
dim(DS_trans) <- c(na,nh,num_int)
dim(D_death) <- c(na,nh,num_int)
dim(D_age) <- c(na,nh,num_int)

DA_trans[,,] <- if(D[i,j,k]*dt*rD < 0) 0 else D[i,j,k]*dt*rD 
DS_trans[,,] <- if(D[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else D[i,j,k]*alpha_smc_array[i,j,k]#smc Infection Clearing
D_death[,,] <- if(D[i,j,k]*dt*eta < 0) 0 else D[i,j,k]*dt*eta
D_age[,,] <- if(D[i,j,k]*dt*age_rate[i] < 0) 0 else D[i,j,k]*dt*age_rate[i]

## A ##
dim(AT_rate) <- c(na,nh,num_int)
dim(AD_rate) <- c(na,nh,num_int)

AT_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * ft
AD_rate[,,] <- phi[i,j,k] * FOI_smc[i,j,k] * (1-ft)

dim(AU_trans) <- c(na,nh,num_int)
dim(AT_trans) <- c(na,nh,num_int)  
dim(AD_trans) <- c(na,nh,num_int)
dim(AS_trans) <- c(na,nh,num_int) #SMC Infection Clearing
dim(A_age) <- c(na,nh,num_int)
dim(A_death) <- c(na,nh,num_int)

AU_trans[,,] <- if(A[i,j,k]*dt*rA < 0) 0 else A[i,j,k]*dt*rA
AT_trans[,,] <- if(A[i,j,k]*dt*AT_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AT_rate[i,j,k]
AD_trans[,,] <- if(A[i,j,k]*dt*AD_rate[i,j,k] < 0) 0 else A[i,j,k]*dt*AD_rate[i,j,k]
AS_trans[,,] <- if(A[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else A[i,j,k]*alpha_smc_array[i,j,k]#smc Infection Clearing
A_death[,,] <- if(A[i,j,k]*dt*eta < 0) 0 else A[i,j,k]*dt*eta
A_age[,,] <- if(A[i,j,k]*dt*age_rate[i] < 0) 0 else A[i,j,k]*dt*age_rate[i]

## U ##
dim(UA_rate) <- c(na,nh,num_int)
dim(UD_rate) <- c(na,nh,num_int)
dim(UT_rate) <- c(na,nh,num_int)

UA_rate[,,] <- (1-phi[i,j,k]) * FOI_smc[i,j,k]
UD_rate[,,] <- phi[i,j,k] * (1-ft) * FOI_smc[i,j,k]
UT_rate[,,] <- phi[i,j,k] * ft * FOI_smc[i,j,k]

dim(UA_trans) <- c(na,nh,num_int)
dim(UD_trans) <- c(na,nh,num_int)
dim(UT_trans) <- c(na,nh,num_int)
dim(US_trans) <- c(na,nh,num_int) 
dim(US_trans_smc) <- c(na,nh,num_int) #SMC Infection Clearing
dim(U_age) <- c(na,nh,num_int)
dim(U_death) <- c(na,nh,num_int)

US_trans[,,] <- if(U[i,j,k]*dt*rU < 0) 0 else U[i,j,k]*dt*rU
UA_trans[,,] <- if(U[i,j,k]*dt*UA_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UA_rate[i,j,k]
UD_trans[,,] <- if(U[i,j,k]*dt*UD_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UD_rate[i,j,k]
UT_trans[,,] <- if(U[i,j,k]*dt*UT_rate[i,j,k] < 0) 0 else U[i,j,k]*dt*UT_rate[i,j,k]
US_trans_smc[,,] <- if(U[i,j,k]*alpha_smc_array[i,j,k] < 0) 0 else U[i,j,k]*alpha_smc_array[i,j,k]#SMC Infection Clearing
U_age[,,] <- if(U[i,j,k]*dt*age_rate[i] < 0) 0 else U[i,j,k]*dt*age_rate[i]
U_death[,,] <- if(U[i,j,k]*dt*eta < 0) 0 else U[i,j,k]*dt*eta

## P ##
dim(PS_trans) <- c(na,nh,num_int)
dim(P_death) <- c(na,nh,num_int)
dim(P_age) <- c(na,nh,num_int)

PS_trans[,,] <-if(P[i,j,k]*dt*rP < 0) 0 else P[i,j,k]*dt*rP
P_death[,,] <- if(P[i,j,k]*dt*eta < 0) 0 else P[i,j,k]*dt*eta
P_age[,,] <- if(P[i,j,k]*dt*age_rate[i] < 0) 0 else P[i,j,k]*dt*age_rate[i]




################################## TRANSITIONS ######################################
update(S[1,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_smc[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + births[1,j,k]
update(S[2:na,1:nh,1:num_int]) <- S[i,j,k] + PS_trans[i,j,k] + US_trans[i,j,k] + US_trans_smc[i,j,k] + AS_trans[i,j,k] + DS_trans[i,j,k] - ST_trans[i,j,k] - SD_trans[i,j,k] - SA_trans[i,j,k] - S_death[i,j,k] - S_age[i,j,k] + S_age[i-1,j,k]

update(T[1,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k]
update(T[2:na,1:nh,1:num_int]) <- T[i,j,k] + AT_trans[i,j,k] + ST_trans[i,j,k] + UT_trans[i,j,k] - TP_trans[i,j,k] - T_age[i,j,k] - T_death[i,j,k] + T_age[i-1,j,k]

update(D[1,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k]
update(D[2:na,1:nh,1:num_int]) <- D[i,j,k] + SD_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k] - DA_trans[i,j,k] - DS_trans[i,j,k] - D_death[i,j,k] - D_age[i,j,k] + D_age[i-1,j,k]

update(A[1,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k]
update(A[2:na,1:nh,1:num_int]) <- A[i,j,k] + SA_trans[i,j,k] + DA_trans[i,j,k] + UA_trans[i,j,k] - AT_trans[i,j,k] - AD_trans[i,j,k] - AU_trans[i,j,k] - AS_trans[i,j,k] - A_death[i,j,k] - A_age[i,j,k] + A_age[i-1,j,k]

update(U[1,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_smc[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k]
update(U[2:na,1:nh,1:num_int]) <- U[i,j,k] + AU_trans[i,j,k] - UD_trans[i,j,k] - UT_trans[i,j,k] - US_trans[i,j,k] - US_trans_smc[i,j,k] - UA_trans[i,j,k] - U_age[i,j,k] - U_death[i,j,k] + U_age[i-1,j,k]

update(P[1,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k]
update(P[2:na,1:nh,1:num_int]) <- P[i,j,k] + TP_trans[i,j,k] - PS_trans[i,j,k] - P_death[i,j,k] - P_age[i,j,k] + P_age[i-1,j,k]

dim(clin_inc) <- c(na,nh,num_int)
clin_inc[,,] <- ST_trans[i,j,k] + SD_trans[i,j,k] + AT_trans[i,j,k] + AD_trans[i,j,k] + UD_trans[i,j,k]  + UT_trans[i,j,k] 

dim(clin_inc_out) <- c(na,nh,num_int)
initial(clin_inc_out[,,]) <- init_P[i,j,k]
update(clin_inc_out[,,]) <- clin_inc[i,j,k] / dt

##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

# fitted immunity parameters:
dCM <- user() # decay of maternal immunity
uCA <- user() # scale parameter (see Supplementary mats. 3.1.2)
dCA <- user() # decay for clinical immunity
dB <- user() # decay for infection blocking immunity
uB <- user() # scale param for IB immunity
dID <- user() # decay for detection immunity
uD <- user() # scale param for ID immunity
x_I[] <- user() # intermediate variable for calculating immunity functions
dim(x_I) <- na
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user() # factor calculated in equilibrium solution
PM <- user() # immunity constant

# ICM - maternally acquired immunity
init_ICM[,,] <- user()
dim(init_ICM) <- c(na,nh,num_int)
initial(ICM[,,]) <- init_ICM[i,j,k]
dim(ICM) <- c(na,nh,num_int)
dim(init_ICM_pre) <- c(nh,num_int)
init_ICM_pre[1:nh,1:num_int] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))

update(ICM[1, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i])
update(ICM[2:na, 1:nh, 1:num_int]) <- ICM[i,j,k] + dt*(-1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i])

# ICA - exposure driven immunity
init_ICA[,,] <- user()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

update(ICA[1, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i])
update(ICA[2:na, 1:nh, 1:num_int]) <- ICA[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i])

dim(FOI_smc_out) <- c(na,nh,num_int)
initial(FOI_smc_out[,,]) <- init_ICA[i,j,k]
update(FOI_smc_out[,,]) <- FOI_smc[i,j,k]

dim(P_smc_out) <- c(na,nh,num_int)
initial(P_smc_out[,,]) <- init_ICA[i,j,k]
update(P_smc_out[,,]) <- P_smc_array[i,j,k]
# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,num_int)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,num_int)
phi[1:na,1:nh,1:num_int] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)


##SMC


# dim(phi_out) <- c(na,nh,num_int)
# initial(phi_out[,,]) <- init_ICA[i,j,k]
# update(phi_out[,,]) <- phi[i,j,k]

# IB - infection blocking immunity
init_IB[,,] <- user()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

update(IB[1, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i])
update(IB[2:na, 1:nh, 1:num_int]) <- IB[i,j,k] + dt*(EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i])

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh,num_int)
b[1:na, 1:nh, 1:num_int] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)

# detection immunity
init_ID[,,] <- user()
dim(init_ID) <- c(na,nh,num_int)
initial(ID[,,]) <- init_ID[i,j,k]
dim(ID) <- c(na,nh,num_int)

update(ID[1, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i])
update(ID[2:na, 1:nh, 1:num_int]) <- ID[i,j,k] + dt*(FOI_smc[i,j,k]/(FOI_smc[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i])

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- user()
fD0 <- user()
gammaD <- user()
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:na] <- 1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,num_int)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,num_int)
FOI_lag[1:na, 1:nh, 1:num_int] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection.
lag_rates <- user()

dim(FOI_eq) <- c(na,nh)
FOI_eq[,] <- user()

dim(FOI_XL) <- c(na,nh,num_int,lag_rates)
initial(FOI_XL[,,,]) <- FOI_eq[i,j]

update(FOI_XL[,,,1]) <- FOI_XL[i,j,k,1] + dt*((lag_rates/dE)*FOI_lag[i,j,k] - (lag_rates/dE)*FOI_XL[i,j,k,1])
update(FOI_XL[,,,2:lag_rates]) <- FOI_XL[i,j,k,l] + dt*((lag_rates/dE)*FOI_XL[i,j,k,l-1] - (lag_rates/dE)*FOI_XL[i,j,k,l])

dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- FOI_XL[i,j,k,lag_rates]

# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
dim(EIR) <- c(na,nh,num_int)
EIR[,,] <- av_human[k] * rel_foi[j] * foi_age[i] * Iv/omega

dim(EIR_out) <- c(na,nh,num_int)
initial(EIR_out[,,]) <- init_ID[i,j,k]
update(EIR_out[,,]) <- EIR[i,j,k] 

initial(Ivout) <- 0
initial(omega_out) <- 0
update(Ivout) <- Iv
update(omega_out) <- omega

dim(rel_foi_out) <- nh
initial(rel_foi_out[]) <- rel_foi[i]
update(rel_foi_out[]) <- rel_foi[i]

##------------------------------------------------------------------------------
#####################
## MOSQUITO STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# Sv - Susceptible mosquitoes
# Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
# Iv - Infectious mosquitoes

# initial state values:
init_Sv <- user()
init_Ev <- user()
init_Iv <- user()
initial(Sv) <- init_Sv * mv0
initial(Ev) <- init_Ev * mv0
initial(Iv) <- init_Iv * mv0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1 <- user() # fitted value of gamma1 characterises cA function
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
lag_ratesMos <- lag_rates
FOIv_eq <- user()
dim(FOIv) <- lag_ratesMos
initial(FOIv[]) <- FOIv_eq*delayGam/lag_ratesMos

update(FOIv[1]) <- FOIv[1] + dt*(lag_FOIv - (lag_ratesMos/delayGam)*FOIv[1])
update(FOIv[2:lag_ratesMos]) <- FOIv[i] + dt*((lag_ratesMos/delayGam)*FOIv[i-1] - (lag_ratesMos/delayGam)*FOIv[i])


dim(FOIvijk) <- c(na,nh,num_int)
omega <- user() #normalising constant for biting rates
FOIvijk[1:na, 1:nh, 1:num_int] <- ((cT*smc_rel_c_mask[i,j,k]*T[i,j,k] + cD*smc_rel_c_mask[i,j,k]*D[i,j,k] + cA[i,j,k]*smc_rel_c_mask[i,j,k]*A[i,j,k] + cU*smc_rel_c_mask[i,j,k]*U[i,j,k])/H) *
  rel_foi[j] * av_mosq[k]*foi_age[i]/omega ## For discrete human compartments
lag_FOIv=sum(FOIvijk)

initial(lag_FOIv_out) <- 0
update(lag_FOIv_out) <- lag_FOIv

ince <- FOIv[lag_ratesMos] * lag_ratesMos/delayGam * Sv

initial(ince_delay[]) <- FOIv_eq*init_Sv*mv0*delayMos_use/lag_ratesMos
dim(ince_delay) <- lag_ratesMos

update(ince_delay[1]) <- ince_delay[1] + dt*(ince - (lag_ratesMos/delayMos_use)*ince_delay[1])
update(ince_delay[2:lag_ratesMos]) <- ince_delay[i] + dt*((lag_ratesMos/delayMos_use)*ince_delay[i-1] -
                                                            (lag_ratesMos/delayMos_use)*ince_delay[i])

incv <- ince_delay[lag_ratesMos]*lag_ratesMos/delayMos_use *surv

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # Lag from parasites to infectious gametocytes
delayMos <- user() # Extrinsic incubation period.
delayMos_use <- if(EIP_tempsens == 1) 111 / (temp - 16) else delayMos
# Number of mosquitoes that become infected at each time point
surv <- exp(-mu*delayMos_use)

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa <- 0.5*PL/dPL

update(Sv) <- if(Sv + dt*(-ince - mu*Sv + betaa) < 0) 0 else Sv + dt*(-ince - mu*Sv + betaa)
update(Ev) <- if(Ev + dt*(ince - incv - mu*Ev) < 0) 0 else Ev + dt*(ince - incv - mu*Ev)
update(Iv) <- if(Iv + dt*(incv - mu*Iv) < 0) 0 else Iv + dt*(incv - mu*Iv)

# Total mosquito population
mv = Sv+Ev+Iv

##------------------------------------------------------------------------------
###################
## LARVAL STATES ##
###################
##------------------------------------------------------------------------------

# Model by White et al.
# (https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-4-153)

# EL - early larval instar stage
# LL - late larval instar stage
# PL - pupal stage

# mean carrying capacity from initial mosquito density:
dLL <- user() # development time of larvae
dPL <- user() #development time of pupae
dEL <- user() #development time of early stage
muLL <- user() #daily density dep. mortality rate of larvae
muEL <- user() #daily den. dep. mortality rate of early stage
muPL <- user() #daily den. dep. mortality rate of pupae
gammaL <- user() # eff. of den. dep. on late stage relative to early stage

# fitted entomological parameters:
mv0 <- user() # initial mosquito density
mu0 <- user() # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
# p10 <- user() # prob of surviving 1 feeding cycle
# p2 <- user() #prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# Entomological variables:
Pa <- -8.3E-4*temp^2 + 0.037*temp + 0.522 #Daily mosquito survival probability. Quadratic function taken from Mordecai 2013
mu0_use <- if(mu_tempsens == 1) -log(Pa) else mu0
p10 <- exp(-mu0_use * tau1)  # probability of surviving one feeding cycle
p2 <- exp(-mu0_use * tau2)  # probability of surviving one resting cycle
eov <- betaL/mu*(exp(mu/fv)-1)
beta_larval <- eov*mu*exp(-mu/fv)/(1-exp(-mu/fv)) # Number of eggs laid per day
b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL)
lambda <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval*muLL*dEL/(2*muEL*mu0_use*dLL*(1+dPL*muPL)))
K0 <- 2*mv0*dLL*mu0_use*(1+dPL*muPL)*gammaL*(lambda+1)/(lambda/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
initial(rain_input) <- 1
update(rain_input) <- daily_rain_input[day]

KL <- K0*rain_input
fv <- 1/( tau1/(1-zbar) + tau2 ) # mosquito feeding rate (zbar from intervention param.)
mu <- -fv*log(p1*p2) # mosquito death rate

# finding equilibrium and initial values for EL, LL & PL
init_PL <- user()
initial(PL) <- init_PL
init_LL <- user()
initial(LL) <- init_LL
init_EL <- user()
initial(EL) <- init_EL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
update(EL) <- if(EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL) < 0) 0 else (EL + dt*(beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL))
# egg hatching - den. dep. mortality - maturing larvae
update(LL) <- if(LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL) < 0) 0 else (LL + dt*(EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL))
# pupae - mortality - fully developed pupae
update(PL) <- if(PL + dt*(LL/dLL - muPL*PL - PL/dPL) < 0) 0 else (PL + dt*(LL/dLL - muPL*PL - PL/dPL))


##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ##
########################
##------------------------------------------------------------------------------

################## SMC ####################### 
#See supplementary materials of [Thompson, 2022] - https://doi.org/10.1016/S2214-109X(22)00416-8 
#drug_efficacy<- user()
max_smc_cov <- user() #Highest proportion of individuals within eligible age range receiving SMC over simulation period
dim(eff_smc_prop) <- n_days
eff_smc_prop[] <- user() #Proportion of individuals in the SMC compartment who are currently receiving each SMC dose = daily_SMC / max_SMC

##Parameters relevant to clearance of existing infections by SMC treatment
dim(alpha_smc) <- n_ts
alpha_smc[] <- user() ##Proportion of existing infections cleared by SMC. 
dim(alpha_smc_array) <-  c(na,nh,num_int)
alpha_smc_array[,,] <- smc_mask[i,j,k] * alpha_smc[iteration] * eff_smc_prop[day] #Multiply by proportion of individuals in SMC compartment currently receiving SMC

##Parameters relating to prophylaxis effect of SMC
dim(P_smc) <- n_days
P_smc[] <- user() #Probability that an individual is protected from clinical malaria at time t. Wanes according to Weibull distribution from last dose.
dim(smc_mask) <- c(na,nh,num_int)
smc_mask[,,] <- user() #Binary array. 1 for those receiving SMC treatment. Else 0. 
dim(P_smc_array) <- c(na,nh,num_int)
P_smc_array[,,] <- smc_mask[i,j,k] * P_smc[day] * eff_smc_prop[day] #* drug_efficacy
dim(FOI_smc) <- c(na,nh,num_int)
FOI_smc[,,] <- FOI[i,j,k] * (1 - P_smc_array[i,j,k])

##Parameters relating to reduced infectivity of existing infections upon SMC dose
dim(rel_c_days) <- n_days
rel_c_days[] <- user() #Equal to SMC_rel_c on SMC_days until infection clearance (default 5 days). Else = 1.
#dim(SMC_rel_c_intermediate) <- c(na,nh,num_int)
dim(smc_rel_c_mask) <- c(na,nh,num_int)
#SMC_rel_c_intermediate[,,] <- SMC_mask[i,j,k] * (1-rel_c_days[day])
smc_rel_c_mask[,,] <- 1 - (smc_mask[i,j,k] * (1-rel_c_days[day])) #Value = SMC_rel_c for those scheduled to receive SMC on rel_c_days. Else = 1. 


################## ITN #######################
# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
max_itn_cov <- user() #Highest proportion of individuals scheduled to be protected by ITNs during simulation period.
 

Q0 <- user() # proportion of anthropophagy
bites_Bed <- user() # endophagy in bed

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN
# d - probability of dying after hitting ITN
# s - probability of successful feed after hitting ITN

# The maximum (and then minimum) r and d values for ITN on day 0 before they decay
r_itn0 <- user()
d_itn0 <- user()
r_itn1 <- user()

dim(itn_decay_daily) <- n_days #Bednet insecticide decay
dim(itn_eff_cov_daily) <- n_days #Proportion of individuals in the ITN compartment who currently have ITNs.= daily_ITN / max(daily_ITN)
itn_decay_daily[] <- user()
itn_eff_cov_daily[] <- user()

d_itn <- d_itn0*itn_decay_daily[day]*itn_eff_cov_daily[day]
r_itn <- (r_itn1 + (r_itn0 - r_itn1)*itn_decay_daily[day])*itn_eff_cov_daily[day] #Bednet repellency does not decay to zero. 
s_itn <- 1 - d_itn - r_itn

################## GENERAL INTERVENTION PARAMETERS #######################
num_int <- user() # Which interventions: None, ITN only, SMC only, both

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-max_itn_cov)*(1-max_smc_cov)  # {No intervention}
cov_[2] <- max_itn_cov*(1-max_smc_cov) # 	   {ITN only}
cov_[3] <- (1-max_itn_cov)*max_smc_cov	#      {SMC only}
cov_[4] <- max_itn_cov*max_smc_cov #	   {Both itn and SMC}
cov[] <- cov_[i]
dim(cov) <- num_int
dim(cov_out) <- num_int

# probability that mosquito bites and survives for each intervention category
dim(w_) <- 4
w_[1] <- 1
w_[2] <- 1 - bites_Bed + bites_Bed*s_itn
w_[3] <- 1 
w_[4] <- 1 - bites_Bed + bites_Bed*s_itn
w[] <- w_[i]
dim(w) <- num_int

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy_) <- 4
yy_[1] <- 1
yy_[2] <- w_[2]
yy_[3] <- 1
yy_[4] <- w_[4]
yy[] <- yy_[i]
dim(yy) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z_) <- 4
z_[1] <- 0
z_[2] <- bites_Bed*r_itn
z_[3] <- 0
z_[4] <- bites_Bed*r_itn
z[] <- z_[i]
dim(z) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi) <- num_int
dim(whi) <- num_int
zhi[1:num_int] <- cov[i]*z[i]
whi[1:num_int] <- cov[i]*w[i]
zh <- sum(zhi)
wh <- sum(whi)
# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar <- Q0*zh
# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar <- 1 - Q0 + Q0*wh

# p1 is the updated p10 given that interventions are now in place:
p1 <- wbar*p10/(1-zbar*p10)
Q <- 1-(1-Q0)/wbar # updated anthropophagy given interventions
av <- fv*Q # biting rate on humans
dim(av_mosq) <- num_int
av_mosq[1:num_int] <- av*w[i]/wh # rate at which mosquitoes bite each int. cat.
dim(av_human) <- num_int
av_human[1:num_int] <- av*yy[i]/wh # biting rate on humans in each int. cat.

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
initial(Sout) <- 0
initial(Tout) <- 0
initial(Dout) <- 0
initial(Aout) <- 0
initial(Uout) <- 0
initial(Pout) <- 0

update(Sout) <- sum(S[,,])
update(Tout) <- sum(T[,,])
update(Dout) <- sum(D[,,])
update(Aout) <- sum(A[,,])
update(Uout) <- sum(U[,,])
update(Pout) <- sum(P[,,])

# Outputs for clinical incidence and prevalence on a given day
# population densities for each age category
den[] <- user()
dim(den) <- na
# index of the age vector above 59 months
age59 <- user(integer=TRUE)
# index of the age vector above 5 years
age05 <- user(integer=TRUE)

dim(prev0to59) <- c(age59,nh,num_int)
prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]
initial(prev) <- 0
update(prev) <- (sum(prev0to59[,,])/sum(den[1:age59]))/H

dim(detected) <- c(na,nh,num_int)
detected[,,] <- T[i,j,k] + D[i,j,k] + A[i,j,k]*p_det[i,j,k]
dim(detected_out) <- c(na,nh,num_int)
initial(detected_out[,,]) <- init_P[i,j,k]
update(detected_out[,,]) <- detected[i,j,k]

dim(all) <- c(na,nh,num_int)
all[,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]
dim(all_out) <- c(na,nh,num_int)
initial(all_out[,,]) <- init_P[i,j,k]
update(all_out[,,]) <- all[i,j,k]

# ##Estimating PfPR_2_10 (to match malariasimulation)
upto_age02 <- user(integer=TRUE) #index of age vector up to 2 years
upto_age10 <- user(integer=TRUE) #index of age vector up to 10 years

#
dim(detect_0_2) <- c(upto_age02,nh,num_int) #Detectable malaria infections between children aged 2-10
detect_0_2[1:upto_age02,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]

dim(detect_0_10) <- c(upto_age10,nh,num_int) #Detectable malaria infections between children aged 2-10
detect_0_10[1:upto_age10,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]

initial(n_detect_2_10) <- 0
update(n_detect_2_10) <- sum(detect_0_10) - sum(detect_0_2)


#
dim(total_0_2) <- c(upto_age02,nh,num_int)
total_0_2[1:upto_age02,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]

dim(total_0_10) <- c(upto_age10,nh,num_int)
total_0_10[1:upto_age10,,] <- S[i,j,k] + T[i,j,k] + D[i,j,k] + A[i,j,k] + U[i,j,k] + P[i,j,k]


initial(n_2_10) <- 0
update(n_2_10) <- sum(total_0_10) - sum(total_0_2)

initial(PfPR_2_10) <- 0
update(PfPR_2_10) <- (sum(detect_0_10) - sum(detect_0_2)) / (sum(total_0_10) - sum(total_0_2))

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]



initial(inc05) <- 0
initial(incUnder5) <- 0
initial(inc) <- 0
update(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
update(inc) <- sum(clin_inc[,,]) / dt
update(incUnder5) <- sum(clin_inc0to5[,,]) / dt

initial(incAdult) <- 0
update(incAdult) <- (sum(clin_inc[,,]) - sum(clin_inc0to5[,,])) / dt

# Param checking outputs
initial(mu_out) <- 0
initial(beta_larval_out) <- 0
initial(KL_out) <- 0
initial(mv_out) <- 0
initial(Q_out) <- 0
initial(d_itn_out) <- 0
initial(r_itn_out) <- 0
initial(s_itn_out) <- 0
initial(cov_out[]) <- TRUE
initial(K0_out) <- 0

update(mu_out) <- mu
update(beta_larval_out) <- beta_larval
update(KL_out) <- KL
update(mv_out) <- mv
update(Q_out) <- Q
update(d_itn_out) <- d_itn
update(r_itn_out) <- r_itn
update(s_itn_out) <- s_itn
update(cov_out[]) <- TRUE
update(K0_out) <- K0

initial(cov1) <- 0
initial(cov2) <- 0
initial(cov3) <- 0
initial(cov4) <- 0

update(cov1) <- cov[1]
update(cov2) <- cov[2]
update(cov3) <- cov[3]
update(cov4) <- cov[4]


initial(zh_out) <- 0
initial(wh_out) <- 0
update(zh_out) <- zh
update(wh_out) <- wh

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph


deaths <- sum(S_death[,,]) + sum(T_death[,,]) + sum(D_death[,,]) + sum(A_death[,,]) + sum(U_death[,,]) + sum(P_death[,,])
initial(deaths_out) <- 0
update(deaths_out) <- deaths
initial(births_out) <- 0
update(births_out) <- sum(births[,,])

initial(H_out) <- 0
update(H_out) <- H

initial(alpha_smc_out) <- 0
update(alpha_smc_out) <- alpha_smc[iteration]