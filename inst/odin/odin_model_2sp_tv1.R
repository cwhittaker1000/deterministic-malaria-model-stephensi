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

# S - SUSCEPTIBLE
init_S[,,] <- user()
dim(init_S) <- c(na,nh,num_int)
initial(S[,,]) <- init_S[i,j,k]
dim(S) <- c(na,nh,num_int)

deriv(S[1, 1:nh, 1:num_int]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] +
  cov[k]*eta*H*het_wt[j] - (eta+age_rate[i])*S[i,j,k]
deriv(S[2:na, 1:nh, 1:num_int]) <- -FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] -
  (eta+age_rate[i])*S[i,j,k] + age_rate[i-1]*S[i-1,j,k]

# T- SUCCESSFULLY TREATED
init_T[,,] <- user()
dim(init_T) <- c(na,nh,num_int)
initial(T[,,]) <- init_T[i,j,k]
dim(T) <- c(na,nh,num_int)

deriv(T[1, 1:nh, 1:num_int]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k]
deriv(T[2:na, 1:nh, 1:num_int]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k] + age_rate[i-1]*T[i-1,j,k]

# D - CLEAR DISEASE
init_D[,,] <- user()
dim(init_D) <- c(na,nh,num_int)
initial(D[,,]) <- init_D[i,j,k]
dim(D) <- c(na,nh,num_int)

deriv(D[1, 1:nh, 1:num_int]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k]
deriv(D[2:na, 1:nh, 1:num_int]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k] + age_rate[i-1]*D[i-1,j,k]

# A - ASYMPTOMATIC DISEASE
init_A[,,] <- user()
dim(init_A) <- c(na,nh,num_int)
initial(A[,,]) <- init_A[i,j,k]
dim(A) <- c(na,nh,num_int)

deriv(A[1, 1:nh, 1:num_int]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k]
deriv(A[2:na, 1:nh, 1:num_int]) <- (1-phi[i,j,k])*FOI[i,j,k]*Y[i,j,k] - FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k] + age_rate[i-1]*A[i-1,j,k]

# U - SUBPATENT DISEASE
init_U[,,] <- user()
dim(init_U) <- c(na,nh,num_int)
initial(U[,,]) <- init_U[i,j,k]
dim(U) <- c(na,nh,num_int)

deriv(U[1, 1:nh, 1:num_int]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k]
deriv(U[2:na, 1:nh, 1:num_int]) <- rA*A[i,j,k] - FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k] + age_rate[i-1]*U[i-1,j,k]

# P - PROPHYLAXIS
init_P[,,] <- user()
dim(init_P) <- c(na,nh,num_int)
initial(P[,,]) <- init_P[i,j,k]
dim(P) <- c(na,nh,num_int)

deriv(P[1, 1:nh, 1:num_int]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k]
deriv(P[2:na, 1:nh, 1:num_int]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k] +
  age_rate[i-1]*P[i-1,j,k]

# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,num_int)
Y[1:na, 1:nh, 1:num_int] <- S[i,j,k]+A[i,j,k]+U[i,j,k]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,num_int)
clin_inc[1:na, 1:nh, 1:num_int] <- phi[i,j,k]*FOI[i,j,k]*Y[i,j,k]

# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph

output(all_pop) <- H

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

deriv(ICM[1, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] + (init_ICM_pre[j,k]-ICM[i,j,k])/x_I[i]
deriv(ICM[2:na, 1:nh, 1:num_int]) <- -1/dCM*ICM[i,j,k] - (ICM[i,j,k]-ICM[i-1,j,k])/x_I[i]

# ICA - exposure driven immunity
init_ICA[,,] <- user()
dim(init_ICA) <- c(na,nh,num_int)
initial(ICA[,,]) <- init_ICA[i,j,k]
dim(ICA) <- c(na,nh,num_int)

deriv(ICA[1, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] -ICA[i,j,k]/x_I[i]
deriv(ICA[2:na, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k] * uCA + 1) - 1/dCA*ICA[i,j,k] - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i]

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

# IB - infection blocking immunity
init_IB[,,] <- user()
dim(init_IB) <- c(na,nh,num_int)
initial(IB[,,]) <- init_IB[i,j,k]
dim(IB) <- c(na,nh,num_int)

deriv(IB[1, 1:nh, 1:num_int]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i]
deriv(IB[2:na, 1:nh, 1:num_int]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i]

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

deriv(ID[1, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i]
deriv(ID[2:na, 1:nh, 1:num_int]) <- FOI[i,j,k]/(FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i]

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
dim(FOI) <- c(na,nh,num_int)
FOI[,,] <- delay(FOI_lag[i,j,k],dE)

# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
# infectious mosquitoes
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()
dim(EIR) <- c(na,nh,num_int)

# CHECK - EIR for each individual species
EIR_species1[,,] <- av_human1[k] * rel_foi[j] * foi_age[i] * Iv1/omega
dim(EIR_species1) <- c(na,nh,num_int)
EIR_species2[,,] <- av_human2[k] * rel_foi[j] * foi_age[i] * Iv2/omega
dim(EIR_species2) <- c(na,nh,num_int)
EIR[,,] <- EIR_species1[i,j,k] + EIR_species2[i,j,k]

output(Ivout1) <- Iv1
output(Ivout2) <- Iv2
output(EIR_species1) <- EIR_species1
output(EIR_species2) <- EIR_species2
output(EIR) <- EIR

output(omega) <- omega
##------------------------------------------------------------------------------
##########################
## SEASONALITY FUNCTION ##
##########################
##------------------------------------------------------------------------------

# Seasonality is added into the model using a Fourier series that was fit to rainfall at every admin 1 level
pi <- user() # weird quirk, need to pass pi

# The parameters for the fourier series
ssa0 <- user()
ssa1 <- user()
ssa2 <- user()
ssa3 <- user()
ssb1 <- user()
ssb2 <- user()
ssb3 <- user()
theta_c <- user()

# Recreation of the rainfall function (used for Species 1, which are the endemic African vectors)
theta_species1 <- if(ssa0 == 0 && ssa1  == 0 && ssa2  == 0 && ssb1  == 0 && ssb2  == 0 && ssb3  == 0 && theta_c  == 0)
  1 else max((ssa0+ssa1*cos(2*pi*t/365)+ssa2*cos(2*2*pi*t/365)+ssa3*cos(3*2*pi*t/365)+ssb1*sin(2*pi*t/365)+ssb2*sin(2*2*pi*t/365)+ ssb3*sin(3*2*pi*t/365) ) /theta_c,0.001)
output(theta_species1) <- theta_species1

## Custom seasonality function for species 2 (which is An.stephensi for us)
custom_seasonality[] <- user()
dim(custom_seasonality) <- time_length
theta_species2 <- if(as.integer(t) == 0) custom_seasonality[as.integer(1)] else custom_seasonality[as.integer(t)]
output(theta_species2) <- theta_species2

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

initial(Sv1) <- init_Sv * mv0
initial(Ev1) <- init_Ev * mv0
initial(Iv1) <- init_Iv * mv0

initial(Sv2) <- init_Sv * mv0 * 0.01 ## CHECK - this is a bit of a fudge, but I think this should be fine as the lower carrying capacity will rapidly bring it down, but might want to start close to 0
initial(Ev2) <- init_Ev * mv0 * 0.01 ## CHECK - this is a bit of a fudge, but I think this should be fine as the lower carrying capacity will rapidly bring it down, but might want to start close to 0
initial(Iv2) <- init_Iv * mv0 * 0.01 ## CHECK - this is a bit of a fudge, but I think this should be fine as the lower carrying capacity will rapidly bring it down, but might want to start close to 0

# cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
# by age/het/int category, infectiousness depends on p_det which depends on detection immunity
cU <- user() # infectiousness U -> mosq
cD <- user() # infectiousness D -> mosq
cT <- user() # T -> mosq
gamma1 <- user() # fitted value of gamma1 characterises cA function
dim(cA) <- c(na,nh,num_int)
cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1

# Force of infection from humans to mosquitoes
dim(FOIvijk_species1) <- c(na,nh,num_int)
dim(FOIvijk_species2) <- c(na,nh,num_int)
omega <- user() #normalising constant for biting rates across diff age-groups

## CHECK
FOIvijk_species1[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq1[k]*foi_age[i]/omega
lag_FOIv_species1 <- sum(FOIvijk_species1)
FOIvijk_species2[1:na, 1:nh, 1:num_int] <- (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] * av_mosq2[k]*foi_age[i]/omega
lag_FOIv_species2 <- sum(FOIvijk_species2)

# Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
delayGam <- user() # Lag from parasites to infectious gametocytes
delayMos <- user() # Extrinsic incubation period.
FOIv_species1 <- delay(lag_FOIv_species1, delayGam)
FOIv_species2 <- delay(lag_FOIv_species2, delayGam)

# Number of mosquitoes that become infected at each time point
surv1 <- exp(-mu1*delayMos)
surv2 <- exp(-mu2*delayMos)

ince1 <- FOIv_species1 * Sv1 ## CHECK
lag_incv1 <- ince1 * surv1
incv1 <- delay(lag_incv1, delayMos)

ince2 <- FOIv_species2 * Sv2 ## CHECK
lag_incv2 <- ince2 * surv2
incv2 <- delay(lag_incv2, delayMos)

# Number of mosquitoes born (depends on PL, number of larvae), or is constant outside of seasonality
betaa1 <- 0.5*PL1/dPL
betaa2 <- 0.5*PL2/dPL

deriv(Sv1) <- -ince1 - mu1*Sv1 + betaa1
deriv(Ev1) <- ince1 - incv1 - mu1*Ev1
deriv(Iv1) <- incv1 - mu1*Iv1

deriv(Sv2) <- -ince2 - mu2*Sv2 + betaa2
deriv(Ev2) <- ince2 - incv2 - mu2*Ev2
deriv(Iv2) <- incv2 - mu2*Iv2

# Total mosquito population
Iv <- Iv1 + Iv2
output(Iv) <- TRUE
mv1 <-  Sv1+Ev1+Iv1
mv2 <- Sv2+Ev2+Iv2
mv <- mv1 + mv2
sp_rate_species1 <- Iv1/mv1
sp_rate_species2 <- Iv2/mv2
output(sporozoite_rate_species1) <- sp_rate_species1
output(sporozoite_rate_species2) <- sp_rate_species2

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
muPL <- user() #daily den. dep. mortality rate of pupae
muEL <- user() #daily den. dep. mortality rate of early stage
gammaL <- user() # eff. of den. dep. on late stage relative to early stage

# fitted entomological parameters:
mv0 <- user() # initial mosquito density
mu0 <- user() # baseline mosquito death rate
tau1 <- user() # duration of host-seeking behaviour
tau2 <- user() # duration of resting behaviour
p10 <- user() # prob of surviving 1 feeding cycle
p2 <- user() # prob of surviving one resting cycle
betaL <- user() # maximum number of eggs per oviposition per mosq

# parameters for species 2 increasing in abundance over time (NOTE this is distinct and on top of the custom seasonality)
density_vec[] <- user()
dim(density_vec) <- time_length
density_vec_sp1[] <- user()
dim(density_vec_sp1) <- time_length
time_length <- user()

# Entomological variables:
eov1 <- betaL/mu1*(exp(mu1/fv1)-1)
beta_larval1 <- eov1*mu1*exp(-mu1/fv1)/(1-exp(-mu1/fv1)) # Number of eggs laid per day

eov2 <- betaL/mu2*(exp(mu2/fv2)-1)
beta_larval2 <- eov2*mu2*exp(-mu2/fv2)/(1-exp(-mu2/fv2)) # Number of eggs laid per day

b_lambda <- (gammaL*muLL/muEL-dEL/dLL+(gammaL-1)*muLL*dEL) ## Parameters assumed to be the same for Species 1 and Species 2
lambda_species1 <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval1*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))
lambda_species2 <- -0.5*b_lambda + sqrt(0.25*b_lambda^2 + gammaL*beta_larval2*muLL*dEL/(2*muEL*mu0*dLL*(1+dPL*muPL)))

K0_species1 <- if(as.integer(t) == 0) 2*density_vec_sp1[as.integer(1)]*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda_species1+1)/(lambda_species1/(muLL*dEL)-1/(muLL*dLL)-1) else 2*density_vec_sp1[as.integer(t)]*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda_species1+1)/(lambda_species1/(muLL*dEL)-1/(muLL*dLL)-1)
## density_vec is defined above and is a monotonically increasing function designed to mimic and simulate invasion and establishment of An. stephensi
K0_species2 <- if(as.integer(t) == 0) 2*density_vec[as.integer(1)]*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda_species2+1)/(lambda_species2/(muLL*dEL)-1/(muLL*dLL)-1) else 2*density_vec[as.integer(t)]*dLL*mu0*(1+dPL*muPL)*gammaL*(lambda_species2+1)/(lambda_species2/(muLL*dEL)-1/(muLL*dLL)-1)

# Seasonal carrying capacity KL = base carrying capacity K0 * effect for time of year theta:
KL_species1 <- K0_species1 * theta_species1
KL_species2 <- K0_species2 * theta_species2

fv1 <- 1/( tau1/(1-zbar1) + tau2 ) # mosquito feeding rate for species 1 (zbar is a derived intervention param.)
fv2 <- 1/( tau1/(1-zbar2) + tau2 ) # mosquito feeding rate species 2 (zbar is a derived intervention param.)

mu1 <- -fv1*log(p1_1*p2) # mosquito death rate - species 1
mu2 <- -fv2*log(p1_2*p2) # mosquito death rate - species 2

# finding equilibrium and initial values for EL, LL & PL
init_PL <- user()
init_LL <- user()
init_EL <- user()

initial(PL1) <- init_PL
initial(LL1) <- init_LL
initial(EL1) <- init_EL

initial(PL2) <- init_PL * 0.01 ## CHECK - this is a bit of a fudge, but aim here is to initialise with a small, non-zero number (to reflect absence of stephensi at start of model run) so I think hopefully okay!
initial(LL2) <- init_LL * 0.01 ## CHECK - this is a bit of a fudge, but aim here is to initialise with a small, non-zero number (to reflect absence of stephensi at start of model run) so I think hopefully okay!
initial(EL2) <- init_EL * 0.01 ## CHECK - this is a bit of a fudge, but aim here is to initialise with a small, non-zero number (to reflect absence of stephensi at start of model run) so I think hopefully okay!

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL1) <- beta_larval1*mv1-muEL*(1+(EL1+LL1)/KL_species1)*EL1 - EL1/dEL # Note assumption of no overlap in breeding sites and no competition between Species 1 and 2

# egg hatching - den. dep. mortality - maturing larvae
deriv(LL1) <- EL1/dEL - muLL*(1+gammaL*(EL1+LL1)/KL_species1)*LL1 - LL1/dLL # Note assumption of no overlap in breeding sites and no competition between Species 1 and 2

# pupae - mortality - fully developed pupae
deriv(PL1) <- LL1/dLL - muPL*PL1 - PL1/dPL

# (beta_larval (egg rate) * total mosquito) - den. dep. egg mortality - egg hatching
deriv(EL2) <- beta_larval2*mv2-muEL*(1+(EL2+LL2)/KL_species2)*EL2 - EL2/dEL # Note assumption of no overlap in breeding sites and no competition between Species 1 and 2

# egg hatching - den. dep. mortality - maturing larvae
deriv(LL2) <- EL2/dEL - muLL*(1+gammaL*(EL2+LL2)/KL_species2)*LL2 - LL2/dLL # Note assumption of no overlap in breeding sites and no competition between Species 1 and 2

# pupae - mortality - fully developed pupae
deriv(PL2) <- LL2/dLL - muPL*PL2 - PL2/dPL

##------------------------------------------------------------------------------
########################
## INTERVENTION MODEL ## requires CHECK in big way
########################
##------------------------------------------------------------------------------

# See supplementary materials S2 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# general parameters
ITN_IRS_on <- user() # days after which interventions begin
num_int <- user() # number of intervention categories, ITN only, IRS only, neither, both
itn_cov <- user() # proportion of population covered by ITN
irs_cov <- user() # proportion of population covered by IRS

# cov is a vector of coverages for each intervention category:
dim(cov_) <- 4
cov_[1] <- (1-itn_cov)*(1-irs_cov)  # {No intervention}
cov_[2] <- itn_cov*(1-irs_cov) # 	   {ITN only}
cov_[3] <- (1-itn_cov)*irs_cov	#      {IRS only}
cov_[4] <- itn_cov*irs_cov #	   {Both ITN and IRS}
cov[] <- cov_[i]
dim(cov) <- num_int

IRS_interval <- user() # how long IRS lasts
ITN_interval <- user() # how long ITN lasts

## Bionomics Species 1
chi_species1 <- user() # proportion of vector endophily
Q0_species1 <- user() # proportion of anthropophagy
bites_Bed_species1 <- user() # endophagy in bed
bites_Indoors_species1 <- user() # endophagy indoors

## Bionomics Species 2
chi_species2 <- user() # proportion of vector endophily
Q0_species2 <- user() # proportion of anthropophagy
bites_Bed_species2 <- user() # endophagy in bed
bites_Indoors_species2 <- user() # endophagy indoors

# General intervention model terminology:
# r - probability of trying to repeat feed after hitting ITN/IRS
# d - probability of dying after hitting ITN/IRS
# s - probability of successful feed after hitting ITN/IRS

# The maximum (and then minimum) r and d values for ITN/IRS on day 0 before they decay
# r_ITN0 <- user()
# d_ITN0 <- user()
# d_IRS0 <- user()
# r_IRS0 <- user()
# r_ITN1 <- user() # min
# irs_loss <- user()
# itn_loss <- user()

r_ITN0_1 <- user() # max repellency by nets (Species 1)
d_ITN0_1 <- user() # max deaths after hitting net (Species 1)
r_ITN1 <- user() # non-zero level of repellency (Species 1)

r_ITN0_2 <- user() # max repellency by nets (Species 2)
d_ITN0_2 <- user() # max deaths after hitting net (Species 2)
r_ITN2 <- user() # non-zero level of repellency (Species 2)

r_IRS0_1 <- user() # max repellency by IRS
d_IRS0_1 <- user() # max death after IRS

r_IRS0_2 <- user() # max repellency by IRS
d_IRS0_2 <- user() # max death after IRS

irs_loss <- user()
itn_loss <- user()

# Calculates decay for ITN/IRS
ITN_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%ITN_interval) * itn_loss)
IRS_decay = if(t < ITN_IRS_on) 0 else exp(-((t-ITN_IRS_on)%%IRS_interval) * irs_loss)

# The r,d and s values turn on after ITN_IRS_on and decay accordingly
d_ITN_1 <- if(t < ITN_IRS_on) 0 else d_ITN0_1*ITN_decay
r_ITN_1 <- if(t < ITN_IRS_on) 0 else r_ITN1 + (r_ITN0_1 - r_ITN1)*ITN_decay
s_ITN_1 <- if(t < ITN_IRS_on) 1 else 1 - d_ITN_1 - r_ITN_1

d_ITN_2 <- if(t < ITN_IRS_on) 0 else d_ITN0_2*ITN_decay
r_ITN_2 <- if(t < ITN_IRS_on) 0 else r_ITN2 + (r_ITN0_2 - r_ITN2)*ITN_decay
s_ITN_2 <- if(t < ITN_IRS_on) 1 else 1 - d_ITN_2 - r_ITN_2

r_IRS_1 <- if(t < ITN_IRS_on) 0 else r_IRS0_1*IRS_decay
d_IRS_1 <- if(t < ITN_IRS_on) 0 else chi_species1*d_IRS0_1*IRS_decay
s_IRS_1 <- if(t < ITN_IRS_on) 1 else 1 - d_IRS_1

r_IRS_2 <- if(t < ITN_IRS_on) 0 else r_IRS0_2*IRS_decay
d_IRS_2 <- if(t < ITN_IRS_on) 0 else chi_species2*d_IRS0_2*IRS_decay
s_IRS_2 <- if(t < ITN_IRS_on) 1 else 1 - d_IRS_2

# probability that mosquito bites and survives for each intervention category
dim(w1_) <- 4
w1_[1] <- 1
w1_[2] <- 1 - bites_Bed_species1 + bites_Bed_species1*s_ITN_1
w1_[3] <- 1 - bites_Indoors_species1 + bites_Indoors_species1*(1-r_IRS_1)*s_IRS_1
w1_[4] <- 1 - bites_Indoors_species1 + bites_Bed_species1*(1-r_IRS_1)*s_ITN_1*s_IRS_1 + (bites_Indoors_species1 - bites_Bed_species1)*(1-r_IRS_1)*s_IRS_1
w1[] <- w1_[i]
dim(w1) <- num_int

dim(w2_) <- 4
w2_[1] <- 1
w2_[2] <- 1 - bites_Bed_species2 + bites_Bed_species2*s_ITN_2
w2_[3] <- 1 - bites_Indoors_species2 + bites_Indoors_species2*(1-r_IRS_2)*s_IRS_2
w2_[4] <- 1 - bites_Indoors_species2 + bites_Bed_species2*(1-r_IRS_2)*s_ITN_2*s_IRS_2 + (bites_Indoors_species2 - bites_Bed_species2)*(1-r_IRS_2)*s_IRS_2
w2[] <- w2_[i]
dim(w2) <- num_int

# probability that mosq feeds during a single attempt for each int. cat.
dim(yy1_) <- 4
yy1_[1] <- 1
yy1_[2] <- w1_[2]
yy1_[3] <- 1 - bites_Indoors_species1 + bites_Indoors_species1*(1-r_IRS_1)
yy1_[4] <- 1 - bites_Indoors_species1 + bites_Bed_species1*(1-r_IRS_1)*s_ITN_1 + (bites_Indoors_species1 - bites_Bed_species1)*(1-r_IRS_1)
yy1[] <- yy1_[i]
dim(yy1) <- num_int

dim(yy2_) <- 4
yy2_[1] <- 1
yy2_[2] <- w2_[2]
yy2_[3] <- 1 - bites_Indoors_species2 + bites_Indoors_species2*(1-r_IRS_2)
yy2_[4] <- 1 - bites_Indoors_species2 + bites_Bed_species2*(1-r_IRS_2)*s_ITN_2 + (bites_Indoors_species2 - bites_Bed_species2)*(1-r_IRS_2)
yy2[] <- yy2_[i]
dim(yy2) <- num_int

# probability that mosquito is repelled during a single attempt for each int. cat.
dim(z1_) <- 4
z1_[1] <- 0
z1_[2] <- bites_Bed_species1*r_ITN_1
z1_[3] <- bites_Indoors_species1*r_IRS_1
z1_[4] <- bites_Bed_species1*(r_IRS_1 + (1-r_IRS_1)*r_ITN_1) + (bites_Indoors_species1 - bites_Bed_species1)*r_IRS_1
z1[] <- z1_[i]
dim(z1) <- num_int

dim(z2_) <- 4
z2_[1] <- 0
z2_[2] <- bites_Bed_species2*r_ITN_2
z2_[3] <- bites_Indoors_species2*r_IRS_2
z2_[4] <- bites_Bed_species2*(r_IRS_2+ (1-r_IRS_2)*r_ITN_2) + (bites_Indoors_species2 - bites_Bed_species2)*r_IRS_2
z2[] <- z2_[i]
dim(z2) <- num_int

# Calculating Z (zbar) and W (wbar) - see Supplementary materials 2 for details
dim(zhi1) <- num_int
dim(whi1) <- num_int
zhi1[1:num_int] <- cov[i]*z1[i]
whi1[1:num_int] <- cov[i]*w1[i]
zh1 <- if(t < ITN_IRS_on) 0 else sum(zhi1)
wh1 <- if(t < ITN_IRS_on) 1 else sum(whi1)

dim(zhi2) <- num_int
dim(whi2) <- num_int
zhi2[1:num_int] <- cov[i]*z2[i]
whi2[1:num_int] <- cov[i]*w2[i]
zh2 <- if(t < ITN_IRS_on) 0 else sum(zhi2)
wh2 <- if(t < ITN_IRS_on) 1 else sum(whi2)

# Z (zbar) - average probability of mosquito trying again during single feeding attempt
zbar1 <- Q0_species1*zh1
zbar2 <- Q0_species2*zh2

# W (wbar) - average probability of mosquito successfully feeding during single attempt
wbar1 <- 1 - Q0_species1 + Q0_species1*wh1
wbar2 <- 1 - Q0_species2 + Q0_species2*wh2

# p1 is the updated p10 given that interventions are now in place:
p1_1 <- wbar1*p10/(1-zbar1*p10)
p1_2 <- wbar2*p10/(1-zbar2*p10)

Q_1 <- 1-(1-Q0_species1)/wbar1 # updated anthropophagy given interventions
Q_2 <- 1-(1-Q0_species2)/wbar2 # updated anthropophagy given interventions

av1 <- fv1*Q_1 # biting rate on humans by species 1
av2 <- fv2*Q_2 # biting rate on humans by species 2

## Species-specific FOIs (calculated above, but which require species-specific av_mosq and av_human) - we then add the FOI on human ones up later when calculating EIR
dim(av_mosq1) <- num_int
av_mosq1[1:num_int] <- av1*w1[i]/wh1 # rate at which mosquitoes of species 1 bite each int. cat. - a bit unclear here why we multiply by w1/wh1 - this isn't consistent with Jamie's 2010 paper
dim(av_human1) <- num_int
av_human1[1:num_int] <- av1*yy1[i]/wh1 # biting rate by species 1 on humans of in each int. cat. - this is consistent with Jamie's 2010 paper in Protocol S2

dim(av_mosq2) <- num_int
av_mosq2[1:num_int] <- av2*w2[i]/wh2 # rate at which mosquitoes of species 2 bite each int. cat. - a bit unclear here why we multiply by w1/wh1 - this isn't consistent with Jamie's 2010 paper
dim(av_human2) <- num_int
av_human2[1:num_int] <- av2*yy2[i]/wh2 # biting rate by species 2 on humans in each int. cat. - this is consistent with Jamie's 2010 paper in Protocol S2

##------------------------------------------------------------------------------
###################
## MODEL OUTPUTS ##
###################
##------------------------------------------------------------------------------

# Outputs for each compartment across the sum across all ages, biting heterogeneities and intervention categories
output(Sout) <- sum(S[,,])
output(Tout) <- sum(T[,,])
output(Dout) <- sum(D[,,])
output(Aout) <- sum(A[,,])
output(Uout) <- sum(U[,,])
output(Pout) <- sum(P[,,])

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
output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])

dim(prev0to5) <- c(age05,nh,num_int)
prev0to5[1:age05,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]
output(prev05) <- sum(prev0to5[,,])/sum(den[1:age05])

# slide positivity in 0 -5 year age bracket
dim(clin_inc0to5) <- c(age05,nh,num_int)
clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
output(inc) <- sum(clin_inc[,,])

# Param checking outputs
output(mu1) <- mu1
output(mu2) <- mu2

output(beta_larval1) <- beta_larval1
output(beta_larval2) <- beta_larval2

output(mv) <- mv
output(mv1) <- mv1
output(mv2) <- mv2

output(wh1) <- wh1
output(wh2) <- wh2

output(d_ITN_1) <- d_ITN_1
output(r_ITN_1) <- r_ITN_1
output(s_ITN_1) <- s_ITN_1

output(d_ITN_2) <- d_ITN_2
output(r_ITN_2) <- r_ITN_2
output(s_ITN_2) <- s_ITN_2

output(d_IRS_1) <- d_IRS_1
output(r_IRS_1) <- r_IRS_1
output(s_IRS_1) <- s_IRS_1

output(d_IRS_2) <- d_IRS_2
output(r_IRS_2) <- r_IRS_2
output(s_IRS_2) <- s_IRS_2

output(p1_1) <- p1_1 # updated prob of surviving 1 feeding cycle given int of pop1
output(p1_2) <- p1_2 # updated prob of surviving 1 feeding cycle given int of pop2

output(wbar1) <- wbar1
output(wbar2) <- wbar2

output(zbar1) <- zbar1
output(zbar2) <- zbar2

output(Q_1) <- Q_1 # updated anthropophagy given int of pop1
output(Q_2) <- Q_2 # updated anthropophagy given int of pop2

output(av1) <- av1 # biting rate of pop1
output(av2) <- av2 # biting rate of pop2
#output(av) <- av # weighted biting rate

output(cov[]) <- TRUE
output(K0_species1) <- K0_species1
output(K0_species2) <- K0_species2
output(KL_species1) <- KL_species1
output(KL_species2) <- KL_species2
output(FOIv_species1) <- FOIv_species1
output(FOIv_species2) <- FOIv_species2
