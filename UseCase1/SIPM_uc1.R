##################################################################
### Spatialy Implicit Population Model (SIPM) : use case 1   #####
##################################################################

#######################
### SIPM parameters ###
#######################

## domain related parameters
#  0 < t < Tmax
Tmax=63

# sampling dates
spl.d=seq(0,56,7)

# number of traps
n_x=12

# landscape composition

# for an example subdomain (from doi:10.1016/j.ecocom.2016.07.003)
promega1= c('pg'= 0.37,'pc'=0.47,'ps'=0.07)

# simulated 'beta' diversity i.e. sampling scheme in composition 
# with 4 landcovers (rs,rc,rg and 0 (e.g. no growth rate in urban areas))

beta_design=matrix(runif(n_x*4),n_x,4)
tot_p=rowSums(beta_design)
for (i in 1:n_x){ beta_design[i,]=beta_design[i,]/tot_p[i] }
beta_design=beta_design[,-4]
colnames(beta_design)=c('pg','pc','ps')

# consistency check 
solve(t(beta_design) %*% beta_design)

## population dynamics parameters

P0=2*10^4
parms=list(rs=0.155,rg=0.304,rc=0.385,beta=0.123/6.3,mu=0.210)

# parameters boundaries

# 10 times more or less offsprings per capita (related to landcover)
parMin=parMax=parms[c('rs','rg','rc')]
tenTM=log(10)/Tmax
parMax[c('rs','rg','rc')]=max(unlist(parMax)+tenTM)
parMin[c('rs','rg','rc')]=min(unlist(parMin)-tenTM)

# livespan half as long to twice as long 
parMax$mu=parms$mu*2
parMin$mu=parms$mu/2

# +/- 10% in birth decay 
parMax$beta=parms$beta*1.1
parMin$beta=parms$beta*0.9

######################
### SIPM functions ###
######################

# SIPM function e.g. P(t)
P<-function(t,par,prop,iniP)
{
  pg=prop['pg']
  pc=prop['pc']
  ps=prop['ps']
  with(par, { return(P0 * exp ( (1/beta) * (1 - exp(-beta*t)) * (ps*rs + pg*rg + pc*rc ) - mu*t )) })  
}

# SIPM with measurement equation

SIPM_Obs_eq<-function(sampling.times=seq(0,56,7),params=parms,propr=promega1,tau=5,alpha=10^-4,iniP=P0)
{
  pred=c()
  for (da in sampling.times)
  {
    pred=c(pred,integrate(P,da,da+5,par=params,prop=propr,iniP=iniP)$value)
  }
  return(alpha*pred)
}

#### SIPM with measurement equation and poisson data model 

SIPM_Samples<-function(sampling.times=seq(0,56,7),params=parms,propr=promega1,tau=5,alpha=10^-4,iniP=P0)
{
  mean_sipm_pred=SIPM_Obs_eq(sampling.times,params,propr,tau,alpha,iniP)
  spl_sipm_pred=rpois(length(mean_sipm_pred),mean_sipm_pred)
  return(spl_sipm_pred)
}

# SIPM's observations at the landscape scale
# that is to say 'n' dynamics in 'n' compositional contexts

SIPM_Obs_Land<-function(sampling.times=seq(0,56,7),params=parms,tau=5,alpha=10^-4,iniP=P0,X=beta_design)
{
  pred_land=c()
  for (i in 1:dim(X)[1])
  {
    propr= X[i,]
    pred_land=c(pred_land,SIPM_Obs_eq(sampling.times,params,propr,tau,alpha,iniP)  )
  }
  return(pred_land)
}

