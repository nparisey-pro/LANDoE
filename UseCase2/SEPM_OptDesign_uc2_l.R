####################################################################
### SEPM local exact D-optimal design by excursion : use case 2  ###
####################################################################

library(here)

source(here('UseCase2/SEPM_uc2.R'))

require(numDeriv)

ForFxGrad=function(x,p=params,where=1:32,timezone=Tfin) 
{
  x=as.list(x)
  x$mu=p$mu
  return( obs_eq_horses2D(params=x,ix.s = where,ix.t=timezone) )
}

##############################################
## On top of the PDE domain, we want an     ##
## ensemble of potential sampling transects ##
##############################################

# 8 out of 32
nb=8

##############################################
# We can compute the model's Jacobians #######
# for a given vector of parameter      #######
##############################################
# This takes a few minutes (on my computer)

a=proc.time()
x_nominal=unlist(params)[!names(params) == "mu" ]
resJac=jacobian(ForFxGrad,x=x_nominal,where=1:32) 
resJacLoc=resJac
b=proc.time()
print(b-a)

### now we can enhance the sampling map 

tr.uFIM=cbind(1:32,resJac)

####################################################################
# for a given set of points on the sampling map, compute the FIM ###
####################################################################

map2FIM <- function(map=tr.uFIM,candidateset=idxS)
{
  # computing FIM from sampled points 
  FIM=matrix(0,2,2)
  for (i in candidateset)
  {
    FIM=FIM+as.vector(map[i,2:3]) %*% t(as.vector(map[i,2:3]))
  }
  return(FIM)
}

map2FIM(candidateset=sample(1:32,nb))


############################################################
# We solve the exact design on the sampling map using an ###  
# excursion algorithm e.g. an overly simplified Fedorov  ###
############################################################

# 8 out of 32
nb=8

repetitions=10
qss=list()
for (repetition in 1:repetitions)
{
  
e=Sys.time()
# set.seed(55)
idxT=1:32
idxS=sample(idxT,nb)
idxSO=idxS
idxC=idxT[!(idxT %in% idxS)]

# first det value
detcur=det(map2FIM(candidateset=idxS))

exch=T

while (exch)
{
  exch=F
  detStep=list()
  for (i in 1:length(idxS))
  {
    for (j in 1:length(idxC))
    {
      idxTmp=idxS
      idxTmp[i]=idxC[j]
      # computing FIM from current sampled points 
      FIM=map2FIM(candidateset=idxTmp)
      detStep[[paste0(i,"_",j)]] = det(FIM)
    }
  }
  if ( (max(unlist(detStep))/detcur) > 1.01 )
  {
    q=which.max(unlist(detStep))
    q=as.integer(unlist(strsplit(names(q),"_")))
    idxTmp=idxS
    idxS[q[1]]=idxC[q[2]]
    idxC[q[2]]=idxTmp[q[1]]
    detcur=max(unlist(detStep))
    exch=T 
  }
}
f=Sys.time()
print(f-e)
qss[[repetition]]=list(idxS,det(map2FIM(candidateset=idxS)))
}

indicOpt=which.max(unlist(qss)[seq((nb+1),length(unlist(qss)),(nb+1))])

idxS=qss[[indicOpt]][[1]]

# let's compare (random vs localy D optimal)
print(det(map2FIM(candidateset=idxSO)))
print(det(map2FIM(candidateset=idxS)))

##############################################################
## Creation of several simple plots to assess the results  ###
##############################################################

######## Equivalent to figure 4

par(mfrow=c(2,2))
# initial conditions
plot(xini,legend=F,main="Initial Densities")
# all transects
plot(xini,legend=F,main="Original Transects")
plot(transects,col='black',legend=F,add=T)
# half but random
plot(xini,legend=F,main="Random sampling of Transects")
plot(transects[ transects$id %in% idxSO,],col='blue',legend=F,add=T)
# half but localy D optimal 
plot(xini,legend=F,main="Localy D-optimal")
plot(transects[ transects$id %in% idxS,],col='red',legend=F,add=T)

dev.copy2pdf(file=here("data/SEPM_uc2_fig4.pdf"))
dev.off()

