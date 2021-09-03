#######################################################################
### SEPM non-local exact D-optimal design by excursion : use case 2 ###
#######################################################################

library(here)

source(here('UseCase2/SEPM_OptDesign_uc2_l.R'))

## min and max \Theta
## as in https://doi.org/10.1111/emr.12350

parmin=c("b"=0.16,"D"=2)
parmax=c("b"=0.27,"D"=30)

#############################################
## sample \Theta via a space filing design ##
#############################################

require(lhs)
# \in [0,1] 

## Load a predefine space-filling design as an equivalent to a predefine seed 

load(here("inputs/others/nlhs_uc2.RData"))

## uncomment to generate a new space-filling design instead 

# nlhs=randomLHS(10,2)
# 
# for (i in 1:length(parmin))
# {
#   nlhs[,i]=parmin[i]+(nlhs[,i]*(parmax[i]-parmin[i]))
# }

####### This is a unitary test, with only a few samples for the space-filling design
####### Computation takes around 30' on a 2018's desktop. I memoized the jacobian to speed up the D-optimal search


startC=proc.time()
require(doParallel)

# cluster definition
ncor=3
clust <- makeCluster(ncor)
registerDoParallel(clust, cores = ncor)

dataDesign = foreach(i = 1:dim(nlhs)[1]) %dopar% {

  try({

    # avoid ram usage concurrency at the begining (maybe cause of dead running processes)
    Sys.sleep((i-1) * 10)

    require(maptools)
    require(rgeos)
    require(raster)
    require(rgdal)
    require(spdep)
    require(spatstat)
    library(Matrix)
    library(gstat)
    library(geospt)
    library(numDeriv)
    library(deSolve)
    library(ReacTran)

    partmp=unlist(params)[!names(params) == "mu" ]
    ntmp=names(partmp) 
    partmp=as.vector(nlhs[i,])
    names(partmp)=ntmp
      resJac=jacobian(ForFxGrad,x=partmp,where=1:32)
      tr.uFIM=cbind(1:32,resJac)
      tr.uFIM
  })
}

stopCluster(clust)

save(dataDesign,nlhs,file=here("data/lhs_uc2_nl.RData"))
stopC=proc.time()
print(stopC-startC)

#### reloading the memoized computations

# load("/mnt/stockage/Boulot/Codes/DoE_next/LANDOE/data/lhs_uc2_nl.RData")
# load("data/lhs_2.RData")

#########################################
# We solve the non-local exact design ###
#########################################

# 8 out of 32 transects
nb=8

repetitions=10
qss=list()
for (repetition in 1:repetitions)
{
e=Sys.time()
# set.seed(75)
idxT=1:32
idxS=sample(idxT,nb)
idxSO=idxS
idxC=idxT[!(idxT %in% idxS)]

# first det value 
detcur=0
for (m in dataDesign)
{
  detcur=detcur+det(map2FIM(map=m,candidateset=idxS))
}

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
      detemp=0
      for (m in dataDesign)
      {
        detemp=detemp+det(map2FIM(map=m,candidateset=idxTmp))
      }
      detStep[[paste0(i,"_",j)]] = detemp
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

##############################################################
## Creation of several simple plots to assess the results  ###
##############################################################

######## Equivalent to figure 5b

# half but non-localy D optimal 
plot(xini,legend=F,main='Non-local D-optimal')
plot(transects[ transects$id %in% idxS,],col='red',legend=F,add=T)

dev.copy2pdf(file=here("data/SEPM_uc2_fig5b.pdf"))
dev.off()

####################
##### ellipses #####
####################

######## Equivalent to figure 5c

idxT=1:32
# fixed
idxSO=c(28, 31, 32, 10, 22,  9,  1, 17)
# or random
# idxSO=sample(idxT,nb)

par(mfrow=c(1,1))

require(ellipse)

tr.uFIM=cbind(1:32,resJacLoc)
ori=ellipse(solve(map2FIM(map=tr.uFIM,candidateset=idxSO))*(sig_uc2)**2,which=1:2,npoints=1000)
opt=ellipse(solve(map2FIM(map=tr.uFIM,candidateset=idxS))*(sig_uc2)**2,which=1:2,npoints=1000)
gro=rbind(ori,opt)
plot(ori,type='l',col='blue',ylim=c(min(gro[,2]),max(gro[,2])),xlim=c(min(gro[,1]),max(gro[,1])),main="b & D") #,xlab=na[de],ylab=na[fi])
lines(opt,col='red')

dev.copy2pdf(file=here("data/SEPM_uc2_fig5c.pdf"))
dev.off()

