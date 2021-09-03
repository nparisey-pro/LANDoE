#########################################################################
### SEPM non local exact D-optimal design by excursion for use case 1 ###
#########################################################################

library(here)

source(here('UseCase1/SEPM_OptDesign_uc1_l.R'))

# from mean value, let's set mins and maxs
# half to twice as many descendants 
parmin=parmax=partmp=unlist(parmsPDE)
parmax[1:3]=partmp[1:3]+log(2)/Tmax
parmin[1:3]=partmp[1:3]-log(2)/Tmax
 
# D, beta and mu are *not* to be treated the same
# D : 2 times faster to 2 times slower 
parmin['D']=partmp['D']*(1/2)^2
parmax['D']=partmp['D']*2^2
# mu : live half as long to twice as long 
parmax['mu']=partmp['mu']*(1/0.5)
parmin['mu']=partmp['mu']*(1/2)
# beta : 
parmax['beta']=partmp['beta']*1.1
parmin['beta']=partmp['beta']*0.9

parmin=parmin[!names(parmin)=="alpha"]
parmax=parmax[!names(parmax)=="alpha"]


####### This is a unitary test, with only a few samples for the space-filling design
####### Computation takes around 10' on a 2018's desktop. I memoized the jacobian to speed up the D-optimal search

## Load a predefine space-filling design as an equivalent to a predefine seed 

load(here("inputs/others/nlhs_uc1.RData"))

## uncomment to generate a new space-filling design instead 

# nbS=10
# require(lhs)
# # \in [0,1] 
# nlhs=randomLHS(nbS,6)
# 
# for (i in 1:length(parmin))
# {
#   nlhs[,i]=parmin[i]+(nlhs[,i]*(parmax[i]-parmin[i]))
# }


####### 

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

    partmp=unlist(parmsPDE)[!names(parmsPDE) == "alpha" ]
    ntmp=names(partmp) 
    partmp=as.vector(nlhs[i,])
    names(partmp)=ntmp
    resJac=jacobian(ForFxGrad,x=partmp,where=spl.pts.map)
    partmp=as.list(partmp)
    partmp$alpha=parmsPDE[["alpha"]]
    vv=obs_eq_CC2D(params=partmp,ix.s = spl.pts.map)

    spl.pts.map.nlFIM=cbind(cbind(rep(spl.pts.map[,1],9),rep(spl.pts.map[,2],9)),
                            resJac,
                            unname(vv))

    spl.pts.map.nlFIM
  })
}

stopCluster(clust)
# serialization (see below)
save(dataDesign,nlhs,file=here("data/lhs_uc1_nl.RData"))
stopC=proc.time()
print(stopC-startC)

#### if you want to re-used the serialized objets
# rm(nlhs,dataDesign)
# load(here("data/lhs_uc1_nl.RData"))

##################################################
# Now, we can solve the non-local exact design ###
##################################################

# Computation takes around 1' on a 2018's desktop.

# let's start from the local optimal design (i.e. idxS) 
# save positions of local optimum
idxSl=idxS

repetitions=10
qss=list()
for (repetition in 1:repetitions)
{
idxT=1:dim(spl.pts.map)[1]
idxS=sample(idxT,6)

idxC=idxT[!(idxT %in% idxS)]

# first det value
idxFIM=rep(F,dim(spl.pts.map)[1]) ; idxFIM[idxS]=T ; idxFIM=rep(idxFIM,9)

FIMs=list()
l=1
for (spl.pts.map.nlFIM in dataDesign)
{
  FIM=matrix(0,6,6)
  for (i in which(idxFIM))
  {
    FIM=FIM+as.vector(spl.pts.map.nlFIM[i,3:8]) %*% t(as.vector(spl.pts.map.nlFIM[i,3:8])) * spl.pts.map.nlFIM[i,9]^{-1}
  }
  FIMs[[l]]=FIM
  l=l+1
}
detcur=0
for (fi in FIMs)
{
  detcur=detcur + det(fi) # detcur=det(FIM)  
}

# # detcur=det(spl.pts.map2FIM(candidateset=idxS))
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
       idxFIM=rep(F,dim(spl.pts.map)[1]) ; idxFIM[idxTmp]=T ; idxFIM=rep(idxFIM,9)
       FIMs=list()
       l=1
       for (spl.pts.map.nlFIM in dataDesign)
       {
       FIM=matrix(0,6,6)
       for (k in which(idxFIM))
       {
         FIM=FIM+as.vector(spl.pts.map.nlFIM[k,3:8]) %*% t(as.vector(spl.pts.map.nlFIM[k,3:8])) * spl.pts.map.nlFIM[k,9]^{-1}
       }
       FIMs[[l]]=FIM
       l=l+1
       }
       # FIM=spl.pts.map2FIM(candidateset=idxTmp)
       detemp=0
       for (fi in FIMs)
       {detemp=detemp+det(fi)}
       detStep[[paste0(i,"_",j)]] = detemp #det(FIM)
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
# f=Sys.time()
# print(f-e)
q=idxS
qss[[repetition]]=list(q,detcur)
}

qIdx=which.max(unlist(qss)[seq(from=7,to=length(unlist(qss)),by=7)])

idxS=qss[[qIdx]][[1]]

##############################################################
## Creation of several simple plots to assess the results  ###
##############################################################

######## Equivalent to figure 2a 

par(mfrow=c(1,4))

img=load.image(here('inputs/gis/figure_OptimalDesign.tiff'))
plot(0,type='n', main="Non-local D-optimal design", xlab="long", ylab="lat",xlim=c(0,1),ylim=c(0,1))
lim <- par()
rasterImage(img, 0, 0, 1, 1)
grid()
# non-local
PP=spl.pts.map[idxS,]
### visual correction to take into account the coarcer map used for computation compared to the display map
PP[round(PP[,'x'],2) == 0.45,'x']=0.45-0.015 # -0.015

points(PP[,2],(1-PP[,1]),pch='*',col='black',cex=1.7)

FIMOpt=spl.pts.map2FIM(candidateset=idxS)

xx=combn(3,2)
for (k in 1:dim(xx)[2])
{
  i=xx[1,k]
  j=xx[2,k]
  ori=ellipse(solve(FIMCLASS),which=c(i,j),npoints=1000)
  opt=ellipse(solve(FIMOpt),which=c(i,j),npoints=1000)
  gro=rbind(ori,opt)
  plot(ori,type='l',col='blue',ylim=c(min(gro[,2]),max(gro[,2])),xlim=c(min(gro[,1]),max(gro[,1])),main=paste0("growth rates ",i,"&",j)) #,xlab=na[de],ylab=na[fi])
  lines(opt,col='red')
}


dev.copy2pdf(file=here("data/SEPM_uc1_fig2a.pdf"))
dev.off()
