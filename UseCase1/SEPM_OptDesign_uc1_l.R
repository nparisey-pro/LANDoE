#####################################################################
### SEPM local exact D-optimal design by excursion for use case 1 ###
#####################################################################

library(here)

source(here('UseCase1/SEPM_uc1.R'))

require(numDeriv)

# to compute observation points
ForFxGrad=function(x,p=parmsPDE,where=spl.pts) 
{
  x=as.list(x)
  x$alpha=p$alpha
  return( obs_eq_CC2D(params=x,ix.s = where) )
}

###################################################
## On top of the PDE domain, we want an ensemble ##
## of potential sampling points (a sampling map) ##
###################################################

# resolution of the sampling map
step=10

# construction of the sampling map 
r <- raster(vals=step^2,extent(0,1,0,1),resolution=1/step)
spl.pts.map=coordinates(r)

### we can filter out certain landcover(s) and/or field(s) 
# e.g. here we don't want to sample on the roads (landcover '3')
rtmp=raster(as.matrix(land))
extract(rtmp,SpatialPoints(spl.pts.map)) != 3

tmp=spl.pts.map[,1]
spl.pts.map[,1]=spl.pts.map[,2]
spl.pts.map[,2]=(1-tmp)

spl.pts.map=spl.pts.map[extract(rtmp,SpatialPoints(spl.pts.map)) != 3,]

tmp1=spl.pts.map[,1]
tmp2=spl.pts.map[,2]

spl.pts.map[,2]=tmp1
spl.pts.map[,1]=1-tmp2

###################################################
# We can compute the model's Jacobians and variances 
# for a given vector of parameter (useful for a local design)
################################################### 
# This takes a few minutes (less than 3' on my computer)

x_nominal=unlist(parmsPDE)[!names(parmsPDE) == "alpha" ]
resJac=jacobian(ForFxGrad,x=x_nominal,where=spl.pts.map)
vv=obs_eq_CC2D(params=parmsPDE,ix.s = spl.pts.map)

### now we can enhance the sampling map 

spl.pts.map.uFIM=cbind(cbind(rep(spl.pts.map[,1],9),rep(spl.pts.map[,2],9)),
                        resJac,
                        unname(vv))

####################################################################
# for a given set of points on the sampling map, compute the FIM ###
####################################################################

spl.pts.map2FIM <- function(map=spl.pts.map.uFIM,candidateset=idxS,idxparms=3:8,idxvar=9,repTim=9)
{
  lg=dim(unique(map[,1:2]))[1]
  # computing FIM from sampled points 
  idxFIM=rep(F,lg) ; idxFIM[candidateset]=T ; idxFIM=rep(idxFIM,repTim)
  FIM=matrix(0,length(idxparms),length(idxparms))
  for (i in which(idxFIM))
  {
    FIM=FIM+as.vector(map[i,idxparms]) %*% t(as.vector(map[i,idxparms])) * map[i,idxvar]^{-1}
  }
return(FIM)
}

############################################################
# We solve the exact design on the sampling map using an ###  
# excursion algorithm e.g. an overly simplified Fedorov  ###
############################################################

# This should take ...
# 40''(on my computer) using spl.pts.map2FIM(xx)
# 6'' (on my computer) with 'copy pasta' code
# My best guess : the function's call are killing the perf with overheads
# Copy pasta it is then... 

repetitions=10
ListDoptLoc=list()
for (repetition in 1:repetitions)
{
e=Sys.time()
idxT=1:dim(spl.pts.map)[1]
idxS=sample(idxT,6)
idxSO=idxS
idxC=idxT[!(idxT %in% idxS)]

# first det value
idxFIM=rep(F,dim(spl.pts.map)[1]) ; idxFIM[idxS]=T ; idxFIM=rep(idxFIM,9)
FIM=matrix(0,6,6)
for (i in which(idxFIM))
{
  FIM=FIM+as.vector(spl.pts.map.uFIM[i,3:8]) %*% t(as.vector(spl.pts.map.uFIM[i,3:8])) * spl.pts.map.uFIM[i,9]^{-1}
}
detcur=det(FIM)

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
      FIM=matrix(0,6,6)
      for (k in which(idxFIM))
      {
        FIM=FIM+as.vector(spl.pts.map.uFIM[k,3:8]) %*% t(as.vector(spl.pts.map.uFIM[k,3:8])) * spl.pts.map.uFIM[k,9]^{-1}
      }
      
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

ListDoptLoc[[repetition]]=idxS
}

oussa=0
detCourant=0
for (k in 1:length(ListDoptLoc))
{
  detHere=det(spl.pts.map2FIM(candidateset=ListDoptLoc[[k]]))
  if (detHere > detCourant) {detCourant=detHere ; oussa=k}
}

idxS=ListDoptLoc[[k]]

## chosed 'random' for graphical purpose
idxSO=c(70,74,23,75,86,36)
##

##############################################################
## Creation of several simple plots to assess the results  ###
##############################################################

######## Equivalent to figure 1 

# 4 panels
# initial conditions
# clustered design
# randomly positioned traps
# D-optimal 

require(imager)
img=load.image(here('inputs/gis/figure_OptimalDesign.tiff'))
par(mfrow=c(2,2))

# initial conditions
plot(uini/1000,main="Initial conditions") 

cols=c("black","black","blue")
pchs=c('x','+','*')
cexs=c(1.2,1.2,1.7)
titles=c("Clustered design","Random design","D-optimal design")
i=1
for (PP in list(spl.pts,spl.pts.map[idxSO,],spl.pts.map[idxS,]))
{
  plot(0,type='n', main=titles[i], xlab="long", ylab="lat",xlim=c(0,1),ylim=c(0,1))
  lim <- par()
  rasterImage(img, 0, 0, 1, 1)
  grid()
  points(PP[,2],(1-PP[,1]),pch=pchs[i],col=cols[i],cex=cexs[i])
  i=i+1
}

dev.copy2pdf(file=here("data/SEPM_uc1_fig1.pdf"))
dev.off()


######## Equivalent to figure 3.a

#################################################################
## Computing D-efficiency and Traveling Salesman Tour Length  ###
#################################################################

x_nominal=unlist(parmsPDE)[!names(parmsPDE) == "alpha" ]
resJacCLASS=jacobian(ForFxGrad,x=x_nominal,where=spl.pts)
vvCLASS=obs_eq_CC2D(params=parmsPDE,ix.s = spl.pts)

FIMCLASS=matrix(0,6,6)
for (i in 1:dim(resJacCLASS)[1])
{
  FIMCLASS=FIMCLASS+as.vector(resJacCLASS[i,]) %*% t(as.vector(resJacCLASS[i,])) * vvCLASS[i]^{-1}
}
detCLASS=det(FIMCLASS)

require(TSP)

set.seed(42) 

# clustered design tour length and D crit
etsp <- ETSP(data.frame(x = spl.pts[,1], y = spl.pts[,2]))
tour <- solve_TSP(etsp)

long=tour_length(tour)
dets=c(detCLASS)

# optimal design tour length and D crit
etsp <- ETSP(data.frame(x = spl.pts.map[idxS,1], y = spl.pts.map[idxS,2]))
tour <- solve_TSP(etsp)
long=c(long,tour_length(tour))

FIMOpt=spl.pts.map2FIM(candidateset=idxS)
dets=c(dets,det(FIMOpt))

##
llPP=list()
##

for (j in 1:100)
{
  idxRd=sample(1:dim(spl.pts.map)[1],6)
  llPP[[j]]=idxRd
  # Tour Length 
  etsp <- ETSP(data.frame(x = spl.pts.map[idxRd,1], y = spl.pts.map[idxRd,2]))
  tour <- solve_TSP(etsp)
  long=c(long,tour_length(tour))
  # D criterion
  FIMRd=spl.pts.map2FIM(candidateset=idxRd)
  dets=c(dets,det(FIMRd))
}

require(rPref)

m=6
test=data.frame(D_efficiency= (dets/max(dets))^{1/m},tour_length=long)
testbu=test
save(testbu,llPP,file=here("data/randpts.RData"))

show_front <- function(pref) {
  plot(test$D_efficiency,test$tour_length,xlab="Efficiency",ylab="Tour Length (km)")
  sky <- psel(test, pref)
  plot_front(test, pref, col = rgb(1, 0, 0))
  points(sky$D_efficiency, sky$tour_length, lwd = 3)
  points(test$D_efficiency[1],test$tour_length[1],pch='X',col='black',cex=1.5)
}


# do this for all four combinations of Pareto compositions
show_front(high(D_efficiency)  * low(tour_length))

dev.copy2pdf(file=here("data/SEPM_uc1_fig3a.pdf"))
dev.off()

######## Equivalent to figures 2.b-2.d but for the local D-optimal design (instead of the non-local)

####################
##### ellipses #####
####################

require(ellipse)

par(mfrow=c(1,3))
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

dev.copy2pdf(file=here("data/SEPM_uc1_Doptl_ellipses.pdf"))
dev.off()
