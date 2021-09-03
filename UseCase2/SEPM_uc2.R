###########################################################
### Spatialy Explicit Population Model : use case 2   #####
###########################################################

library(here)

# for solving PDE numericaly
require(ReacTran)
# for GIS related data 
require(raster)
require(maptools)
gpclibPermit()

######################
### SEPM parameters ##
######################

## Map for carrying capacity and (also) initial conditions ##
## as in https://doi.org/10.1111/emr.12350

land=raster(here("inputs/gis/wh_map_f.tif"))
# resize for computational purpose
tr=raster(ncol=100,nrow=100,xmn=0,xmx=1,ymn=0,ymx=1)
land=resample(land,tr,meth='ngb')

# Map for transects 
transects=readShapeLines(here("inputs/gis/wh_transects_f.shp"))
# cleaning up
transects@data=data.frame(id=1:length(transects@lines))

## nominal parameters ##
## as in https://doi.org/10.1111/emr.12350 (low growth and dispersal scenario)

params=list("b"=0.16,"D"=2,"mu"=0.1)
sig_uc2=10**(-3)
## Population initial conditions ##
## as in https://doi.org/10.1111/emr.12350 (supp mat)

klowmed=0.141
khigh=2.462
  
xini=land
xini[xini == 0]<-klowmed
xini[xini == 1]<-khigh
K=xini
xini=.9*xini ## modif 0.9*

## Nominal experimental Design ## 
# as in https://doi.org/10.1111/emr.12350 

# Number of transects : 32 
# Goal : choose only half 

transects_r=list()
for (i in 1:length(transects$id))
{
  temp=transects[ transects$id == i,]
  transects_r[[i]]  = rasterize(temp,land)
}

######################
### SEPM functions ###
######################

Nx <- Ny <- dim(xini)[1]
xgrid <- setup.grid.1D(x.up = 0, L = Nx, N = Nx); x <- xgrid$x.mid
ygrid <- setup.grid.1D(x.up = 0, L = Ny, N = Ny); y <- ygrid$x.mid

horses2D <- function(t, C, parms, K) {
  X <- matrix(nrow = Nx, ncol = Ny, data = C[1 : (Nx*Ny)]);
  dX <- tran.2D (C = X/K, D.x = parms$D, D.y = parms$D, 
                 dx = xgrid, dy = ygrid)$dC +  (parms$b * (1 - X/K) - parms$mu ) * X ;
  list(c(dX))
}

Tfin=7.0 

# ## unit test ##
# a=proc.time()
# out <- ode.2D (y = c(as.matrix(xini)), times = 1:Tfin, parms = params, 
#                func = horses2D, names = c("X"), dimens = c(Nx, Ny), method = "ode45",
#                K=as.matrix(K))
# b=proc.time()

set.seed(42)
ix.sC=sample(1:32,16)

obs_eq_horses2D<-function(params,xi=xini,ix.t=Tfin,ix.s=ix.sC,cap=K,tr_r=transects_r)
{ ## modif times = 1:ix.t
  out <- ode.2D (y = c(as.matrix(xi)), times = seq(0,ix.t,0.05), parms = params, 
                 func = horses2D, names = c("X"), dimens = c(Nx, Ny), method = "ode45",
                 K=as.matrix(cap))
  # X(t=7,x)
  Xtf=subset(out,which="X", subset = time %in% ix.t)
  
  # use raster syntaxic sugars to get observations 
  rXtf=raster(matrix(Xtf,Nx,Ny))
  obs=c()
  for (j in ix.s)
  {
    obs=c(obs,sum(as.vector(rXtf * tr_r[[j]]),na.rm=T))
  }
  return(obs)
}

# unit test 
observations=obs_eq_horses2D(params=params)
