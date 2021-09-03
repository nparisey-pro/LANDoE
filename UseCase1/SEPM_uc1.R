###########################################################
### Spatialy Explicit Population Model : use case 1   #####
###########################################################

library(here)

source(here('UseCase1/SIPM_uc1.R'))

# for solving PDE numericaly
require(ReacTran)
# for GIS related data 
require(raster)

######################
### SEPM parameters ##
######################

# hash table (for mapping landcovers to growth rates)
hashtab=c("rc"=1,"rs"=2,"ru"=3,"rg"=4)

## landcover example ##
## as in doi:10.1007/s10144-013-0415-0 

land=raster(here('inputs/gis/landcover1km.grd'))

## nominal parameters ##
## as in doi:10.1016/j.ecocom.2016.07.003

Dp=0.0001896148
Dh=100^2 * Dp / 6.3 
parmsPDE=parms
parmsPDE$alpha=10^-4
parmsPDE$D=Dh

## Population initial conditions ##
# Created using packages NLMR & landscapetools 
# Packages not included because of there many dependencies 

uini=raster(here("inputs/gis/initcond.grd"))
vini <-uini*0; 

## Nominal experimental Design ## 
# as in doi:10.1016/j.ecocom.2016.07.003
# Number of traps : 3 traps per field ; we choose here 2 fields of different type
# Sampling times : 9 samplings; traps open on mondays, sampled fridays, closed the week-ends

spl.d=sort(c(seq(0,61,7),seq(5,61,7)))
# sampled : upper left field and lower right field 
set.seed(57)
spl.pts=cbind(c(runif(3,0.15,0.3),runif(3,0.8,1)),
              c(runif(3,0.04,0.15),runif(3,0.2,0.4)))

######################
### SEPM functions ###
######################

Nx <- Ny <- dim(uini)[1]
xgrid <- setup.grid.1D(x.up = 0, L = Nx, N = Nx); x <- xgrid$x.mid
ygrid <- setup.grid.1D(x.up = 0, L = Ny, N = Ny); y <- ygrid$x.mid

carabid2D <- function(t, C, parms,R) {
  u <- matrix(nrow = Nx, ncol = Ny, data = C[1 : (Nx*Ny)]);
  v <- matrix(nrow = Nx, ncol = Ny, data = C[(Nx*Ny+1) : (2*Nx*Ny)]);
  du <- tran.2D (C = u, D.x = parms$D, D.y = parms$D, 
                 dx = xgrid, dy = ygrid)$dC + ( R * exp(-1*parms$beta*t) - parms$mu) * u  ;
  dv <- parms$alpha * u; 
  list(c(du,dv))
}

contextcarab2D<-function (params,ui=uini,vin=vini,
                ix.t=spl.d,landcov=land,hash=hashtab)
{
  # mapping lancover to growth rates 
  values(landcov)[ values(landcov) == hash['ru'] ] <- 0
  for (what in c('rc','rs','rg'))
  {    values(landcov)[ values(landcov) == hash[what] ] <- params[[what]]   }
  R=as.matrix(landcov)
  
  out <- ode.2D (y = c(as.matrix(ui),as.matrix(vin)), times = ix.t, parms = params, 
                 func = carabid2D, names = c("u","v"), dimens = c(Nx, Ny), method = "ode45",
                 R=R)
  return(out)
}
 
# ## unit test ##
# a=proc.time()
# out=contextcarab2D(params=parmsPDE,ix.t=spl.d)
# print(proc.time()-a)  
# # some visuals 
# image(out, which = "u", subset = time %in% seq(5,61,7), grid = list(x, y), mfrow = c(3,3), ask = FALSE,zlim=c(0,0.5*10^5))

obs_eq_CC2D<-function(params,ui=uini,vin=vini,
                      ix.t=spl.d,ix.s=spl.pts,landcov=land,hash=hashtab)
{
  out<-contextcarab2D(params,ui,vin,ix.t,landcov,hash)
  
  # extract results for times 't' and 't-tau'
  ts_mtau=subset(out,which="v", subset = time %in% spl.d[seq(1,length(spl.d),2)])
  ts=subset(out,which="v", subset = time %in% spl.d[seq(2,length(spl.d),2)])
  
  # compute definite integral int^{t}_{t-\tau} over the whole spatial domain 
  obs_eq=list()
  val=c()
  for (i in 1:dim(ts)[1])
  {
    obs_eq[[i]]=matrix(as.vector(ts[i,]),Nx,Ny) - matrix(as.vector(ts_mtau[i,]),Nx,Ny) 
    
    # extract values at sampling points (x,y)
    for (j in 1:dim(ix.s)[1])
    {
      val=c(val,obs_eq[[i]][round(ix.s[j,1]*Nx),round(ix.s[j,2]*Ny)])
    }  
  }
  return(val)
}

# # unit test
# a=proc.time()
# wut=obs_eq_CC2D(params=parmsPDE)
# print(proc.time()-a) 
