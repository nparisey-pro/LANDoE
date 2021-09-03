This archive contains code and data related to the article "Optimal spatial monitoring of populations described by reaction-diffusion models" by N. Parisey, M. Leclerc and K. Adamczyk-ChauvatK.

We recommend you to explore this code using the associated rstudio project ('LANDoE.Rproj') notably because relative paths are handled by the package 'here' and the project's structure.

The archive is structured as follow :

(0) Files 'LICENCE' adn README.md refer, respectively, to the GNU GPL Licence under which these scripts are distributed and this present guide ;

(i) Repository 'inputs/' contains geomatic maps ('inputs/gis/') used to describe spatial domains, initial conditions and traps positions for use cases and serialized objets ('inputs/others/') ;

(ii) Repository 'UseCase1/' contains all R scripts necessary to solve use case 1, notably codes for 
(ii).1 simulation of the spatialy implict population model (SIPM_uc1.R) and related observations ;
(ii).2 simulation of the spatialy explicit population model (SEPM_uc1.R) and related observations ;
(ii).3 resolution of the local D-optimal design ('SEPM_OptDesign_uc1_l.R') ; 
(ii).4 resolution of the non-local D-optimal design ('SEPM_OptDesign_uc1_nl.R').

(iii) Repository 'UseCase2/' contains all R scripts necessary to solve use case 2, notably code for
(iii).1 simulation of the spatialy explicit population model (SEPM_uc2.R) and related observations ;
(iii).2 resolution of the local D-optimal design ('SEPM_OptDesign_uc2_l.R') ;
(iii).3 resolution of the non-local D-optimal design ('SEPM_OptDesign_uc2_nl.R').

(iv) Repository 'data/' is for simple graphics and serialized objects (notably Fisher Matrices) generated during code run, it starts empty.

To test the code, one can source the scripts for resolution of local D-optimal designs as their computation is quite fast (a few minutes at worst). The non-local D-optimal designs take longer, a few dozen minutes at most (including solving the local designs as well, by default).


Required R packages:
deSolve, doParallel, ellipse, gstat, geospt, here, imager, lhs, maptools, Matrix, numDeriv, raster, ReacTran, rgdal, rgeos, rPref, spatstat, spdep, TSP

For information, code ran sucessfully on several Ubuntu desktops (16.04 and 18.04).
It was developped using :
R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
deSolve   	1.28
doParallel 	1.0.15
ellipse		0.4.1
gstat		2.0-6
geospt		1.0-2
here		1.0.1
imager		0.42.1
lhs 		1.0.2
maptools	1.0-1
Matrix		1.2-18
numDeriv	1.1
raster		3.0-12
ReacTran	1.4.3.1
rgdal		1.4-8
rgeos		0.5-3
rPref		1.3
spatstat	1.64-1
spdep		1.1-8
TSP		1.1-9

Contact : N. Parisey, nicolas.parisey@inrae.fr
