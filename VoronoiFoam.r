#rm(list=ls(all=TRUE))
library(TDA)
library(scatterplot3d)
source("voronoi3dfct.R")
# ########################################
# ## Generate VF Data
# ########################################
## INPUT of voronoi3d function
##############################
## Boxlim: limit of the side of the box
## resolution: grid space for approximating the voronoi cells.
## perturb: sd of gaussian noise to perturb the grid
## Ncells: number of voronoi cells
## N: total number of particles
## percClutter: percentage of clutter noise
## percWall: percentage of particles on the walls
## percFil: percentage of particles on the filaments
## percClust: percentage of particles on the clusters


Boxlim=c(0,50)
Xlim=Boxlim
Ylim= Boxlim
Zlim= Boxlim
resolution=0.5  ## grid space for approximating the voronoi cells.
perturb=1 ## variance around the filaments
N=5000   ## number of particles

percFil = 0.7
vf1 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil, percClust=0.02)

diag1 = gridDiag(vf1, dtm, lim = cbind(Xlim,Ylim,Zlim), by = resolution, sublevel=T, printProgress=T, m0=0.001)

scatterplot3d(vf1, pch = 19, cex.symbol = .5, xlab = "", ylab = "", zlab = "")
plot(diag1$diagram)










