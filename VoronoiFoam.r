#rm(list=ls(all=TRUE))
library(TDA)
library(scatterplot3d)
source("voronoi3dfct.r")
# -------------------------------------
# Generate VF Data
# -------------------------------------
# INPUT of voronoi3d function
# -------------------------------------
# Boxlim: limit of the side of the box
# resolution: grid space for approximating the voronoi cells.
# perturb: sd of gaussian noise to perturb the grid
# Ncells: number of voronoi cells
# N: total number of particles
# percClutter: percentage of clutter noise
# percWall: percentage of particles on the walls
# percFil: percentage of particles on the filaments
# percClust: percentage of particles on the clusters

# Example of a call to voronoi.
voronoi_example <- function() {
  Boxlim <- c(0,20)
  Xlim <- Boxlim
  Ylim <- Boxlim
  Zlim <- Boxlim
  resolution <- 0.5  ## grid space for approximating the voronoi cells.
  perturb <- 1 ## variance around the filaments
  N <- 1000   ## number of particles

  percFil <- 0.1
  vf1 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil, percClust=0.02)

  diag1 <- gridDiag(vf1, dtm, lim=cbind(Xlim,Ylim,Zlim), by=resolution, sublevel=T, printProgress=T, m0=0.001)

  scatterplot3d(vf1, pch=19, cex.symbol=.5, xlab="", ylab="", zlab="")
  plot(diag1$diagram)
}

# Given a voronsoi diagram, remove the first one of the 0th homologies, because it is always going to be infinity and thereby distracting. 
# This somehow messes up plotting.
clean <- function(diag) {
  infinity <- diag[1,]
  if (infinity[[1]] == 0) { diag <- diag[2:length(diag[,1]),] }
  return(diag)
}

# Function to create diagonal matrixes. 
voronoi_set <- function(percFil, N=1000, G=15, res=0.5, err=1, boxlim=c(0,10)) {
  set <- sapply(seq(1:G), function(i) {
    vf <- voronoi3d(boxlim, res, err, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil, percClust=0.02)
    diag <- gridDiag(vf, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    return(diag$diagram)
  })
}

voronoi_compilation <- function() {
  # Set a pretty big scope. 
  Boxlim <- c(0,50)
  resolution <- 0.5  # grid space for approximating the voronoi cells.
  perturb <- 1 # variance around the filaments.
  N <- 10000   # number of particles.
  groupN <- 15 # size of each set.
  # Do something with it.
  percFils <- seq(from=0.1, to=0.9, by=0.1)
  numSet <- length(percFils)
  storage <- vector("list", numSet)
  for (i in 1:numSet) {
    currSet <- voronoi_set(percFils[i], N, groupN, resolution, perturb, Boxlim)
    storage[[i]] = currSet
  }
  saveRDS(storage, "./voronoifoamfull.rds")
}

voronoi_baseline <- function() {
  Boxlim <- c(0, 50)
  resolution <- 0.5
  perturb <- 1
  N <- 10000
  groupN <- 15
  baseline <- voronoi_set(0.1, N, groupN, resolution, perturb, Boxlim)
  saveRDS(baseline, "./voronoibaseline.rds")
}








