#rm(list=ls(all=TRUE))
library(TDA)
library(scatterplot3d)
source("Voronoi3Dfct.r")
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
voronoi_example <- function(
  boxlim, percFil=0.1, res=0.5, perturb=1, N=10000
) {
  Xlim <- boxlim
  Ylim <- boxlim
  Zlim <- boxlim
  vf<- voronoi3d(
    boxlim, res, perturb, Ncells=64, N, percClutter=0,
    percWall=1-0.02-percFil, percFil=percFil, percClust=0.02
  )
  # scatterplot3d(
  #   vf, pch = 19, cex.symbol = .5,
  #   xlab = "", ylab = "", zlab = "", main=main
  # )
  diag <- gridDiag(
    vf, dtm, lim=cbind(Xlim,Ylim,Zlim), by=res,
    sublevel=T, printProgress=T, m0=0.001
  )
  return(diag)
}

voronoi_plot <- function(
  boxlim, percFil=0.1, res=0.5, perturb=1, N=10000, main="", path=""
) {
  Xlim <- boxlim
  Ylim <- boxlim
  Zlim <- boxlim
  vf<- voronoi3d(
    boxlim, res, perturb, Ncells=64, N, percClutter=0,
    percWall=1-0.02-percFil, percFil=percFil, percClust=0.02
  )
  png(filename=path)
  scatterplot3d(
    vf, pch = 19, cex.symbol = .5,
    xlab = "", ylab = "", zlab = "", main=main
  )
  dev.off()
}

# Function to create diagonal matrixes.
voronoi_set <- function(percFil, N=1000, G=15, res=0.5, err=1, boxlim=c(0,10)) {
  set <- sapply(seq(1:G), function(i) {
    vf <- voronoi3d(
      boxlim, res, err, Ncells=64, N, percClutter=0,
      percWall=1-0.02-percFil, percFil=percFil, percClust=0.02
    )
    diag <- gridDiag(
      vf, dtm, lim=cbind(boxlim, boxlim, boxlim),
      by=res, sublevel=T, printProgress=T, m0=0.001
    )
    return(diag$diagram)
  })
}

# Function to generate a foam.
voronoi_compilation <- function(
  N, Boxlim, res=0.5, perturb=1, groupN=15,
  baseline=0.1, nameId=1, folder='.'
) {
  percFils <- seq(from=0.1, to=0.3, by=0.05)
  numSet <- length(percFils)
  storage <- vector("list", numSet)
  for (i in 1:numSet) {
    currSet <- voronoi_set(percFils[i], N, groupN, res, perturb, Boxlim)
    storage[[i]] = currSet
  }
  if (substr(folder, nchar(folder), nchar(folder)+1) != "/") {
    folder <- paste(folder, "/", sep="")
  }
  saveRDS(storage, paste(folder, "foam", toString(nameId), ".rds", sep=""))
  # Run and save the baseline.
  baseline <- voronoi_set(baseline, N, groupN, res, perturb, Boxlim)
  saveRDS(baseline, paste(folder, "baseline", toString(nameId), ".rds", sep=""))
}

# ---------------------------------------------------

voronoi_only_set <- function(
  percFil, N=1000, G=15, res=0.5, err=1, boxlim=c(0,10)
) {
  set <- sapply(seq(1:G), function(i) {
    vf <- voronoi3d(
      boxlim, res, err, Ncells=64, N, percClutter=0,
      percWall=1-0.02-percFil, percFil=percFil, percClust=0.02
    )
    return(vf)
  })
}

voronoi_only_compilation <- function(
  N, Boxlim, res=0.5, perturb=1, groupN=15,
  baseline=0.1, nameId=1, folder='.'
) {
  percFils <- seq(from=0.1, to=0.3, by=0.05)
  numSet <- length(percFils)
  storage <- vector("list", numSet)
  for (i in 1:numSet) {
    currSet <- voronoi_only_set(percFils[i], N, groupN, res, perturb, Boxlim)
    storage[[i]] = currSet
  }
  if (substr(folder, nchar(folder), nchar(folder)+1) != "/") {
    folder <- paste(folder, "/", sep="")
  }
  saveRDS(storage, paste(folder, "foam", toString(nameId), ".rds", sep=""))
  # Run and save the baseline.
  baseline <- voronoi_only_set(baseline, N, groupN, res, perturb, Boxlim)
  saveRDS(baseline, paste(folder, "baseline", toString(nameId), ".rds", sep=""))
}





