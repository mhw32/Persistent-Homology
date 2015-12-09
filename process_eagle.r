# Tests with EAGLE simulation data.
library(rhdf5)
library(rgl)
library(scatterplot3d)
library(TDA)

load_CDM <- function() {
  d <- t(h5read("./simulations/Output_Eagle_Volume.hdf5", "P1/SubhaloPositions"))
  return(d)
}

# Increase a processing step to ensure only genuine objects.
load_WDM <- function() {
  d <- t(h5read("./simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloPositions"))
  HMspher <- t(h5read("./simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloHMSpher"))
  MaxMass <- t(h5read("./simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloMaxMass"))
  selection <- (MaxMass > 2.2e8) & (HMspher > 0.2)
  reframe <- d[selection,]
  return(reframe)
}

# we can't just take a sample. Instead we have to divide the cube itself into 27 smaller sets.
# A cube is cubic, so n=3 --> 3 in each dimension = 27 slices.
slice_cube <- function(cube) {
  slices <- vector('list', 27)
  slicenum <- 1
  h <- c(1, 33, 66, 100)
  for (i in 2:4) {
    for (j in 2:4) {
      for (k in 2:4) {
        dim1 <- h[i-1:i]
        dim2 <- h[j-1:j]
        dim3 <- h[k-1:k]
        logic <- (cube[,1] >= dim1[1] & cube[,1] <= dim1[2]) & (cube[,2] > dim2[1] & cube[,2] <= dim2[2]) & (cube[,3] > dim3[1] & cube[,3] <= dim3[2])
        slices[[slicenum]] <- cube[logic,]
        slicenum <- slicenum + 1
      }
    }
  }
  return(slices)
}


# With the set, create persistence diagrams from each one.
persistify_set <- function(sampleset) {
  diagrams <- lapply(sampleset, function(set) {
    res <- 0.5
    boxlim <- c(0,10)
    diag <- gridDiag(set, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    return(diag$diagram)
  })
  return(diagrams)
}
