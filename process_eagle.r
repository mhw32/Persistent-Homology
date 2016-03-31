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
  h <- c(0, 100/3, 2*100/3, 100)
  for (i in 2:4) {
    for (j in 2:4) {
      for (k in 2:4) {
        dim1 <- c(h[i-1], h[i])
        dim2 <- c(h[j-1], h[j])
        dim3 <- c(h[k-1], h[k])
        logic <- (cube[,1] > dim1[1] & cube[,1] <= dim1[2]) & (cube[,2] > dim2[1] & cube[,2] <= dim2[2]) & (cube[,3] > dim3[1] & cube[,3] <= dim3[2])
        newcube <- cube[logic,]
        # Renormalize everything
        newcube[,1] <- newcube[,1] - (h[i-1])
        newcube[,2] <- newcube[,2] - (h[j-1])
        newcube[,3] <- newcube[,3] - (h[k-1])
        slices[[slicenum]] <- newcube
        slicenum <- slicenum + 1
      }
    }
  }
  return(slices)
}

slice_cube_robust <- function(cube, n) {
  # This is the final slices.
  slices <- vector('list', n^3)
  slicenum <- 1
  # initialize some item/counters
  item <- 0
  counter <- 1
  # create an array of splits
  h <- rep(NA, n+1)
  while (item < 99) {
    h[counter] <- item
    item <- item + 100/n
    counter <- counter + 1
  }
  h[counter] <- 100
  # Loop through each dimension
  for (i in seq(2,n+1)) {
    for (j in seq(2,n+1)) {
      for (k in seq(2,n+1)) {
        dim1 <- c(h[i-1], h[i])
        dim2 <- c(h[j-1], h[j])
        dim3 <- c(h[k-1], h[k])
        # Get a boolean array of each splice.
        logic <- (cube[,1] > dim1[1] & cube[,1] <= dim1[2]) & (cube[,2] > dim2[1] & cube[,2] <= dim2[2]) & (cube[,3] > dim3[1] & cube[,3] <= dim3[2])
        newcube <- cube[logic,]
        # Renormalize everything (0, 100/n)
        newcube[,1] <- newcube[,1] - (h[i-1])
        newcube[,2] <- newcube[,2] - (h[j-1])
        newcube[,3] <- newcube[,3] - (h[k-1])
        slices[[slicenum]] <- newcube
        slicenum <- slicenum + 1
      }
    }
  }
  return(slices)
}

# With the set, create persistence diagrams from each one.
persistify_set <- function(sampleset, n) {
  diagrams <- lapply(sampleset, function(set) {
    res <- 1.5
    boxlim <- c(0,100/n)
    diag <- gridDiag(set, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    return(diag$diagram)
  })
  return(diagrams)
}
