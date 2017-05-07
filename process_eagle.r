# Tests with EAGLE simulation data.
library(rhdf5)
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

load_WDM_vanilla <- function() {
  d <- t(h5read("./simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloPositions"))
  return(d)
}

load_CDM_masses <- function() {
  d <- t(h5read("./simulations/Output_Eagle_Volume.hdf5", "P1/SubhaloMasses"))
  return(d)
}

load_WDM_masses <- function() {
  d <- t(h5read("./simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloMasses"))
  return(d)
}

# remove the least massive particles
load_downsampled_CDM_by_mass <- function() {
  cdm <- load_CDM()
  wdm <- load_WDM()
  cdm_masses <- load_CDM_masses()
  # get sizes
  num_wdm <- dim(wdm)[1]
  # sort the cdm by size
  indexes <- sort(cdm_masses, index.return=TRUE)
  cdm <- cdm[rev(indexes[["ix"]])[1:num_wdm],]
  return(cdm)
}

load_downsampled_CDM_by_mass_cuts <- function(threshold) {
  stopifnot(!(threshold %in% c(1, 0.5, 0.1)))
  if (threshold == 1) {
    wdm_mass_cut <- 1e8
    cdm_mass_cut <- 8e8
  } else if (threshold == 0.5) {
    wdm_mass_cut <- 1e9
    cdm_mass_cut <- 2e8
  } else {
    wdm_mass_cut <- 1e10
    cdm_mass_cut <- 1e10
  }

  cdm <- load_CDM()
  wdm <- load_WDM()
  cdm_max_masses <- load_CDM_max_masses()
  wdm_max_masses <- load_WDM_max_masses()
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
  indexmap <- vector('list', n^3)
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
        indexmap[[slicenum]] <- rbind(dim1, dim2, dim3)
        slicenum <- slicenum + 1
      }
    }
  }
  return(slices)
}

regroup_cube_robust <- function(slices, indexmap, n) {
  num_row <- 0
  for (i in 1:64) {
    num_row <- num_row + dim(slices[[i]])[1]
  }
  counter <- 1
  cube <- matrix(NA, num_row, 3)
  for (i in 1:64) {
    slice <- slices[[i]]
    map <- indexmap[[i]]
    num_in_slice <- dim(slice)[1]
    for (j in 1:3) {
      slice[,j] <- slice[,j] + map[j,1]
    }
    cube[counter:(counter+num_in_slice-1),] <- slice
    counter <- counter + num_in_slice
  }
  return(cube)
}


# With the set, create persistence diagrams from each one.
persistify_set <- function(sampleset, n, res=1.5) {
  diagrams <- lapply(sampleset, function(set) {
    boxlim <- c(0, 100/n)
    diag <- gridDiag(
      set, dtm, lim=cbind(boxlim, boxlim, boxlim),
      by=res, sublevel=T, printProgress=T, m0=0.001
    )
    return(diag$diagram)
  })
  return(diagrams)
}
