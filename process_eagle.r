# Tests with EAGLE simulation data.
library(rhdf5)
library(rgl)
library(scatterplot3d)
library(TDA)

load_eagle <- function() {
  d = t(h5read("./simulations/Output_Eagle_Volume.hdf5", "P1/SubhaloPositions"))
  return(d)
}

# Generate a single dataset.
generate_one_sample <- function(dataset, samplenum) {
  idx = sample(1:nrow(dataset), samplenum, replace=FALSE)
  return(dataset[idx,])
}

# Create a list of setnum of samples from the eagle structure.
generate_sample_set <- function(dataset, setnum, samplenum) {
  # Now we can generate a bunch of them.
  set <- lapply(1:setnum, function(i) { generate_one_sample(dataset, samplenum) })
  return(set)
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
