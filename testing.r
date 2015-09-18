# Author: Mike Wu
# Playing around a little bit more with persistent homology.
library(TDA)

# -------------------------------------------------------------------------
# TRIAL ONE : Given N simulations of n particle datasets, generate and overlay persistent diagrams. 

# We need a function to generate data from a uniform circle.
# Inputs : n [number to sample], r [radius of the circle], c [center of the circle].
circle <- function(num=n, rad=r, center=c(0, 0)) {
  # Define center of the radius.
  x0     <- center[1] 
  y0     <- center[2]
  # Sample uniform points and put them on a point.
  u      <- 2 * pi * runif(num)
  # Return the circles. 
  result <- rad * cbind(x = cos(u) + x0, y = sin(u) + y0)
  return(result)
}

# Run the RIPS diagram for N independent simulation
multiRIPS <- function(maxdimension, maxscale, N, K) {
  for (i in 1:N) {
    Circle <- circle(K, radius, center)
    Diag   <- ripsDiag(Circle, maxdimension, maxscale, library = "GUDHI", printProgress = FALSE)
    # Use points so that we can overlay these RIPS diagrams.
    par(new = TRUE)
    plot(Diag$diagram, main = paste("RIPS Diagram", "( N =", toString(N), ", n =", toString(K), ")"))
  }
}

# Adds gaussian noise
makeNoise <- function(mu, sigma, num) {
  result <- rnorm(num, mean = mu, sd = sigma)
  return(result)
}

noisyRIPS <- function(maxdimension, maxscale, N, K, emu, esig) {
  for (i in 1:N) {
    Circle <- circle(K, radius, center)
    # Add some gaussian noise to the circle.
    Diag   <- ripsDiag(Circle, maxdimension, maxscale, library = "GUDHI", printProgress = FALSE)
    # Use points so that we can overlay these RIPS diagrams.
    par(new = TRUE)
    plot(Diag$diagram, main = paste("RIPS Diagram", "( N =", toString(N), ", n =", toString(K), "noise = N(", toString(emu), ",", toString(esig), ")", ")"
))
  }
}

# Example running setting
maxdimension <- 1          # components and loops
maxscale     <- 5          # limit of the filtration
numparticle  <- 10         # number of samples to draw per circle
radius       <- 2          # radius of circle
center       <- c(0, 0)    # center of circle
multiRIPS(maxdimension, maxscale, 10, 250)

# ---------------------------------------------------------
# TRIAL TWO : Do the same thing with "Kernel Density" RIPS.

# Before we can do any distance estimates, we need to construct a grid of points over which we will evaluate the functions.
# Inputs : X is the sampled units from a circle
makeGrid <- function(X, Xlim, Ylim, by) {
  Xseq <- seq(from = Xlim[1], to = Xlim[2], by = by)
  Yseq <- seq(from = Ylim[1], to = Ylim[2], by = by)
  Grid <- expand.grid(Xseq, Yseq)
  return(Grid)
}

# Example of creating a grid.
numparticle <- 10
radius      <- 2      
center      <- c(0, 0)
X           <- circle(numparticle, radius, center)
Xlim        <- c(center[1] - radius - 1, center[2] + radius + 1)
Ylim        <- c(center[1] - radius - 1, center[2] + radius + 1)
by          <- 0.1
Grid        <- makeGrid(X, Xlim, Ylim, by)

multiKDE <- function(N, K, Xlim, Ylim, by, h=0.3) {
  for (i in 1:N) {
    Circle <- circle(K, radius, center)
    # Use gridDiag with the Kernal Density Estimator.
    Diag <- gridDiag(X = Circle, FUN = kde, lim = cbind(Xlim, Ylim), by = by, 
                     sublevel = FALSE, library = "Dionysus", printProgress = FALSE, h = h)
    par(new = TRUE)
    plot(Diag$diagram, main = paste("KDE Diagram", "( N =", toString(N), ", n =", toString(K), ")"))
  }
}

# ----------------------------------------------------------------
# TRIAL THREE : Do the same thing with "Distance to Measure" RIPS.
multiDTM <- function(N, K, Xlim, Ylim, by, m0=0.1) {
  for (i in 1:N) {
    Circle <- circle(K, radius, center)
    # Use gridDiag with the Kernal Density Estimator.
    Diag <- gridDiag(X = Circle, FUN = dtm, lim = cbind(Xlim, Ylim), by = by, 
                     sublevel = FALSE, library = "Dionysus", printProgress = FALSE, m0 = m0)
    par(new = TRUE)
    plot(Diag$diagram, main = paste("DTM Diagram", "( N =", toString(N), ", n =", toString(K), ")"))
  }
}

