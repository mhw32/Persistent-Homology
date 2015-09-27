# Author: Mike Wu
# Playing around a little bit more with persistent homology.
library(TDA)

# -------------------------------------------------------------------------
# TRIAL ONE : Given N simulations of n particle datasets, generate and overlay persistent diagrams. 

# We need a function to generate data from a uniform circle.
# Inputs : n [number to sample], r [radius of the circle], c [center of the circle].
genCircle <- function(num=n, rad=r, center=c(0, 0)) {
  # Define center of the radius.
  x0     <- center[1] 
  y0     <- center[2]
  # Sample uniform points and put them on a point.
  u      <- 2 * pi * runif(num)
  # Return the circles. 
  result <- rad / 2 * cbind(x = cos(u) + x0, y = sin(u) + y0)
  return(result)
}

# Run the RIPS diagram for N independent simulation
multiRIPS <- function(maxdimension, maxscale, N, K) {
  for (i in 1:N) {
    Circle <- genCircle(K, radius, center)
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

# First way of adding noise. This shifts the actual points on the circle a little bit. 
noisyRIPS <- function(maxdimension, maxscale, N, K, emu, esig) {
  for (i in 1:N) {
    Circle <- genCircle(K, radius, center)
    Noise1 <- makeNoise(emu, esig, K)
    Noise2 <- makeNoise(emu, esig, K)
    # Combine the vectors with the noise
    Circle[, 1] = Circle[, 1] + Noise1
    Circle[, 2] = Circle[, 2] + Noise2
    # Add some gaussian noise to the circle.
    Diag   <- ripsDiag(Circle, maxdimension, maxscale, library = "GUDHI", printProgress = FALSE)
    # Use points so that we can overlay these RIPS diagrams.
    par(new = TRUE)
    plot(Diag$diagram, main = paste("RIPS Diagram", "( N =", toString(N), ", n =", toString(K), "noise = Gaussian(", toString(emu), ",", toString(esig), ")", ")"
))
  }
}

# This adds noise a second way. This adds additional points themselves to the circle.
# But it lets circles maintain its structure.
noisy2RIPS <- function(maxdimension, maxscale, N, K, E) {
  for (i in 1:N) {
    Circle <- genCircle(K, radius, center)
    Noise1 <- makeNoise(center[1], radius/2, E)
    Noise2 <- makeNoise(center[1], radius/2, E)
    Noise <- cbind(Noise1, Noise2)
    # Combine the vectors with the noise
    Circle <- rbind(Circle, Noise)
    # Add some gaussian noise to the circle.
    Diag   <- ripsDiag(Circle, maxdimension, maxscale, library = "GUDHI", printProgress = FALSE)
    # Use points so that we can overlay these RIPS diagrams.
    par(new = TRUE)
    plot(Diag$diagram, main = paste("RIPS Diagram", "( N =", toString(N), ", n =", toString(K), "noise = Gaussian(", toString(center[1]), ",", toString(radius/2), ")", ")"
    ))
  }
}

# Example running setting
maxdimension <- 1          # components and loops
maxscale     <- 5          # limit of the filtration
numparticle  <- 10         # number of samples to draw per circle
radius       <- 5          # radius of circle
center       <- c(0, 0)    # center of circle
multiRIPS(maxdimension, maxscale, 10, 250)

# ----------------------------------------------------------------
# TRIAL THREE : Do the same thing with "Distance to Measure" RIPS.
multiDTM <- function(N, K, Xlim, Ylim, by, m0=0.1) {
  for (i in 1:N) {
    Circle <- genCircle(K, radius, center)
    # Use gridDiag with the Kernal Density Estimator.
    Diag <- gridDiag(X = Circle, FUN = dtm, lim = cbind(Xlim, Ylim), by = by, 
                     sublevel = FALSE, library = "Dionysus", printProgress = FALSE, m0 = m0)
    par(new = TRUE)
    plot(Diag$diagram, main = paste("DTM Diagram", "( N =", toString(N), ", n =", toString(K), ")"))
  }
}

# Reformatting Function after reading TDA tutorial more carefully.
# ----------------------------------------------------------------

singleKDE <- function(Circle, Grid, h=0.3) {
  band=bootstrapBand(X=Circle, FUN=kde, Grid=Grid, B=100, parallel=FALSE, alpha=0.1, h=h)
  Diag = gridDiag( X=X, FUN=kde, h=0.3, lim=cbind(Xlim,Ylim), by=by, sublevel=FALSE, library="Dionysus", printProgress=FALSE )$diagram
  return(list(first=band, second=Diag))
}

# Before we can do any distance estimates, we need to construct a grid of points over which we will evaluate the functions.
# Inputs : X is the sampled units from a circle
makeGrid <- function(X, Xlim, Ylim, by) {
  Xseq <- seq(from = Xlim[1], to = Xlim[2], by = by)
  Yseq <- seq(from = Ylim[1], to = Ylim[2], by = by)
  Grid <- expand.grid(Xseq, Yseq)
  return(Grid)
}

multiKDE <- function(N, K, Xlim=NULL, Ylim=NULL, by=0.065, rad=5, center=c(0,0), h=0.3, main="KDE Diagram") {
  par(new = TRUE)
  for (i in 1:N) {
    Circle <- genCircle(K, rad, center)
    if (is.null(Xlim)) {
      Xlim <- c(center[1] - rad - 1, center[1] + rad + 1)
    }
    if (is.null(Ylim)) {
      Ylim <- c(center[2] - rad - 1, center[2] + rad + 1) 
    }
    grid <- makeGrid(Circle, Xlim, Ylim, by)
    result <- singleKDE(Circle, grid, h)
    plot(result$second, band=2*result$first$width, main=main)
  }
}

