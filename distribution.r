# Direct Distribution Functions.
# Mike Wu

# Method: Parametrize the data. Because I can't really find a good multivariate non-parametric test, I can try to project the 2-dimension data into 1-dimension. I want to preserve relative distances between points. We note that by construction, the points in the persistence diagram follow a diagonal line. Therefore, we can rotate the plot and put the diagonal in the x-axis. (birth + death) / 2 and (birth - death) / 2.
# This is actually a matrix transformation. We also need to scale the standard deviation properly: A <- 1/sqrt(2)*matrix(c(1,1,1,-1),2,2). This should preserve all distances and essentially equals SVD. Then we can remove the 2nd dimension and be left with a scaled dimension. This is the first principle component!

source('tools.r')
source('distance.r')
library(MASS)

# 2-to-1 Dimensionality reduction (SVD).
reduce <- function(mat) {
  A <- 1/sqrt(2)*matrix(c(1,1,1,-1),2,2)
  z <- mat %*% A # Reparametrization
  return(z[,1]) # z --> The diagonal (should represent all the info)
}

# -------------------------------------------------
# [DEPRECATED] Given a persistence diagram. Pick out a 
# single dimension. Reduce dimension. Run a Manning Whitney 
# U test (nonparametric 1-d).
distribDimStat <- function(set, dim) {
  setnum <- length(set)
  reducedset <- lapply(1:length(set), function(i) {
    # Reduce specific dimension(s).
    reduce(sliceDim(set[[i]], dim))
  })
  return(reducedset)
}

# -------------------------------------------------
# Contour Test. This returns the density estimate.
contourDimStat <- function(set, dim) {
  setnum <- length(set)
  sliceset <- lapply(1:length(set), function(i) {
    input <- sliceDim(set[[i]], dim)
    inputx <- input[,1]
    inputy <- input[,2]
    # 2D KDE to get the density.
    return(kde2d(inputx, inputy)$z)
  })
  return(sliceset)
}

# -------------------------------------------------
# Persistence Intensity Functions
# -- uses weighting (based on height)
intensityDimStat <- function(set, dim, tau) {
  setnum <- length(set)
  sliceset <- lapply(1:length(set), function(i) {
    input <- sliceDim(set[[i]], dim)
    return(intensityDiagFunc(input, dim, tau))
  })
  return(sliceset)
}

intensityDiagFunc <- function(diag, dim, tau) {
  input <- sliceDim(diag, dim)
  xvec <- input[,1]
  yvec <- input[,2]
  xrange <- seq(min(xvec), max(xvec), by=0.01)
  yrange <- seq(min(yvec), max(yvec), by=0.01)
  numx <- length(xrange)
  numy <- length(yrange)

  stats <- as.matrix(rep(0, numx*numy), nrow=numx, ncol=numy)
  for (i=1:numx) {
    for (j=1:numy) {
      x <- xrange(i)
      y <- yrange(j)
      stats[i,j] <- intensityeqn(x, y, xvec, yvec, tau)
    }
  }
  return(stats)
}

intensityeqn <- function(x, y, births, deaths, tau) {
  z <- rep(0, length(x))
  f <- function(x, y, b, d, tau) {
    ans <- (d - b) * (1 / tau^2) * gaussianKernel((x - b) / tau) * gaussianKernel((y - d) / tau)
    return(ans)
  }
  loopans <- mapply(f, x=x, y=y, b=births, d=deaths, tau=tau)
  return(sum(loopans))
}

gaussianKernel <- function(x, y, h=1) {
  k <- exp(-euclidDist(x, y)^2 / h^2)
  return(k)
}

# -------------------------------------------------
# Persistence Images
# -- uses weighting (based on height)
pimageDimStat <- function(set, dim) {
}

# -------------------------------------------------
# [DEPRECATED] Global KDE Test.
globalDimStat <- function(set, dim) {
  setnum <- length(set)
  sliceset <- lapply(1:length(set), function(i) {
    input <- sliceDim(set[[i]], dim)
  })
  return(sliceset)
}
