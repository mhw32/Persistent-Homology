# Separate distance function.
# Mike Wu
# To be used in the permutation tests.
source('NHST.r')
source('multiassign.r')
source('tools.r')
source('distance.r')

# The code below is based on the link:
# https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/

# Permutation Test
# ---------------------------------------------
# This is a wrapper function for the kernel and the energy test.
# It estimates the p-value through randomization. Asymptotically correct.

# We need a global sort of permutation test.
# (A generation of nhst in NHST.r)
# N = number of permutations to do.
permutationTest <- function(N, X, L, h=NULL) {
  # If no prior on h, find MLE.
  if (h == NULL)
    h <- findBestParam(X, L)
  Z <- 0 # Initialization
  loss_orig <- kernelStat(X, L, h)
  # Preserve some amount of order.
  order <- c(1:length(X))
  filler <- lapply(c(1:N), function(i) {
    random <- sample(order)
    # Only vary the order of the labels!
    loss_new <- kernelStat(X, L[random], h)
    if (loss_new <= loss_orig)
      Z <<- Z + 1
  })
  # Average the outputs.
  Z <- Z / N
  return(Z)
}

# Kernel Test
# ---------------------------------------------
# d is the distance between x and y.
# Use bottleneck or wasserstein.
# P, Q are two persistence diagrams
# DISTFUN is a distance function
# h is a smoothing parameter.
gaussianKernel <- function(P, Q, distfunc, h) {
  k <- exp(-distfunc(P, Q)^2 / h^2)
  return(k)
}

# X are the merged persistence diagrams.
# L are the merged persistence labels.
# h is the smoothing parameter
kernelStat <- function(X, L, h) {
  # Split the persistence diagrams into 2 groups.
  G1 <- X[L == 0]
  G2 <- X[L == 1]
  n <- length(G1)
  m <- length(G2)
  # To avoid loops, create a grid to loop over.
  g(nn_grid, nm_grid, mm_grid) %=% list(permutate(n, n), permutate(n, m), permutate(m, m))
  g(nn_grid_count, nm_grid_count, mm_grid_count) %=% list(length(nn_grid[[1]]), length(nm_grid[[1]]), length(mm_grid[[1]]))
  # There are three parts, calculate them separately.
  g(sum1, sum2, sum3) %=% list(0, 0, 0)
  # For each of the double for loops, loop through the grid.
  tmp <- lapply(seq(1:nn_grid_count), function(i) {
    sum1 <<- sum1 + gaussianKernel(G1[[nn_grid[i, 1]]], G1[[nn_grid[i, 2]]], contourDist, h)
  })
  tmp <- lapply(seq(1:nm_grid_count), function(i) {
    sum2 <<- sum2 + gaussianKernel(G1[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], contourDist, h)
  })
  tmp <- lapply(seq(1:mm_grid_count), function(i) {
    sum3 <<- sum3 + gaussianKernel(G2[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], contourDist, h)
  })
  # Now that we have all out sums, calculate the Tstat.
  T <- 1/(n^2)*sum1 - 2/(m*n)*sum2 + 1/(m^2)*sum3
  return(T)
}

# What is the best value for parameter h?
# IDK, so do a dumb search over a predefined space and take the suprenum.
findBestParam <- function(X, L, min=0, max=5, by=0.1) {
  hspace <- seq(min, max, by=by)
  v <- sapply(hspace, function(i) {
    t <- kernelStat(X, L, i)
    return(c(i, t))
  })
  # Suprenum is the max(abs())
  maxidx <- which.max(abs(v[2,]))
  maxh <- v[1,][maxidx]
  return(maxh)
}
