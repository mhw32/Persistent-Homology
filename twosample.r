# Separate distance function. 
# Mike Wu
# To be used in the permutation tests. 
source('/Users/grub/Desktop/Cisewski-Lab/NHST.r')
source('/Users/grub/Desktop/Cisewski-Lab/multiassign.r')
# The code below is based on the link:
# https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/

# Kernel Test
# ---------------------------------------------

# d is the distance between x and y. 
# Use bottleneck or wasserstein.
# P, Q are two persistence diagrams
# DISTFUN is a distance function
# h is a smoothing parameter.
gaussianKernel <- function(P, Q, DISTFUN, h) {
  k <- exp(-DISTFUN(P, Q)^2 / h^2)
  return(k)
}

# X are the merged persistence diagrams.
# L are the merged persistence labels.
# KERNEL is the kernel we are using
kernelStat <- function(X, L, KERNEL, DISTFUN, h) {
  # Split the persistence diagrams into 2 groups.
  G1 <- X[L == 0]
  G2 <- X[L == 1]
  n <- length(G1)
  m <- length(G2)
  # To avoid loops, create a grid to loop over.
  nn_grid <- permutate(n, n)
  nm_grid <- permutate(n, m)
  mm_grid <- permutate(m, m)
  nn_grid_count <- length(nn_grid[[1]])
  nm_grid_count <- length(nm_grid[[1]])
  mm_grid_count <- length(mm_grid[[1]])
  # There are three parts, calculate them separately.
  g(sum1, sum2, sum3) %=% list(0, 0, 0)
  # For each of the double for loops, loop through the grid.
  tmp <- lapply(seq(1:nn_grid_count), function(i) {
    sum1 <<- sum1 + KERNEL(G1[[nn_grid[i, 1]]], G1[[nn_grid[i, 2]]], DISTFUN, h)
  })
  tmp <- lapply(seq(1:nm_grid_count), function(i) {
    sum2 <<- sum2 + KERNEL(G1[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, h)
  })
  tmp <- lapply(seq(1:mm_grid_count), function(i) {
    sum3 <<- sum3 + KERNEL(G2[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, h)
  })
  # Now that we have all out sums, calculate the Tstat. 
  T <- 1/(n^2)*sum1 - 2/(m*n)*sum2 + 1/(m^2)*sum3
  return(T)
}

# What is the best value for parameter h?
# IDK, so do a dumb search over a predefined space and take the suprenum.
findBestParam <- function(X, L, KERNEL, DISTFUN, min=0, max=5, by=0.1) {
  hspace <- seq(min, max, by=by)
  v <- sapply(hspace, FUN=function(i) {
    t <- kernelStat(X, L, KERNEL, DISTFUN, i)
    return(c(i, t))
  })
  # Suprenum is the max(abs())
  maxidx <- which.max(abs(v[2,])) 
  maxh <- v[1,][maxidx]
  return(maxh)
}

# Energy Test
# ---------------------------------------------

energyKernel <- function(P, Q, DISTFUN, alpha=1) {
  k <- DISTFUN(P, Q)^alpha
  return(k)
}

# This is super similar to a kernel test and probably can be collapsed into it.
kernelStat <- function(X, L, ENERGY, DISTFUN, alpha=1) {
  # ------------------------------------------ 
  # Split the persistence diagrams into 2 groups.
  G1 <- X[L == 0]
  G2 <- X[L == 1]
  n <- length(G1)
  m <- length(G2)
  # To avoid loops, create a grid to loop over.
  nn_grid <- permutate(n, n)
  nm_grid <- permutate(n, m)
  mm_grid <- permutate(m, m)
  nn_grid_count <- length(nn_grid[[1]])
  nm_grid_count <- length(nm_grid[[1]])
  mm_grid_count <- length(mm_grid[[1]])
  # There are three parts, calculate them separately.
  g(sum1, sum2, sum3) %=% list(0, 0, 0)
  # For each of the double for loops, loop through the grid.
  tmp <- lapply(seq(1:nn_grid_count), function(i) {
    sum1 <<- sum1 + ENERGY(G1[[nn_grid[i, 1]]], G1[[nn_grid[i, 2]]], DISTFUN, alpha)
  })
  tmp <- lapply(seq(1:nm_grid_count), function(i) {
    sum2 <<- sum2 + ENERGY(G1[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, alpha)
  })
  tmp <- lapply(seq(1:mm_grid_count), function(i) {
    sum3 <<- sum3 + ENERGY(G2[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, alpha)
  })
  # ------------------------------------------
  # Now that we have all out sums, calculate the energy
  # First half : one-half of the harmonic mean of the sample sizes
  E <- (n*m/(n+m)) * (2/(m*n)*sum2 - 1/(n^2)*sum1 - 1/(m^2)*sum3)
  return(E)
}

# Rossenbaum Test
# ---------------------------------------------

bipartiteMatch <- function() {

}

rossenbaumStat <- function() {

}