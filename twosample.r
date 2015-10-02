# Separate distance function. 
# Mike Wu
# To be used in the permutation tests. 
source('/Users/grub/Desktop/Cisewski-Lab/NHST.r')
source('/Users/grub/Desktop/Cisewski-Lab/multiassign.r')
source('/Users/grub/Desktop/Cisewski-Lab/tools.r')
# The code below is based on the link:
# https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/

# Kernel Test
# ---------------------------------------------
# d is the distance between x and y. 
# Use bottleneck or wasserstein.
# P, Q are two persistence diagrams
# DISTFUN is a distance function
# h is a smoothing parameter.
gaussianKernel <- function(P, Q, DISTFUN, ...) {
  h <- c(...)[1] # This only needs 1 parameter.
  k <- exp(-DISTFUN(P, Q)^2 / h^2)
  return(k)
}

# X are the merged persistence diagrams.
# L are the merged persistence labels.
# KERNEL is the kernel we are using.
# ... is a placeholder for arbitrary number of parameters.
kernelStat <- function(X, L, KERNEL, DISTFUN, ...) {
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
    sum1 <<- sum1 + KERNEL(G1[[nn_grid[i, 1]]], G1[[nn_grid[i, 2]]], DISTFUN, ...)
  })
  tmp <- lapply(seq(1:nm_grid_count), function(i) {
    sum2 <<- sum2 + KERNEL(G1[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, ...)
  })
  tmp <- lapply(seq(1:mm_grid_count), function(i) {
    sum3 <<- sum3 + KERNEL(G2[[mm_grid[i, 1]]], G2[[mm_grid[i, 2]]], DISTFUN, ...)
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
energyKernel <- function(P, Q, DISTFUN, ...) {
  alpha <- c(...)[1]
  k <- DISTFUN(P, Q)^alpha
  return(k)
}

energyStat <- function(X, L, ENERGY, DISTFUN, ...) {
  # This is super similar to a kernel test and can be collapsed into it.  
  T <- kernelStat(X, L, ENERGY, DISTFUN, ...)
  # one-half of the harmonic mean of the sample sizes * -Tstat
  # https://cran.r-project.org/web/packages/energy/energy.pdf
  E <- (n*m/(n+m)) * -T 
  return(E)
}

# Permutation Test
# ---------------------------------------------
# This is a wrapper function for the kernel and the energy test.
# It estimates the p-value through randomization. Asymptotically correct.

# Given two vectors, X, Y, try to do a non-bipartite minimum match
# via a grid formulation.
minimatch <- function(X, Y) {
  g(n, m) %=% list(length(X), length(Y))
  Z <- c(X, Y)
  len <- length(Z)
  # Create a grid over the vectors. 
  grid <- permutate(len1, len2)
  gridlen <- length(grid[[1]])
  # Loop through the grid and calculate distance.
  distances <- sapply(seq(1:gridlen), function(i) {
    d <- DISTFUN(Z[[grid[i, 1]]], Z[[grid[i, 2]]])
    return(d)
  }) 
  # Calculate all distances and add it to the grid table.
  grid['Dist'] = distances
  # Find index of the minimum for each row item 1..len
  indices <- sapply(seq(1:len), function(i) {
    idx <- which.min(grid['Dist'][grid['Var2'] == i])
    return(c(i, idx))
  })
  # Tranpose the indices (this is the index matching).
  indices <- t(indices)
  return(indices)
}

# Rossenbaum Test (Cross-Match)
# ---------------------------------------------
# Unlike the other two tests, this one does not require a permutation wrapper. 
rosenbaumStat <- function(Z, L, DISTFUN) {
  g(n, m) <- list(length(Z[L == 0]), length(Z[L == 1]))
  allL <- c(L[L == 0], L[L == 1])
  # Non-bipartite matching across Z
  indices <- minimatch(Z[,1], Z[,2])
  # Abstract the indices and pull out the labels.
  labelidx <- cbind(allL[indices[,1]], allL[indices[,2]])
  # Count number bigger, same, smaller
  g(zeros, both, ones) %=% list(0, 0, 0)
  tmp <- lapply(seq(1:length(allL)), function(i) {
    if (labelidx[i, 1] == labelidx[i, 2]) {
      if (labelidx[i,1] == 1) {
        ones <<- ones + 1
      } else {
        zeros <<- zeros + 1
      }
    } else {
      both <<- both + 1
    }
  })
  # Define T as the number (0,1), (1,0) pairs.
  T <- both
  # Compute the exact distribution of T under H_0 (the probability)
  N <- zeros + ones + both
  P <- (2^T * factorial(N)) / (choose(N, m) * factorial(zeros) * factorial(both) * factorial(ones))
  return(P)
}




