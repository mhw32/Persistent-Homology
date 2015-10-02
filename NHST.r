# Author: Mike Wu
# Based on http://arxiv.org/pdf/1310.7467v1.pdf
# Randomization-style Null Hypothesis Significance Tests for Persistence Diagrams

library(TDA)
source('/Users/grub/Desktop/Cisewski-Lab/homology.r')

# Create all perms or 1..a, 1..b. 
# Dups are preserved.
permutate <- function(a, b) {
  A <- seq(1:a)
  B <- seq(1:b)
  C <- expand.grid(A, B)
  return(C)
}

# TODO: THIS NEEDS TO BE SPED UP.
# X is a vector of persistence diagrams. 
# L is a vector of labels. 
lossfunc <- function(X, L) {
  # Split the persistence diagrams into 2 groups.
  G <- list(X[L == 0], X[L == 1])
  n <- c(length(G[[1]]), length(G[[2]]))
  sum1 <- 0 # Initialize outer counter.
  fillme <- lapply(seq(1:2), function(m) {
    # Get all permutations for this n[m]
    tups <- permutate(n[m], n[m])
    num <- length(tups[[1]]) # Number of tuples
    P <- G[[m]]
    C <- 1 / (2 * n[m] * (n[m] - 1))
    sum2 <- 0 # Initialize a counter.
    filler <- lapply(seq(1:num), function(i) {
      sum2 <<- sum2 + bottleneck(P[[tups[i,1]]], P[[tups[i,2]]], dimension=1)
      # sum2 <- sum2 + wasserstein(P[[tups[i,1]]], P[[tups[i,2]]], dimension=1, p=2)
    })
    sum1 <<- sum1 + C * sum2
  });
  return(sum1)
}

# Given n1 + n2 persistence diagrams with n1 + n2 labels. Note that n1, and the following n2 diagrams must be disjoint. This is meshed together as X with labels L. Given N (number of iterations).
nhst <- function(X, L, N) {
  Z <- 0 # Initialization
  loss_orig <- lossfunc(X, L)
  # Preserve some amount of order.
  order <- c(1:length(X))
  filler <- lapply(c(1:20), function(i) {
    random <- sample(order)
    # Only vary the order of the labels!
    loss_new <- lossfunc(X, L[random])
    if (loss_new <= loss_orig)
      Z <<- Z + 1
  })
  # Average the outputs.
  Z <- Z / N
  return(Z)
}
# By law of large numbers, the expectation of Z --> P(loss(L_new) <= loss(L_obs)) as n --> infinity.

# Code for an example.
# -----------------------------------------------------------------
circle1 <- function(N, mu, sig) {
  K <- circleUnif(N, r=1)
  K <- addRandomNoise(K, mu, sig)
  filler <- lapply(c(1:20), function(i) {
    K2 <- circleUnif(N, r=1)
    K2 <- addRandomNoise(K2, mu, sig)
    K <<- rbind(K, K2)
  })
  return(K)
}

circle2 <- function(N, mu, sig) {
  L1 <- circleUnif(N/2, 3/5)
  L1[,1] <- L1[,1] - 2/5
  L2 <- circleUnif(N/2, 2/5)
  L2[,1] <- L2[,1] + 3/5
  L <- rbind(L1, L2)
  L <- addRandomNoise(L, mu, sig)
  filler <- lapply(c(1:20), function(i) {
    L1b <- circleUnif(N/2, 3/5)
    L1b[,1] <- L1b[,1] - 2/5
    L2b <- circleUnif(N/2, 2/5)
    L2b[,1] <- L2b[,1] + 3/5
    Lb <- rbind(L1b, L2b)
    Lb <- addRandomNoise(Lb, mu, sig)
    L <<- rbind(L, Lb)
  })
  return(L)
}

# For demonstration, let's just do the example in the book. 
example <- function(num, noisemu=0, noisestd=0) {
  # Define the two circles. 
  N <- 5
  allX1 <- vector("list",N)
  allX2 <- vector("list",N)
  L1 <- rep(0, N)
  L2 <- rep(1, N)
  filler <- lapply(c(1:N), function(i) {
    print(i)
    C1 <- circle1(50, noisemu, noisestd)
    C2 <- circle2(50, noisemu, noisestd)
    Xlim <- c(-1.5, 1.5)
    Ylim <- c(-1.5, 1.5)
    # Use the same grid for both
    grid <- makeGrid(Xlim, Ylim, 0.065)
    X1 <- gridDiag(X=C1, FUN=kde, h=0.3, lim=cbind(Xlim,Ylim), by=0.065, sublevel=FALSE, library="Dionysus", printProgress=FALSE)
    allX1[[i]] <<- X1$diagram
    X2 <- gridDiag(X=C2, FUN=kde, h=0.3, lim=cbind(Xlim,Ylim), by=0.065, sublevel=FALSE, library="Dionysus", printProgress=FALSE)
    allX2[[i]] <<- X2$diagram
  })
  # Then we can combine these. 
  X <- cbind(allX1, allX2)
  L <- cbind(L1, L2)
  proba <- nhst(X, L, num)
  return(proba)
}



