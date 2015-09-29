# Author: Mike Wu
# Based on http://arxiv.org/pdf/1310.7467v1.pdf
# Randomization-style Null Hypothesis Significance Tests for Persistence Diagrams

library(TDA)
source('/Users/grub/Desktop/Cisewski-Lab/homology.r')

# TODO: THIS NEEDS TO BE SPED UP.
# X is a vector of persistence diagrams. 
# L is a vector of labels. 
lossfunc <- function(X, L) {
  # Split the persistence diagrams into 2 groups.
  G <- list(X[L == 0], X[L == 1])
  n <- c(length(G[[1]]), length(G[[2]]))
  sum1 <- 0 # Initialize outer counter.
  for (m in 1:2) {
    P <- G[[m]]
    C <- 1 / (2 * n[m] * (n[m] - 1))
    sum2 <- 0 # Initialize a counter.
    for (i in 1:n[m]) {
      for (j in 1:n[m]) {
        sum2 <- sum2 + bottleneck(P[[i]], P[[j]], dimension=1)
        # sum2 <- sum2 + wasserstein(P[i]$diagram, P[j]$diagram, dimension=1, p=2)
      }
    }
    sum1 <- sum1 + C * sum2
  }
  return(sum1)
}

# Given n1 + n2 persistence diagrams with n1 + n2 labels. Note that n1, and the following n2 diagrams must be disjoint. This is meshed together as X with labels L. Given N (number of iterations).
nhst <- function(X, L, N) {
  Z <- 0 # Initialization
  loss_orig <- lossfunc(X, L)
  # Preserve some amount of order.
  order <- c(1:length(X))
  for (i in 1:N) {
    random <- sample(order)
    loss_new <- lossfunc(X[random], L[random])
    if (loss_new <= loss_orig) {
      Z <- Z + 1
    }
  }
  # Average the outputs.
  Z <- Z / N
  return(Z)
}
# By law of large numbers, the expectation of Z --> P(loss(L_new) <= loss(L_obs)) as n --> infinity.

# -----------------------------------------------------------------

circle1 <- function(N, mu, sig) {
  K <- circleUnif(N, r=1)
  K <- addRandomNoise(K, mu, sig)
  for (i in 1:20) {
    K2 <- circleUnif(N, r=1)
    K <- rbind(K, K2)
    K <- addRandomNoise(K, mu, sig)
  }
  return(K)
}

circle2 <- function(N, mu, sig) {
  L1 <- circleUnif(N/2, 3/5)
  L1[,1] <- L1[,1] - 2/5
  L2 <- circleUnif(N/2, 2/5)
  L2[,1] <- L2[,1] + 3/5
  L <- rbind(L1, L2)
  L <- addRandomNoise(L, mu, sig)
  for (i in 1:20) {
    L1b <- circleUnif(N/2, 3/5)
    L1b[,1] <- L1b[,1] - 2/5
    L2b <- circleUnif(N/2, 2/5)
    L2b[,1] <- L2b[,1] + 3/5
    Lb <- rbind(L1b, L2b)
    L <- rbind(L, Lb)
    L <- addRandomNoise(L, mu, sig)
  }
  return(L)
}

# For demonstration, let's just do the example in the book. 
example <- function() {
  # Define the two circles. 
  N <- 50
  allX1 <- vector("list",5)
  allX2 <- vector("list",5)
  L1 <- c(0, 0, 0, 0, 0)
  L2 <- c(1, 1, 1, 1, 1)
  for (i in 1:5) {
    C1 <- circle1(50, 0, 0)
    C2 <- circle2(50, 0, 0)
    Xlim <- c(-1.5, 1.5)
    Ylim <- c(-1.5, 1.5)
    # Use the same grid for both
    grid <- makeGrid(Xlim, Ylim, 0.065)
    X1[[i]] <- gridDiag(X=C1, FUN=kde, h=0.3, lim=cbind(Xlim,Ylim), by=0.065, sublevel=FALSE, library="Dionysus", printProgress=FALSE)
    allX1[i] = X1$diagram
    X2[[i]] <- gridDiag(X=C2, FUN=kde, h=0.3, lim=cbind(Xlim,Ylim), by=0.065, sublevel=FALSE, library="Dionysus", printProgress=FALSE)
    allX2[i] = X2$diagram
  }
  # Then we can combine these. 
  X <- cbind(allX1, allX2)
  L <- cbind(L1, L2)
  proba <- nhst(X, L, 50)
  return(proba)
}



