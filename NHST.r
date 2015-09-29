# Author: Mike Wu
# Based on http://arxiv.org/pdf/1310.7467v1.pdf
# Randomization-style Null Hypothesis Significance Tests for Persistence Diagrams

library(TDA)

# Define a distance function.
# X1, X2 are two persistence objects. 
# Dim is the order of homological object.
# p is the power needed for wassierstein. It should be set if using wasserstein distance.
distance <- function(P1, P2, dim=1, type='bottleneck', p=NULL) {
  if !(type %in% c('bottleneck', 'wasserstein')) type <- 'botteneck'
  # Return the bottleneck distance.
  if (type == 'botteneck') {
    result <- bottleneck(P1$diagram, P2$diagram, dimension=dim)
  } else { # Return the Wasserstein distance.
    p <- if (is.null(p)) 1 else p
    result <- wasserstein(P1$diagram, P2$diagram, p=p, dimension=dim)
  }
  return(distance)
}

# TODO: THIS NEEDS TO BE SPED UP.
# X is a vector of persistence diagrams. 
# L is a vector of labels. 
lossfunc <- function(X, L, dim=1, type='bottleneck', p=NULL) {
  # Split the persistence diagrams into 2 groups.
  G <- c(X[L == 0], X[L == 1])
  n <- c(length(X1), length(X2))
  sum1 <- 0 # Initialize outer counter.
  for (m in 1:2) {
    P <- G[m]
    C <- 1 / (2 * n[m] * (n[m] - 1))
    sum2 <- 0 # Initialize a counter.
    for (i in 1:n[m]) {
      for (j in 1:n[m]) {
        sum2 <- sum2 + distance(P[i], P[j], dim, type, p)
      }
    }
    sum1 <- sum1 + C * sum2
  }
  return(sum1)
}

# Given n1 + n2 persistence diagrams with n1 + n2 labels. Note that n1, and the following n2 diagrams must be disjoint. This is meshed together as X with labels L. Given N (number of iterations).
nhst <- function(X, L, N, dim=1, type='bottleneck', p=NULL) {
  Z <- 0 # Initialization
  loss_orig <- lossfunc(X, L, dim, type, p)
  # Preserve some amount of order.
  order <- c(1:length(X))
  for (i in 1:N) {
    random <- sample(order)
    loss_new <- lossfunc(X[random], L[random], dim, type, p)
    if (loss_new <= loss_orig) {
      Z <- Z + 1
    }
  }
  # Average the outputs.
  Z <- Z / N
  return(Z)
}
