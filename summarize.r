# Fxns related to landscapes.
library(TDA)

# Landscape / Silhouette Functions
# ---------------------------------------------------------------
# Get area under the curve.
landscapeAUC <- function(diagram, KK=1, dim=1) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=1000)
  land <- landscape(diagram, KK=KK, dimension=dim, tseq)
  return(integrate(tseq, land))
}

# Get area under the curve.
silhouetteAUC <- function(diagram, p=1, dim=1) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=1000)
  silh <- silhouette(diagram, p=p, dimension=dim, tseq)
  return(integrate(tseq, silh))
}

sileuler <- function(diagram, p=1, length=1000) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=length)

  # Calculate silhouette for each dimension.
  s0 <- silhouette(diagram, p=p, dimension=0, tseq)
  s1 <- silhouette(diagram, p=p, dimension=1, tseq)
  s2 <- silhouette(diagram, p=p, dimension=2, tseq)

  # Calculate the alternating 'betti'.
  seuler <- rep(0, length)
  for (i in seq(length)) {
    seuler(i) = s0(i) - s1(i) + s2(i)
  }

  # Integrate the absolute value.
  score <- integrate(tseq, abs(seuler))
  return(score)
}

# This wrapper function returns a fxnized single dimension.
dimWrapper <- function(dim, fxn) {
  inner <- function(diagram) { return(fxn(diagram, dim=dim)) }
  return(inner)
}

# This wrapper function loops through all dimensions.
allWrapper <- function(fxn) {
  inner <- function(diagram, mindim=0, maxdim=2) {
    areas <- sapply(seq(from=mindim, to=maxdim, by=1), function(i) {
      fxn(diagram, dim=i)
    })
    return(areas)
  }
  return(inner)
}








