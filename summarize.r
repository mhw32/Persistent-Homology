# Fxns related to landscapes.
library(TDA)

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

# Wrapper function to do area for all dim
dimWrapper <- function(fxn) {
  calculate <- function(diagram, mindim=0, maxdim=2) {
    areas <- sapply(seq(from=mindim, to=maxdim, by=1), function(i) {
      fxn(diagram, dim=i)
    })
    return(areas)
  }
  return(calculate)
}

