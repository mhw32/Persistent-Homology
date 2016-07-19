# Collection of useful tools to share among files.

# Math functions
# Combination
choose <- function(n, k) {
  factorial(n) / (factorial(k)*factorial(n-k))
}
# Permutation
pick <- function(n, k) {
  factorial(n) / factorial(n-k)
}

hotelling_pval <- function(p) {
  pval <- pf(1/p$statistic[2],p$parameter[2], p$parameter[1])
  return(pval)
}

# A generic function for performing operations across sequences of foam structures.
gridOperation <- function(foam, fxn) {
  # Get dimension of foams.
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  # Do something to the first column.
  matrix <- sapply(seq(1:colnum), function(j) {
    fxn(foam[[1]][[j]])
  })
  maxdim <- if (is.null(dim(matrix))) 1 else length(dim(matrix))
  # Do something to the rest of the columns.
  for (i in 2:setnum) {
    vector <- sapply(seq(1:colnum), function(j) {
      fxn(foam[[i]][[j]])
    })
    matrix <- abind(matrix, vector, along=maxdim+1)
  }
  colnames(matrix) <- NULL
  return(matrix)
}

# Take only a single dimension out of a persistence diagram.
sliceDim <- function(diag, dim) {
  multi <- if (length(dim) > 1) TRUE else FALSE
  logic <- if (multi) diag[,1] %in% dim else diag[,1] == dim
  bd <- cbind(diag[,2][logic], diag[,3][logic])
  colnames(bd) <- c('Birth', 'Death')
  return(bd)
}

# Remove the first 0th dimension component from every single nested persistence diagram in the foam object (nested list).
cleanFoam <- function(foam) {
  # Cleaning a single persistence diagram.
  cleanDiag <- function(diag) { diag[2:nrow(diag),] }
  # Applying cleanDiag to an entire vector.
  cleanVec <- function(vec) {
    vector <- lapply(vec, function(diag) { cleanDiag(diag) })
    return(vector)
  }
  newfoam <- lapply(foam, function(vec) { cleanVec(vec) })
  return(newfoam)
}

# Prior to cleaning
normFoam <- function(foam) {
  normVec <- function(vec) {
    vector <- lapply(vec, function(diag) { normalize(diag) })
    return(vector)
  }
  newfoam <- lapply(foam, function(vec) { normVec(vec) })
  return(newfoam)
}

# ----------------------------------------------------------------------
# Operations/Tools for a single diagram.
cleanDiag <- function(diag) { diag[2:nrow(diag),] }

# Applied to a single diag.
normalize <- function(diag) {
  normdiag <- cleanDiag(diag)

  births <- diag[,2]
  deaths <- diag[,3]
  normbirths <- normdiag[,2]
  normdeaths <- normdiag[,3]

  minbirths <- min(normbirths)
  maxbirths <- max(normbirths)
  mindeaths <- min(normdeaths)
  maxdeaths <- max(normdeaths)

  maxtotal <- max(maxbirths, maxdeaths)
  mintotal <- min(minbirths, mindeaths)

  newbirths <- (births - mintotal) / (maxtotal - mintotal)
  newdeaths <- (deaths - mintotal) / (maxtotal - mintotal)
  diag[,2] <- newbirths
  diag[,3] <- newdeaths
  return(diag)
}

# Often it is easier to make my own plotting function to allow editing.
plotDiag <- function(X){
  # X is the persistence diagram
  plot(X[,2], X[,3], pch = c(X[,1]+1), col = c(X[,1]+1), xlab = "Birth", ylab = "Death", main = "", cex.lab=2.0, cex.axis=2.0, cex.main=2.0, cex.sub=2.0)
  par(mar=c(1.1,1.1,1.1,1.1))
  abline(a = 0, b = 1)
  legend("topleft", c("0","1","2"), pch = c(1,2,3), col = c(1,2,3), cex=2.0)
}

