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

# ----------------------------------------------------------------------
# Operations/Tools for a single diagram.

# Applied to a single diag.
normalize <- function(diag) {
  births <- diag[,2]
  deaths <- diag[,3]
  minbirths <- min(births)
  maxbirths <- max(births)
  mindeaths <- min(deaths)
  maxdeaths <- max(deaths)
  normbirths <- (births - minbirths) / (maxbirths - minbirths)
  normdeaths <- (deaths - mindeaths) / (maxdeaths - mindeaths)
  diag[,2] <- normbirths
  diag[,3] <- normdeaths
  return(diag)
}

# Given a voronsoi diagram, remove the first one of the 0th homologies, because it is always going to be infinity and thereby distracting.
# This somehow messes up plotting.
clean <- function(diag) {
  infinity <- diag[1,]
  if (infinity[[1]] == 0) { diag <- diag[2:length(diag[,1]),] }
  return(diag)
}
