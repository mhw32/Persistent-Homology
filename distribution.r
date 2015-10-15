# Direct Distribution Functions.
# Mike Wu

# Method: Parametrize the data. Because I can't really find a good multivariate non-parametric test, I can try to project the 2-dimension data into 1-dimension. I want to preserve relative distances between points. We note that by construction, the points in the persistence diagram follow a diagonal line. Therefore, we can rotate the plot and put the diagonal in the x-axis. (birth + death) / 2 and (birth - death) / 2. 
# This is actually a matrix transformation. We also need to scale the standard deviation properly: A <- 1/sqrt(2)*matrix(c(1,1,1,-1),2,2). This should preserve all distances and essentially equals SVD. Then we can remove the 2nd dimension and be left with a scaled dimension. This is the first principle component!

source('tools.r')

# 2-to-1 Dimensionality reduction based on 
reduce <- function(mat) {
  A <- 1/sqrt(2)*matrix(c(1,1,1,-1),2,2)
  z <- mat %*% A # Reparametrization
  return(z[,1]) # z --> The diagonal (should represent all the info)
}

# Given a persistence diagram. Pick out a single dimension. Reduce dimension. Run a Manning Whitney U test (nonparametric 1-d).
distribDimStat <- function(set, dim) {
  setnum <- length(set)
  reducedset <- sapply(1:setnum, function(i) {
    # Reduce specific dimension(s).
    reduce(sliceDim(set[[i]], dim)) 
  })
  # Because these are all from a single high dimensional distribution, I will combine them together.
  vector <- unlist(reducedset)
  return(vector)
}


