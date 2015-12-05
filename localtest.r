library(ks)
library(MASS)
library(fields)
source('tools.r')

localsigplot <- function(x1, x2, main, path) {
  # Perform a local significance test.
  loct <- kde.local.test(x1, x2)
  png(filename=path)
  image.plot(loct$fhat1$eval.points[[1]], loct$fhat1$eval.points[[2]],loct$pvalue, main=main)
  dev.off()
}

localdiagplot <- function(diag1, diag2, dim, main, path) {
  diagdim1 <- cbind(diag1[,2][diag1[,1] == dim], diag1[,3][diag1[,1] == dim])
  diagdim2 <- cbind(diag2[,2][diag2[,1] == dim], diag2[,3][diag2[,1] == dim])
  localsigplot(diagdim1, diagdim2, main, path)
}
