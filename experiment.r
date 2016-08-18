# This is to generate the example for the paper. 
# Do 5 peaks at the beginning going down to 300 points. 

library(TDA)
source("tools.r")

loadBinFile <- function(path="~/Desktop/data.bin") {
  infile <- path
  con <- file(infile, "rb")
  dim <- readBin(con, "integer", 2)
  Mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
  close(con)
  return(Mat)
}

x1 <- circleUnif( n=3 , r=1 ) # z = 1 --> 0.7
# x2 <- circleUnif( n=10 , r=1 ) # z = 0.7 --> 0.4
# x3 <- circleUnif( n=50 , r=1 ) # z = 0.4 --> 0.2
x4 <- circleUnif( n=100 , r=1 ) # z = 0.2 --> 0.0

num_total <- 0
for (i in seq(0, 1, 0.001)) {
	if (i <= 0.5) {
		num_total <- num_total + 3
	} else {
		num_total <- num_total + 100
	}
}

count <- 1
data <- matrix(, ncol=2, nrow=num_total)
for (i in seq(0, 1, 0.001)) {
	if (i <= 0.5) {
		numdata <- 3
		data[count:(count+numdata-1),] <- x1
	} else {
		numdata <- 100
		data[count:(count+numdata-1),] <- x4
	}
	count <- count + numdata
}

## Ranges of the grid
Xlim <- c(-1.8, 1.8)
Ylim <- c(-1.8, 1.8)
lim <- cbind(Xlim, Ylim)
by <- 0.05
h <- .3 #bandwidth

diag <- gridDiag(data, kde, lim = lim, by = by, sublevel = FALSE, location = TRUE, printProgress = TRUE, h = h)
diag <- cleanDiag(diag$diagram)




