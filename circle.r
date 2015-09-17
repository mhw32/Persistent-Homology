library(TDA)
# Example with RipsDiag function
XX <- circleUnif(30)  # This creates a circular object (aka point cloud)
# Second arg = max dimension of the homological features (like how high do I want to go in dimension).
# Third arg = max value for rips filtration (What is this?)
Diag <- ripsDiag(XX, 5, 10, printProgress = TRUE)
Plotting <- Diag$diagram
# Plot the actual thing
plot(Plotting)
# You can also plot the barcode
# plot(Diag[["diagram"]], barcode = TRUE, main = "Barcode")

# A second more random Example
## distance matrix for triangle with edges of length: 1,2,4
distX <- matrix(c(0, 1, 2, 1, 0, 4, 2, 4, 0), ncol = 3)
maxscale <- 5
maxdimension <- 1
## note that the input distXX is a distance matrix
DiagTri <- ripsDiag(distX, maxdimension, maxscale, dist = "arbitrary", printProgress = TRUE)
Plotting <- DiagTri$diagram
plot(Plotting)

# Example 3
Circle1 <- circleUnif(60)
Circle2 <- circleUnif(60, r = 2) + 3
Circles <- rbind(Circle1, Circle2)
plot(Circles)
maxscale <- 5
maxdimension <- 1 # (0 for components, 1 for loops, 2 for voids, etc.):





