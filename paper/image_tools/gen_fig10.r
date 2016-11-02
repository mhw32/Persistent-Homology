# -- Figure 10 : 3 cdm4 samples
library(scatterplot3d)
# library(rhdf5)
# library(rgl)

# load_CDM <- function() {
#   d <- t(h5read("../../simulations/Output_Eagle_Volume.hdf5", "P1/SubhaloPositions"))
#   return(d)
# }

# slice_cube_robust <- function(cube, n) {
#   # This is the final slices.
#   slices <- vector('list', n^3)
#   slicenum <- 1
#   # initialize some item/counters
#   item <- 0
#   counter <- 1
#   # create an array of splits
#   h <- rep(NA, n+1)
#   while (item < 99) {
#     h[counter] <- item
#     item <- item + 100/n
#     counter <- counter + 1
#   }
#   h[counter] <- 100
#   # Loop through each dimension
#   for (i in seq(2,n+1)) {
#     for (j in seq(2,n+1)) {
#       for (k in seq(2,n+1)) {
#         dim1 <- c(h[i-1], h[i])
#         dim2 <- c(h[j-1], h[j])
#         dim3 <- c(h[k-1], h[k])
#         # Get a boolean array of each splice.
#         logic <- (cube[,1] > dim1[1] & cube[,1] <= dim1[2]) & (cube[,2] > dim2[1] & cube[,2] <= dim2[2]) & (cube[,3] > dim3[1] & cube[,3] <= dim3[2])
#         newcube <- cube[logic,]
#         # Renormalize everything (0, 100/n)
#         newcube[,1] <- newcube[,1] - (h[i-1])
#         newcube[,2] <- newcube[,2] - (h[j-1])
#         newcube[,3] <- newcube[,3] - (h[k-1])
#         slices[[slicenum]] <- newcube
#         slicenum <- slicenum + 1
#       }
#     }
#   }
#   return(slices)
# }

# cdm <- load_CDM()
# data <- slice_cube_robust(cdm, 4)
# saveRDS(data, 'intermediate/fig_10_cdm4.rds')

data <- readRDS('intermediate/fig_10_cdm4.rds')

nums <- c(51, 17, 34)
for (slice in nums) {
	pdf(paste("figure_10_cdm_slice_", slice, ".pdf", sep=""))
	scatterplot3d(data[[slice]], 
	              xlab='', 
	              ylab='', 
	              zlab='', 
	              pch='.',
	              color=rgb(0, 0, 0, 0.1),
	              tick.marks=FALSE,
              	label.tick.marks=FALSE)
	dev.off()
}
