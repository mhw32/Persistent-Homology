# -- Figure 10 : 3 cdm4 samples
library(scatterplot3d)

data <- readRDS('intermediate/fig_10_cdm4.rds')

nums <- c(51, 17, 34)
for (slice in nums) {
	pdf(paste("figure_10_cdm_slice_", slice, ".pdf", sep=""))
	scatterplot3d(data[[slice]], 
	              xlab='X Axis', 
	              ylab='Y Axis', 
	              zlab='Z Axis', 
	              pch='.',
	              color=rgb(0, 0, 0, 0.01),
	              cex.axis=1.5,
	              cex.lab=2)
	dev.off()
}
