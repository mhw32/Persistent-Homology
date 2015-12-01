

# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
library(rhdf5)
library(rgl)
library(scatterplot3d)



#d1 = t(h5read("Output_Eagle_Volume.hdf5", "P1/ParticlePositions"))
d1 = t(h5read("Output_Eagle_Volume.hdf5", "P1/SubhaloPositions"))


index0 = sample(1:nrow(d1),10000,replace = FALSE)
scatterplot3d(d1[index0,], pch = 19, tick.marks = FALSE, xlab = "", ylab = "",zlab = "", cex.symbols = .5)
title(paste(length(index0), " of ", nrow(d1), " particles", sep = ""))


apply(d1,2,fivenum)
