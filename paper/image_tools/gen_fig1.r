# -- Figure 1 : Eagle simulation images --
source("../../process_eagle.r")
library(scatterplot3d)
library(rhdf5)
library(rgl)

load_CDM <- function() {
  d <- t(h5read("../../simulations/Output_Eagle_Volume.hdf5", "P1/SubhaloPositions"))
  return(d)
}

# Increase a processing step to ensure only genuine objects.
load_WDM <- function() {
  d <- t(h5read("../../simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloPositions"))
  HMspher <- t(h5read("../../simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloHMSpher"))
  MaxMass <- t(h5read("../../simulations/Output_Eagle_VolumeW.hdf5", "P1/SubhaloMaxMass"))
  selection <- (MaxMass > 2.2e8) & (HMspher > 0.2)
  reframe <- d[selection,]
  return(reframe)
}


whole_cdm <- load_CDM()
whole_wdm <- load_WDM()

pdf("figure_1_whole_cdm.pdf")
scatterplot3d(whole_cdm, 
              xlab='X Axis', 
              ylab='Y Axis', 
              zlab='Z Axis', 
              pch='.',
              color=rgb(0, 0, 0, 0.01),
              cex.axis=1.5,
              cex.lab=2)
dev.off()
pdf("figure_1_whole_wdm.pdf")
scatterplot3d(whole_wdm,
              xlab='X Axis',
              ylab='Y Axis',
              zlab='Z Axis',
              pch='.',
              color=rgb(0, 0, 0, 0.01),
              cex.axis=1.5,
              cex.lab=2) 
dev.off()
 
