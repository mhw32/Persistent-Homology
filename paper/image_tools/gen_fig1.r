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

png("figure_1_whole_cdm.png")
scatterplot3d(whole_cdm, 
              xlab='', 
              ylab='', 
              zlab='', 
              pch='.',
              color=rgb(0, 0, 0, 0.01),
              tick.marks=FALSE,
              label.tick.marks=FALSE)
dev.off()
png("figure_1_whole_wdm.png")
scatterplot3d(whole_wdm,
              xlab='', 
              ylab='', 
              zlab='',
              pch='.',
              color=rgb(0, 0, 0, 0.01),
              tick.marks=FALSE,
              label.tick.marks=FALSE) 
dev.off()
 
