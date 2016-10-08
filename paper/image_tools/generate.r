# Figure 1 - Eagle simulation images
source("process_eagle.r")
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
 
