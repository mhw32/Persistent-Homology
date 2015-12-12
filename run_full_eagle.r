source("process_eagle.r")
source("testlib.r")
source("localtest.r")
library(ks)
# Loading the data
cdm <- load_CDM()
wdm <- load_WDM()
# Plot these actual diagrams
png(filename="./saved_states/cdmplot.png")
scatterplot3d(cdm, pch=19, cex.symbols=0.005)
dev.off()
png(filename="./saved_states/wdmplot.png")
scatterplot3d(cdm, pch=19, cex.symbols=0.005)
dev.off()
# Settings for generating diagrams
res <- 2
boxlim <- c(0, 100) # Full size for the box.
# Create the one persistence diagram
cdm_diag <- gridDiag(cdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
wdm_diag <- gridDiag(wdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
# Save these diagrams.
saveRDS(cdm_diag, "./saved_states/full_cdm_diag.rds")
saveRDS(wdm_diag, "./saved_states/full_wdm_diag.rds")
# Run the single test we need  -- global kde test.
pval <- ks::kde.test(cdm_diag, wdm_diag)$pvalue
print(paste("Full global kde test result:",pval,sep=" "))
# Run the local kde test.
localdiagplot(cdm_diag, wdm_diag, 0, "CDM/WDM 0 dim", "./saved_states/localks-0-dim.pdf")
localdiagplot(cdm_diag, wdm_diag, 1, "CDM/WDM 1 dim", "./saved_states/localks-1-dim.pdf")
localdiagplot(cdm_diag, wdm_diag, 2, "CDM/WDM 2 dim", "./saved_states/localks-2-dim.pdf")
