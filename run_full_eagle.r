source("process_eagle.r")
source("testlib.r")
source("localtest.r")
source("tools.r")
library(ks)

# Loading the data
cdm <- load_CDM()
wdm <- load_WDM()

# Plot these actual diagrams
# png(filename="./saved_states/cdmplot.png")
# scatterplot3d(cdm, pch=19, cex.symbols=0.001)
# dev.off()
# png(filename="./saved_states/wdmplot.png")
# scatterplot3d(wdm, pch=19, cex.symbols=0.001)
# dev.off()

# Settings for generating diagrams
res <- 1.5
boxlim <- c(0, 100) # Full size for the box.

# Create the one persistence diagram
cdm_diag <- gridDiag(cdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
wdm_diag <- gridDiag(wdm, dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)

# Save these diagrams.
saveRDS(cdm_diag, "./saved_states/full_cdm_diag.rds")
saveRDS(wdm_diag, "./saved_states/full_wdm_diag.rds")

# Run the single test we need  -- global kde test.
pval <- ks::kde.test(cdm_diag$diagram, wdm_diag$diagram)$pvalue
print(paste("Full global kde test result:",pval,sep=" "))

# Full global kde test result: 0.6944359
# Run the local kde test.
# localdiagplot(cdm_diag$diagram, wdm_diag$diagram, 0, "CDM/WDM 0 dim", "./saved_states/localks-0-dim.png")
# localdiagplot(cdm_diag$diagram, wdm_diag$diagram, 1, "CDM/WDM 1 dim", "./saved_states/localks-1-dim.png")
# localdiagplot(cdm_diag$diagram, wdm_diag$diagram, 2, "CDM/WDM 2 dim", "./saved_states/localks-2-dim.png")
