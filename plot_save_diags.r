source("tools.r")
library(ks)

savepath <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"
normFlag <- FALSE

# Do the big one first
print("Everything was generated now with the SAME stats. 1.5 resolution, 10k particles")
cdmdiag <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/full_cdm_diag.rds")
wdmdiag <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/full_wdm_diag.rds")

pval <- ks::kde.test(cdm_diag$diagram, wdm_diag$diagram)$pvalue
print(paste("GKD Symmtric WDM-CDM Comparison pvalue (all): ", pval, sep=" "))
print("---------------------------------------------------------------------------")

# Do all the slices now
for (i in 2:4) {
  num <- i * i * i
  cdmobj <- readRDS(paste(savepath, "cdm", i, ".rds", sep=""))
  wdmobj <- readRDS(paste(savepath, "wdm", i, ".rds", sep=""))
  for (j in 1:num) {
    cdmdiag <- cdmobj[[j]]
    wdmdiag <- wdmobj[[j]]

    if (normFlag == TRUE) {
      cdmdiag <- normalize(cdmdiag)
      wdmdiag <- normalize(wdmdiag)
    }

    # plot respective persistence diagrams.
    png(filename=paste(savepath, "cdmplot(norm=", normFlag, ",", i, ",", j, ")", ".png", sep=""))
    plotDiag(cdmdiag)
    dev.off()
    png(filename=paste(savepath, "wdmplot(norm=", normFlag, ",", i, ",", j, ")", ".png", sep=""))
    plotDiag(wdmdiag)
    dev.off()

    # calculate the pval.
    pval <- ks::kde.test(cdmdiag, wdmdiag)$pvalue
    print(paste("GKD Symmetric WDM-CDM Comparison pvalue (", i, "cuts, ", j, "th cut): ", pval, sep=""))
  }
}
