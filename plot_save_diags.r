source("tools.r")
library(ks)

savepath <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"

for (i in 2:4) {
  num <- i * i * i
  cdmobj <- readRDS(paste(savepath, "cdm", i, ".rds", sep=""))
  wdmobj <- readRDS(paste(savepath, "wdm", i, ".rds", sep=""))
  for (j in 1:num) {
    cdmdiag <- cdmobj[[j]]
    wdmdiag <- wdmobj[[j]]
    # plot respective persistence diagrams.
    png(filename=paste(savepath, "cdmplot(", i, ",", j, ")", ".png", sep=""))
    plot(cdmdiag)
    dev.off()
    png(filename=paste(savepath, "wdmplot(", i, ",", j, ")", ".png", sep=""))
    plot(wdmdiag)
    dev.off()
    # calculate the pval.
    pval <- ks::kde.test(cdmdiag, wdmdiag)$pvalue
    print(paste("pvalue (", i, "cuts, ", j, "th cut: ", pval, sep=""))
  }
}
