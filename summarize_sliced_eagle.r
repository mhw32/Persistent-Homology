source("process_eagle.r")
source("euler.r")
source("summarize.r")

savepath <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"

for (i in 2:4) {
  num <- i * i * i
  cdmobj <- readRDS(paste(savepath, "cdm", i, ".rds", sep=""))
  wdmobj <- readRDS(paste(savepath, "wdm", i, ".rds", sep=""))
  for (j in 1:num) {
    cdmdiag <- cdmobj[[j]]
    wdmdiag <- wdmobj[[j]]

    eulerPlot(cdmdiag, main=paste('CDM Slice (', j, '/', num, ') Euler Characteristic'), path=savepath)
    eulerPlot(wdmdiag, main=paste('WDM Slice (', j, '/', num, ') Euler Characteristic'), path=savepath)

    for (k in 1:3) {
      silhouettePlot(cdmdiag, dim=k, main=paste('CDM Slice (', j, '/', num, ')  Silhouette Dimension', k), path=savepath)
      silhouettePlot(wdmdiag, dim=k, main=paste('WDM Slice (', j, '/', num, ')  Silhouette Dimension', k), path=savepath)
    }
  }   
}