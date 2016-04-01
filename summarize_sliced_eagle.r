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

    eulerDualPlot(cdmdiag, wdmdiag, main1=paste('CDM Slice (', j, '/', num, ') Euler Characteristic'), main2=paste('WDM Slice (', j, '/', num, ') Euler Characteristic'), path=paste(savepath, "wdm_cdm_euler_", j, "_out_of_", num, ".png", sep=""))

    for (k in 1:3) {
      silhouetteDualPlot(cdmdiag, wdmdiag, dim=k, main1=paste('CDM Slice (', j, '/', num, ')  Silhouette Dimension', k), main2=paste('WDM Slice (', j, '/', num, ')  Silhouette Dimension', k), path=paste(savepath, "cdm_silh_", j, "_out_of_", num, "_dim_", k, ".png", sep=""))
    }
  }   
}
