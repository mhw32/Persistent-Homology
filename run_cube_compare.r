source("process_eagle.r")
library(ks)

cdm_splits <- readRDS("./saved_states/cdm_diags.rds")
wdm_splits <- readRDS("./saved_states/wdm_diags.rds")

resArr <- c(1, 0.5, 0.25, 0.1)
boxlim <- c(0, 33)
len <- length(cdm_diags)

for (res in resArr) {
  for (l in 1:len) {
    cdm_diag <- gridDiag(cdm_splits[[l]], dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    wdm_diag <- gridDiag(wdm_splits[[l]], dtm, lim=cbind(boxlim, boxlim, boxlim), by=res, sublevel=T, printProgress=T, m0=0.001)
    pval <- ks::kde.test(cdm_diag$diagram, wdm_diag$diagram)$pvalue
    print(paste("res: ", res, " idx: ", l, " pval: ", pval, sep=" "))
  }
}
