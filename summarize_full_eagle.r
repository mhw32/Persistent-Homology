source("process_eagle.r")
source("euler.r")
source("summarize.r")

cdm_source <- "./saved_states/full_cdm_diag.rds"
wdm_source <- "./saved_states/full_wdm_diag.rds"

cdm_diag <- readRDS(cdm_source)
wdm_diag <- readRDS(wdm_source)

eulerPlot(cdm_diag, main='CDM Full Euler Characteristic', path='./saved_states/real_data/cut_diags/')
eulerPlot(wdm_diag, main='WDM Full Euler Characteristic', path='./saved_states/real_data/cut_diags/')

for (i in 1:3) {
  silhouettePlot(cdm_diag, dim=i, main=paste('CDM Full Silhouette Dimension', i), path='./saved_states/real_data/cut_diags/')
  silhouettePlot(wdm_diag, dim=i, main=paste('WDM Full Silhouette Dimension', i), path='./saved_states/real_data/cut_diags/')
}