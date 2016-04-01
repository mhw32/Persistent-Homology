source("process_eagle.r")
source("euler.r")
source("summarize.r")

cdm_source <- "./saved_states/full_cdm_diag.rds"
wdm_source <- "./saved_states/full_wdm_diag.rds"

cdm_diag <- readRDS(cdm_source)$diagram
wdm_diag <- readRDS(wdm_source)$diagram

eulerPlot(cdm_diag, main='CDM Full Euler Characteristic', path='./saved_states/real_data/cut_diags/cdm_full_euler.png')
eulerPlot(wdm_diag, main='WDM Full Euler Characteristic', path='./saved_states/real_data/cut_diags/wdm_full_euler.png')

for (i in 1:3) {
  silhouettePlot(cdm_diag, dim=i, main=paste('CDM Full Silhouette Dimension', i), path=paste('./saved_states/real_data/cut_diags/cdm_full_silh_', i, '.png', sep=""))
  silhouettePlot(wdm_diag, dim=i, main=paste('WDM Full Silhouette Dimension', i), path=paste('./saved_states/real_data/cut_diags/wdm_full_silh_', i, '.png', sep=""))
}
