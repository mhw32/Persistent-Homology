source("process_eagle.r")
source("testlib.r")
source("localtest.r")
library(ks)

savepath <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"

for (i in 2:4) {
  # Load Data
  cdm <- load_CDM()
  wdm <- load_WDM()
  # Cut up the data
  cdm_slices <- slice_cube_robust(cdm, i)
  wdm_slices <- slice_cube_robust(wdm, i)
  # Get persistence diagrams
  cdm_diags <- persistify_set(cdm_slices, i)
  wdm_diags <- persistify_set(wdm_slices, i)
  # Save the data
  saveRDS(cdm_diags, paste(savepath, "cdm", i, ".rds", sep=""))
  saveRDS(wdm_diags, paste(savepath, "wdm", i, ".rds", sep=""))
}

