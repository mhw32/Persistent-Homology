sourcepath <- "/home/fas/cisewski/mhw32/scratch/homology/"
source(paste(sourcepath, "process_eagle.r"))
source(paste(sourcepath, "testlib.r"))
source(paste(sourcepath, "localtest.r"))
library(ks)

savepath <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"

for (i in 2:4) {
  # Load Data
  cdm <- load_CDM()
  wdm <- load_WDM()
  # Cut up the data
  cdm_slices <- slice_cube_robust(cdm, i)
  wdm_slices <- slice_cube_robust(wdm, i)
  # Save the data
  saveRDS(cdm_slices, paste(savepath, "cdm", i, ".rds"))
  saveRDS(wdm_slices, paste(savepath, "wdm", i, ".rds"))
}

