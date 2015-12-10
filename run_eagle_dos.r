source("process_eagle.r")

cdm <- load_CDM()
wdm <- load_WDM()

cdm_slices <- slice_cube(cdm, 4)
wdm_slices <- slice_cube(wdm, 4)

cdm_diags <- persistify_set(cdm_slices)
wdm_diags <- persistify_set(wdm_slices)

saveRDS(cdm_diags, "/home/fas/cisewski/mhw32/scratch/homology/saved_states/cdm_diags.rds")
saveRDS(wdm_diags, "/home/fas/cisewski/mhw32/scratch/homology/saved_states/wdm_diags.rds")
