source('testlib.r')
source('process_eagle.r')

# Unnormalized with wdm baseline.
cdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/cdm_diags_2.rds")
wdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/wdm_diags_2.rds")
bigset <- vector('list', 1)
bigset[[1]] <- cdm
test_wrapper(bigset, wdm, "split2-dmtest-wdmbase", FALSE)
# Unnormalized with cdm baseline.
cdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/cdm_diags_2.rds")
wdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/wdm_diags_2.rds")
bigset <- vector('list', 1)
bigset[[1]] <- wdm
test_wrapper(bigset, cdm, "split2-dmtest-cdmbase", FALSE)
# Normalized with wdm baseline.
cdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/cdm_diags_2.rds")
wdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/wdm_diags_2.rds")
bigset <- vector('list', 1)
bigset[[1]] <- cdm
test_wrapper(bigset, wdm, "split2-dmtest-norm-wdmbase", TRUE)
# Normalized with cdm baseline.
cdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/cdm_diags_2.rds")
wdm <- readRDS("/home/fas/cisewski/mhw32/scratch/homology/saved_states/wdm_diags_2.rds")
bigset <- vector('list', 1)
bigset[[1]] <- wdm
test_wrapper(bigset, cdm, "split2-dmtest-norm-cdmbase", TRUE)


