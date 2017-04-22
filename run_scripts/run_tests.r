source('testlib.r')
source('process_eagle.r')

source <- "/home/fas/cisewski/mhw32/scratch/homology/saved_states/real_data/cut_diags/"

for (i in 2:4) {
  cdm <- readRDS(paste(source, "cdm", i, ".rds", sep=""))
  wdm <- readRDS(paste(source, "wdm", i, ".rds", sep=""))
  bigset <- vector('list', 1)
  bigset[[1]] <- cdm
  test_wrapper(bigset, wdm, paste("Split", i, "WDMbaseNormFalse", sep=""), norm=FALSE, run_base=FALSE)

  cdm <- readRDS(paste(source, "cdm", i, ".rds", sep=""))
  wdm <- readRDS(paste(source, "wdm", i, ".rds", sep=""))
  bigset <- vector('list', 1)
  bigset[[1]] <- cdm
  test_wrapper(bigset, wdm, paste("Split", i, "WDMbaseNormTrue", sep=""), norm=TRUE, run_base=FALSE)

  cdm <- readRDS(paste(source, "cdm", i, ".rds", sep=""))
  wdm <- readRDS(paste(source, "wdm", i, ".rds", sep=""))
  bigset <- vector('list', 1)
  bigset[[1]] <- wdm
  test_wrapper(bigset, cdm, paste("Split", i, "CDMbaseNormFalse", sep=""), norm=FALSE, run_base=FALSE)

  cdm <- readRDS(paste(source, "cdm", i, ".rds", sep=""))
  wdm <- readRDS(paste(source, "wdm", i, ".rds", sep=""))
  bigset <- vector('list', 1)
  bigset[[1]] <- wdm
  test_wrapper(bigset, cdm, paste("Split", i, "CDMbaseNormTrue", sep=""), norm=TRUE, run_base=FALSE)
}
