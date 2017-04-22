source("process_eagle.r")
source("testlib.r")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop(paste(
    "At least 2 argument must be supplied ",
    "(outDir[string], res(olution)[float]).\n",
    sep="",
  ), call.=FALSE)
} else {
  outDir <- args[1]
  res <- as.double(args[2])
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

for (i in 2:4) {
  cdm <- load_CDM()
  wdm <- load_WDM()
  # Cut up the data
  cdm_slices <- slice_cube_robust(cdm, i)
  wdm_slices <- slice_cube_robust(wdm, i)
  # Get persistence diagrams
  cdm_diags <- persistify_set(cdm_slices, i, res=res)
  wdm_diags <- persistify_set(wdm_slices, i, res=res)
  # Save the data
  saveRDS(cdm_diags, paste(outDir, "cdm_diags_", i, ".rds", sep=""))
  saveRDS(wdm_diags, paste(outDir, "wdm_diags_", i, ".rds", sep=""))
}

