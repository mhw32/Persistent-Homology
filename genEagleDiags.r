source("process_eagle.r")
source("testlib.r")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop(paste(
    "At least 3 argument must be supplied ",
    "(outDir[string], res(olution)[float], threshold[float], downsample[boolean]).\n",
    sep="",
  ), call.=FALSE)
} else {
  outDir <- args[1]
  res <- as.double(args[2])
  threshold <- as.double(args[3])
  downsample <- as.logical(args[4])
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

for (i in 1:1) {
  if (downsample) {
    cdm <- load_downsampled_CDM(threshold)
    wdm <- load_downsampled_WDM(threshold)
    cdm_diags_name <- paste("cdm_diags_downsampled_", threshold, "_", sep="")
    wdm_diags_name <- paste("wdm_diags_downsampled_", threshold, "_", sep="")
  } else {
    cdm <- load_CDM()
    wdm <- load_WDM()
    cdm_diags_name <- "cdm_diags_"
    wdm_diags_name <- "wdm_diags_"
  }
  # Cut up the data
  cdm_slices <- slice_cube_robust(cdm, i)
  wdm_slices <- slice_cube_robust(wdm, i)
  # Get persistence diagrams
  cdm_diags <- persistify_set(cdm_slices, i, res=res)
  wdm_diags <- persistify_set(wdm_slices, i, res=res)
  # Save the data
  saveRDS(cdm_diags, paste(outDir, cdm_diags_name, i, ".rds", sep=""))
  saveRDS(wdm_diags, paste(outDir, wdm_diags_name, i, ".rds", sep=""))
}

