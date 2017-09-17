source("process_eagle.r")
source("distance.r")
source("tools.r")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop(paste(
    "At least 2 arguments must be supplied ",
    "(norm[boolean], dataDir[string], outDir[string]).\n",
    sep="",
  ), call.=FALSE)
} else {
    norm <- as.logical(args[1])
    dataDir <- args[2]
    outDir <- args[3]
}

if (substr(dataDir, nchar(dataDir), nchar(dataDir)+1) != "/") {
  dataDir <- paste(dataDir, "/", sep="")
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

dists <- array(0,dim=c(4, 64, 3))

for (i in 1:4) {
  print(paste('Executing distance calculation set: ', i))
  cdm_slices <- readRDS(paste(dataDir, '/cdm_diags_', i, '.rds', sep=''))
  wdm_slices <- readRDS(paste(dataDir, '/wdm_diags_', i, '.rds', sep=''))

  for (l in 1:i^3) {
    cdm_diag <- cleanDiag(cdm_slices[[l]])
    wdm_diag <- cleanDiag(wdm_slices[[l]])

    if (norm) { # optionally normalize
      cdm_diag <- normalize(cdm_diag)
      wdm_diag <- normalize(wdm_diag)
    }

    # get the bottle neck distances.
    dists[i,l,1] <- bottleneckDist(cdm_diag, wdm_diag, dimension=0)
    dists[i,l,2] <- bottleneckDist(cdm_diag, wdm_diag, dimension=1)
    dists[i,l,3] <- bottleneckDist(cdm_diag, wdm_diag, dimension=2)
  }
}

if (norm) {
  out_path <- paste(outDir, '/eagle_dists_norm.rds', sep='');
} else {
  out_path <- paste(outDir, '/eagle_dists.rds', sep='');
}
saveRDS(dists, file=out_path)
