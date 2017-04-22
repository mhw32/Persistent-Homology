#!/usr/bin/env Rscript
source('testlib.r')

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop(paste(
    "At least 3 arguments must be supplied ",
    "(numReps[int], dataDir[string], outDir[string]).\n",
    sep="",
  ), call.=FALSE)
} else {
    numReps <- args[1]
    dataDir <- args[2]
    outDir <- args[3]
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

if (substr(dataDir, nchar(dataDir), nchar(dataDir)+1) != "/") {
  dataDir <- paste(dataDir, "/", sep="")
}

for (i in 1:numReps) {
  print(paste('Executing voronoi set [vanilla]: ', i))
  foam <- readRDS(paste(dataDir, "/foam", i, ".rds", sep=""))
  baseline <- readRDS(paste(dataDir, "/baseline", i, ".rds", sep=""))
  test_wrapper(
    foam, baseline, FALSE, FALSE,
    ext=paste(paste('set=['i, ']', sep=''), 'norm=[0]', 'base=[0.1]', sep='_'),
    folder=outDir,
  )

  print(paste('Executing voronoi set [normalized]: ', i))
  foam <- readRDS(paste(dataDir, "/foam", i, ".rds", sep=""))
  baseline <- readRDS(paste(dataDir, "/baseline", i, ".rds", sep=""))
  test_wrapper(
    foam, baseline, TRUE, FALSE,
    ext=paste(paste('set=['i, ']', sep=''), 'norm=[1]', 'base=[0.1]', sep='_'),
    folder=outDir,
  )
}
