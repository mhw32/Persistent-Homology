#!/usr/bin/env Rscript
source('VoronoiFoam.r')

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop(paste(
    "At least 5 arguments must be supplied ",
    "(numReps[integer], res[double], numPerGroup[integer], outDir[string], genDiags[bool]).\n",
    sep="",
  ), call=FALSE)
} else {
  numReps <- as.integer(args[1])
  res <- as.double(args[2])  # default 0.5
  numPerGroup <- as.integer(args[3])  # default 15
  outDir <- args[4]
  genDiags <- as.logical(args[5])
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

N <- 5000
Boxlim <- c(0, 20)
perturb <- 1
baseline <- 0.1

for (i in 1:numReps) {
  print(paste('Generating voronoi set: ', i))
  if (genDiags) {
    voronoi_compilation(
      N=N, Boxlim=Boxlim, res=res, perturb=perturb,
      groupN=numPerGroup, baseline=baseline, nameId=i, folder=outDir
    )
  } else {
      voronoi_only_compilation(
        N=N, Boxlim=Boxlim, res=res, perturb=perturb,
        groupN=numPerGroup, baseline=baseline, nameId=i, folder=outDir
      )
  }
}

