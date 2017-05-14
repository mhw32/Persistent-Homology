#!/usr/bin/env Rscript
source("corrTestSuite.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop(paste(
    "At least 2 arguments must be supplied ",
    "(dataDir[string], outDir[string], norm[boolean]).\n",
    sep="",
  ), call.=FALSE)
} else {
    dataDir <- args[1]
    outDir <- args[2]
    norm <- as.integer(args[3])
}

if (substr(dataDir, nchar(dataDir), nchar(dataDir)+1) != "/") {
  dataDir <- paste(dataDir, "/", sep="")
}

if (substr(outDir, nchar(outDir), nchar(outDir)+1) != "/") {
  outDir <- paste(outDir, "/", sep="")
}

load(paste(dataDir, 'cdm_corr_norm(', norm, ').rds', sep=''))
load(paste(dataDir, 'wdm_corr_norm(', norm, ').rds', sep=''))
tlibs <- corrSimuTestSuite(cdm, wdm)
corrDimTest <- tlibs$`indiv-corr`
corrDimProba <- corrDimTest(cdm, wdm)
saveRDS(
  corrDimProba, 
  file=paste(outDir, 'eagle_proba_norm(', norm, ').rds', sep='')
)
