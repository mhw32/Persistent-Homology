#!/usr/bin/env Rscript
source('eulerlib.r')

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop(paste(
    "At least 4 arguments must be supplied ",
    "(dataDir[string], norm[boolean], plotDir[string].\n",
    sep="",
  ), call=FALSE)
} else {
  dataDir <- args[1]
  norm <- as.logical(args[3])
  plotDir <- args[4]
}

if (substr(dataDir, nchar(dataDir), nchar(dataDir)+1) != "/") {
  dataDir <- paste(dataDir, "/", sep="")
}

if (substr(plotDir, nchar(plotDir), nchar(plotDir)+1) != "/") {
  plotDir <- paste(plotDir, "/", sep="")
}

response <- euler_eagle_directory_test(
  directory=dataDir,
  num_split=4,
  norm=norm
)

cdmAverages <- array(, dim=c(4, 1000))
wdmAverages <- array(, dim=c(4, 1000))
for (i in 1:4) {
    cdmSplit <- response[[paste("cdm", "split", i, sep="_")]]
    wdmSplit <- response[[paste("wdm", "split", i, sep="_")]]
    cdmAverage <- apply(cdmSplit, c(3), function(x) mean(x))
    wdmAverage <- apply(wdmSplit, c(2), function(x) mean(x))
    cdmAverages[i,] <- cdmAverage
    wdmAverages[i,] <- wdmAverage
}

saveRDS(
  cdmAverages,
  paste(
    plotDir,
    paste(
      "euler", "cdm", "averages",
      paste("[norm=", norm, "]", sep=""),
      sep="_"
    ),
    ".rds",
    sep=""
  )
)

saveRDS(
  wdmAverages,
  paste(
    plotDir,
    paste(
      "euler", "wdm", "averages",
      paste("[norm=", norm, "]", sep=""),
      sep="_"
    ),
    ".rds",
    sep=""
  )
)

for (i in 1:4) {
  cdmMat <- response[[paste("cdm", "split", i, sep="_")]][1,,]
  wdmMat <- response[[paste("wdm", "split", i, sep="_")]]
  eulerComparePlot(
    wdmMat, cdmMat,
    nullLabel=paste("Null [wdm (split ", i ")]", sep=""),
    testLabel=paste("Test [cdm (split ", i ")]", sep=""),
    plotTitle=paste("Null vs Test Hypothesis Euler Characteristic"),
    outFile=paste(
      plotDir,
      paste(
        "euler", "eagle", "plot",
        paste("[split=", i, "]", sep=""),
        paste("[norm=", norm, "]", sep=""),
        sep="_"
      ),
      ".png",
      sep=""
    )
  )
}

