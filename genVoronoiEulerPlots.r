#!/usr/bin/env Rscript
source('eulerlib.r')

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop(paste(
    "At least 4 arguments must be supplied ",
    "(dataDir[string], numFiles[integer], norm[boolean], plotDir[string].\n",
    sep="",
  ), call=FALSE)
} else {
  dataDir <- args[1]
  numFiles <- as.integer(args[2])
  norm <- as.logical(args[3])
  plotDir <- args[4]
}

if (substr(dataDir, nchar(dataDir), nchar(dataDir)+1) != "/") {
  dataDir <- paste(dataDir, "/", sep="")
}

if (substr(plotDir, nchar(plotDir), nchar(plotDir)+1) != "/") {
  plotDir <- paste(plotDir, "/", sep="")
}

response <- euler_voronoi_directory_test(
  directory=dataDir,
  num=numFiles,
  norm=norm
)

percFils <- seq(from=0.1, to=0.3, by=0.05)
baselineMat <- response[["baseline"]]
foamMat <- response[["foam"]]

if (length(percFils) != dim(foamMat)[1]) {
  stop('percFil dimension must match foam dimension')
}

# save the averages
stopifnot(dim(baselineMat)[2] == dim(foamMat)[3])
baselineAvg <- apply(baselineMat, c(1, 3), function(x) mean(x))
foamAvg <- apply(foamMat, c(1, 2, 4), function(x) mean(x))

saveRDS(
  baselineAvg,
  paste(
    plotDir,
    paste(
      "euler", "baseline", "averages",
      paste("[norm=", norm, "]", sep=""),
      sep="_"
    ),
    ".rds",
    sep=""
  )
)

saveRDS(
  foamAvg,
  paste(
    plotDir,
    paste(
      "euler", "foam", "averages",
      paste("[norm=", norm, "]", sep=""),
      sep="_"
    ),
    ".rds",
    sep=""
  )
)

for (i in 1:length(percFils)) {
  eulerComparePlot(
    baselineMat, foamMat[i,,,], nullLabel="Null [percfil=0.1]",
    testLabel=paste("Test [percfil=", percFils[i], "]", sep=""),
    plotTitle=paste("Null vs Test Hypothesis Euler Characteristic"),
    outFile=paste(
      plotDir,
      paste(
        "euler", "voronoi", "plot",
        paste("[percfil=", percFils[i], "]", sep=""),
        paste("[norm=", norm, "]", sep=""),
        sep="_"
      ),
      ".png",
      sep=""
    )
  )
}

