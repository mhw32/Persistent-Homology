source('euler.r')
source('tools.r')
library(abind)

euler_tests <- function(foam, baseline, norm=FALSE, run_base=TRUE) {
  foam[[length(foam)+1]] = baseline
  # normalize foam if needed
  if (norm == TRUE) { foam <- normFoam(foam) }
  # no need to remove vestigial for Euler tests
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  basenum <- setnum
  euler_length <- 1000
  # if we don't need to compare base with base
  if (run_base == FALSE) { setnum <- setnum - 1 }

  baseline_storage <- array(0, dim=c(colnum, euler_length))
  foam_storage <- array(0,dim=c(setnum, colnum, euler_length))

  for (i in 1:colnum) {  # loop over each col
    baseline_euler <- eulerVanillaFunction(baseline[[i]])
    baseline_storage[i,] <- baseline_euler
    for (t in 1:setnum) {  # loop over each param
      foam_euler <- eulerVanillaFunction(foam[[t]][[i]])
      foam_storage[t,i,] <- foam_euler
    }
  }

  return(list(
    "baseline"=baseline_storage,
    "foam"=foam_storage[1:dim(foam_storage)[1]-1,,]
  ))
}

# apply euler_tests to an entire directory
euler_directory_test <- function(directory, num=100, norm=FALSE) {
  if (substr(directory, nchar(directory), nchar(directory)+1) != "/") {
    directory <- paste(directory, "/", sep="")
  }
  for (i in 1:num) {
    print(paste('Generating euler set: ', i))
    foam <- readRDS(paste(directory, 'foam', i, '.rds', sep=''))
    baseline <- readRDS(paste(directory, 'baseline', i, '.rds', sep=''))
    response <- euler_tests(foam, baseline, norm=norm)
    if (i == 1) {
      foam_storage <- response[["foam"]]
      baseline_storage <- response[["baseline"]]
    } else {
      foam_storage <- abind(foam_storage, response[["foam"]], along=2)
      baseline_storage <- abind(baseline_storage, response[["baseline"]], along=1)
    }
  }
  return(list(
    "baseline"=baseline_storage,
    "foam"=foam_storage
  ))
}

# plot the average EC function (point wise average) and a
# 1 sd band around the average against the same thing for the null model.
# @baselineMat: (# examples, euler length)
# @foamMat: (# dim, # examples, euler_length)
eulerComparePlot <- function(
  baseline, foam, nullLabel="null", testLabel="test",
  plotTitle="Euler plot", outFile="./euler.png"
) {
  if (dim(baseline)[2] != dim(foam)[2]) {
    stop("dimension (2) of foam and baseline must match")
  }
  euler_length <- dim(baseline)[2]
  avg_baseline <- array(0, dim=(euler_length))
  avg_foam <- array(0, dim=(euler_length))
  minus_1std_baseline <- array(0, dim=(euler_length))
  minus_1std_foam <- array(0, dim=(euler_length))
  plus_1std_baseline <- array(0, dim=(euler_length))
  plus_1std_foam <- array(0, dim=(euler_length))
  for (d in 1:euler_length) {
    avg_baseline[d] <- mean(baseline[,d])
    avg_foam[d] <- mean(foam[,d])
    minus_1std_baseline[d] <- avg_baseline[d] - sd(baseline[,d])
    minus_1std_foam[d] <- avg_foam[d] - sd(foam[,d])
    plus_1std_baseline[d] <- avg_baseline[d] + sd(baseline[,d])
    plus_1std_foam[d] <- avg_foam[d] + sd(foam[,d])
  }
  plot_ymin <- min(min(minus_1std_baseline), min(minus_1std_foam))
  plot_ymax <- max(max(plus_1std_baseline), max(plus_1std_foam))
  lines_lwd <- 4.0
  euler_x <- seq(euler_length)
  png(filename=outFile)
  plot(
    euler_x, avg_baseline, col="orangered1", lwd=lines_lwd,
    type="l", xlab="index", ylab="euler", main=plotTitle,
    ylim=c(plot_ymin, plot_ymax)
  )
  lines(euler_x, avg_foam, col="springgreen3", lwd=lines_lwd)
  legend(
    "topright", c(nullLabel, testLabel), lty=c(1,1),
    lwd=c(lines_lwd,lines_lwd), col=c("orangered1", "springgreen3")
  )
  polygon(
    c(euler_x, euler_x),
    c(minus_1std_baseline, plus_1std_baseline),
    col=adjustcolor("tomato", alpha.f=0.2), border=NA,
  )
  polygon(
    c(euler_x, euler_x),
    c(minus_1std_foam, plus_1std_foam),
    col=adjustcolor("seagreen2", alpha.f=0.2), border=NA
  )
  dev.off()
}
