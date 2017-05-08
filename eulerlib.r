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
  foam_storage <- array(0, dim=c(setnum, colnum, euler_length))

  for (i in 1:colnum) {  # loop over each col
    baseline_euler <- eulerVanillaFunction(baseline[[i]])
    baseline_storage[i,] <- baseline_euler
    for (t in 1:setnum) {  # loop over each param
      foam_euler <- eulerVanillaFunction(foam[[t]][[i]])
      foam_storage[t,i,] <- foam_euler
    }
  }
  foam_reshape <- array(
    foam_storage[1:dim(foam_storage)[1]-1,,],
    c(setnum-1, colnum, euler_length)
  )
  return(list(
    "baseline"=baseline_storage,
    "foam"=foam_reshape
  ))
}

# apply euler_tests to an entire directory
euler_voronoi_directory_test <- function(directory, num=100, norm=FALSE) {
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
      foam_storage <- abind(foam_storage, response[["foam"]], along=4)
      baseline_storage <- abind(baseline_storage, response[["baseline"]], along=3)
    }
  }
  foam_storage <- aperm(foam_storage, c(1, 4, 2, 3))
  baseline_storage <- aperm(baseline_storage, c(3, 1, 2))
  return(list(
    "baseline"=baseline_storage,
    "foam"=foam_storage
  ))
}

euler_eagle_directory_test <- function(directory, num_split=4, norm=FALSE) {
  if (substr(directory, nchar(directory), nchar(directory)+1) != "/") {
    directory <- paste(directory, "/", sep="")
  }
  return_obj <- list()
  for (i in 1:num_split) {
    print(paste('Generating euler set for split: ', i))
    cdm <- readRDS(paste(directory, 'cdm_diags_', i, '.rds', sep=''))
    wdm <- readRDS(paste(directory, 'wdm_diags_', i, '.rds', sep=''))
    cdm_foam <- vector('list', 1)
    cdm_foam[[1]] <- cdm
    response <- euler_tests(cdm_foam, wdm, norm=norm)
    return_obj[[paste("wdm", "split", i, sep="_")]] <- response[["baseline"]]
    return_obj[[paste("cdm", "split", i, sep="_")]] <- response[["foam"]]
  }
  return(return_obj)
}

# plot the average EC function (point wise average) and a
# 1 sd band around the average against the same thing for the null model.
# @baselineMat: (# examples, euler length)
# @foamMat: (# dim, # examples, euler_length)
eulerComparePlot <- function(
  baseline, foam, nullLabel="null", testLabel="test",
  plotTitle="Euler plot", outFile="./euler.png"
) {
  if (dim(baseline)[3] != dim(foam)[3]) {
    stop("dimension (3) of foam and baseline must match")
  }
  euler_length <- dim(baseline)[3]
  avg_baseline <- array(0, dim=(euler_length))
  avg_foam <- array(0, dim=(euler_length))
  minus_1std_baseline <- array(0, dim=(euler_length))
  minus_1std_foam <- array(0, dim=(euler_length))
  plus_1std_baseline <- array(0, dim=(euler_length))
  plus_1std_foam <- array(0, dim=(euler_length))
  for (d in 1:euler_length) {
    avg_baseline[d] <- mean(baseline[,,d])
    avg_foam[d] <- mean(foam[,,d])
    if (length(avg_baseline[d]) > 1) {
      minus_1std_baseline[d] <- avg_baseline[d] - sd(baseline[,,d])
      plus_1std_baseline[d] <- avg_baseline[d] + sd(baseline[,,d])
    } else {
      minus_1std_baseline[d] <- 0
      plus_1std_baseline[d] <- 0
    }

    if (length(avg_foam[d]) > 1) {
      minus_1std_foam[d] <- avg_foam[d] - sd(foam[,,d])
      plus_1std_foam[d] <- avg_foam[d] + sd(foam[,,d])
    } else {
      minus_1std_foam[d] <- 0
      plus_1std_foam[d] <- 0
    }
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
