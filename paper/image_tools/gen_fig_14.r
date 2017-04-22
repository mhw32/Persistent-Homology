#!/usr/bin/r
library(scatterplot3d)

infile <- "intermediate/fig_14_min_data_top_3.bin"
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
minMat <- matrix(readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

infile <- "intermediate/fig_14_max_data_top_3.bin"
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
maxMat <- matrix(readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

cdm_slices <- readRDS('../../saved_states/wdm_cdm_raw/cdm_slices_raw.rds')
wdm_slices <- readRDS('../../saved_states/wdm_cdm_raw/wdm_slices_raw.rds')
indexmap <- readRDS('intermediate/indexmap.rds')

regroup_cube_robust <- function(slices, indexmap, n, min_find, max_find) {
  num_row <- 0
  for (i in 1:64) {
    num_row <- num_row + dim(slices[[i]])[1]
  }
  counter <- 1
  cube <- matrix(NA, num_row, 3)
  colors <- matrix(NA, num_row, 1)
  for (i in 1:64) {
    slice <- slices[[i]]
    map <- indexmap[[i]]
    num_in_slice <- dim(slice)[1]
    for (j in 1:3) {
      slice[,j] <- slice[,j] + map[j,1]
    }
    cube[counter:(counter+num_in_slice-1),] <- slice
    if (i %in% max_find) {
      colors[counter:(counter+num_in_slice-1),] <- rgb(1, 0, 0, 0.02)
    } else if (i %in% min_find) {
      colors[counter:(counter+num_in_slice-1),] <- rgb(0, 0, 1, 0.02)
    } else {
      colors[counter:(counter+num_in_slice-1),] <- rgb(0, 0, 0, 0.01)
    }
    counter <- counter + num_in_slice
  }
  return(list(cube=cube, colors=colors))
}


for (d in 1:2) {
  min_find <- minMat[d,]
  max_find <- maxMat[d,]
  cdm_obj <- regroup_cube_robust(cdm_slices, indexmap, 4, min_find, max_find)
  wdm_obj <- regroup_cube_robust(wdm_slices, indexmap, 4, min_find, max_find)

  cdm_cube <- cdm_obj$cube
  cdm_colors <- cdm_obj$colors
  wdm_cube <- wdm_obj$cube
  wdm_colors <- wdm_obj$colors

  png(paste("figure_14_wdm_colored_top_bot_3_dim_",d,".png", sep=""))
  scatterplot3d(wdm_cube, 
                xlab='', 
                ylab='', 
                zlab='', 
                pch='.',
                color=wdm_colors,
                tick.marks=FALSE,
                label.tick.marks=FALSE)
  dev.off()

  png(paste("figure_14_cdm_colored_top_bot_3_dim_",d,".png", sep=""))
  scatterplot3d(cdm_cube, 
                xlab='', 
                ylab='', 
                zlab='', 
                pch='.',
                color=cdm_colors,
                tick.marks=FALSE,
                label.tick.marks=FALSE)
  dev.off()
}
