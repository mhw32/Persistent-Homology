#!/usr/bin/r

infile <- "intermediate/fig_14_data.bin"
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
Mat <- matrix(readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

cdm <- readRDS('../../saved_states/wdm_cdm_raw/cdm_slices_raw.rds')
wdm <- readRDS('../../saved_states/wdm_cdm_raw/wdm_slices_raw.rds')
Mat_dim <- dim(Mat)

num_cdm_row <- 0
num_wdm_row <- 0
for (grid_i in 1:64) {
  num_wdm_row <- num_wdm_row + dim(wdm[[grid_i]])[1]
  num_cdm_row <- num_cdm_row + dim(cdm[[grid_i]])[1]
}

wdm_mat <- matrix(, nrow=num_wdm_row, ncol=3, byrow=TRUE)
wdm_col <- matrix(, nrow=num_wdm_row, ncol=1, byrow=TRUE)
cdm_mat <- matrix(, nrow=num_cdm_row, ncol=3, byrow=TRUE)
cdm_col <- matrix(, nrow=num_cdm_row, ncol=1, byrow=TRUE)

dim_i <- 1
m <- Mat[dim_i,]
counter_wdm <- 1
counter_cdm <- 1

for (grid_i in 1:64) {
  wdm_slice <- wdm[[grid_i]]
  num_wdm_slice <- dim(wdm[[grid_i]])[1]
  wdm_mat[counter_wdm:(counter_wdm+num_wdm_slice-1), 1:3] <- wdm_slice

  cdm_slice <- cdm[[grid_i]]
  num_cdm_slice <- dim(cdm[[grid_i]])[1]
  cdm_mat[counter_cdm:(counter_cdm+num_cdm_slice-1), 1:3] <- cdm_slice

  if (grid_i %in% m) {
    wdm_col[counter_wdm:(counter_wdm+num_wdm_slice-1), 1] <- rgb(1, 0, 0, 0.06)
    cdm_col[counter_cdm:(counter_cdm+num_cdm_slice-1), 1] <- rgb(1, 0, 0, 0.06)
  } else {
    wdm_col[counter_wdm:(counter_wdm+num_wdm_slice-1), 1] <- rgb(0, 0, 0, 0.03)
    cdm_col[counter_cdm:(counter_cdm+num_cdm_slice-1), 1] <- rgb(0, 0, 0, 0.03)
  }
  counter_wdm <- counter_wdm + num_wdm_slice
  counter_cdm <- counter_cdm + num_cdm_slice
}

library(scatterplot3d)
png("figure_14_wdm_colored.png")
scatterplot3d(wdm_mat, 
              xlab='', 
              ylab='', 
              zlab='', 
              pch='.',
              color=wdm_col,
              tick.marks=FALSE,
              label.tick.marks=FALSE)
dev.off()

png("figure_14_cdm_colored.png")
scatterplot3d(cdm_mat, 
              xlab='', 
              ylab='', 
              zlab='', 
              pch='.',
              color=cdm_col,
              tick.marks=FALSE,
              label.tick.marks=FALSE)
dev.off()
