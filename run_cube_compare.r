source("process_eagle.r")
source("distance.r")
source("tools.r")

dists <- array(0,dim=c(4, 64, 3))
norm <- False

for (i in 1:4) {
  cdm_slices <- readRDS(paste('cdm_diags_', i, '.rds', sep=''))
  wdm_slices <- readRDS(paste('wdm_diags_', i, '.rds', sep=''))
  for (l in 1:i^3) {
    cdm_diag <- cleanDiag(cdm_slices[[l]])
    wdm_diag <- cleanDiag(wdm_slices[[l]])

    if (norm) { # optionally normalize
      cdm_diag <- normalize(cdm_diag)
      wdm_diag <- normalize(wdm_diag)
    }

    # get the bottle neck distances.
    dists[i,l,1] <- bottleneckDist(cdm_diag, wdm_diag, dimension=1)
    dists[i,l,2] <- bottleneckDist(cdm_diag, wdm_diag, dimension=2)
    dists[i,l,3] <- bottleneckDist(cdm_diag, wdm_diag, dimension=3)
  }
}

saveRDS(dists, file='outputs/dists/dists.rds')