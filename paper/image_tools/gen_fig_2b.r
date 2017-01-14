library(TDA)

X <- readRDS('intermediate/fig_2_tmp.rds')

# Construct a grid of points over which we can eval the function
by <- 0.025
Xseq <- seq(-1.5, 1.5, by=by)
Yseq <- seq(-1.5, 1.5, by=by)
Grid <- expand.grid(Xseq, Yseq)

# distance to a measure
m0 <- 0.1
DTM <- dtm(X, Grid, m0)

saveRDS(DTM, 'intermediate/fig_2_dtm.rds')
