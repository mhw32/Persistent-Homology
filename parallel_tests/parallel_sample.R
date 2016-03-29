library(snow)
library(Rmpi)

np <- mpi.universe.size()
cluster <- makeMPIcluster(np)

sayhello <- function() {
    info <- Sys.info()[c("nodename", "machine")]
    paste("Hello from", info[1], "with CPU type", info[2])
}

names <- clusterCall(cluster, sayhello)
stopCluster(cluster)
mpi.exit()
