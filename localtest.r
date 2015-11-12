library(ks)
source('tools.r')

# Read in two sample persistence diagrams.
big <- readRDS('./saved_states/diagbig.rds')
small <- readRDS('./saved_states/diagsmall.rds')
diagbig <- big$diagram
diagsmall <- small$diagram
# Pull out the dimensions of two diagrams.
dim <- 1
dimbig <- diagbig[,2:3][diagbig[,1] == dim]
dimsmall <- diagsmall[,2:3][diagbig[,1] == dim]
# Run the local significance test.
test <- kde.local.test(dimsmall, dimbig) 
# ^ really unsure about this.