# -- Figure 13 : heatmap analysis (euler) --

library(TDA)
library(pracma)

eulerChar <- function(tseq, diagram, maxdimension, threshold=0) {
    eulerFct <- function(t) {
        if (threshold>0) {
            persistence <- diagram[,3]-diagram[,2]
            diagram <- diagram[which(persistence>=threshold),]
        }
        betti <- numeric(maxdimension)
        for (i in 0:maxdimension) {
            betti[i+1]=((-1)^i)*sum(diagram[,1]==i & diagram[,2]<=t & diagram[,3]>t)
        }
        out <- sum(betti)
        return(out)
    }
    out <- numeric(length(tseq))
    count <- 0
    for (t in tseq) {
        count <- count + 1
        out[count] <- eulerFct(t)
    }
    return(out)
}


cdm_diags <- readRDS('intermediate/fig_13_cdm_diags_4.rds')
wdm_diags <- readRDS('intermediate/fig_13_wdm_diags_4.rds')

# max margin

index <- 41
my_cdm <- cdm_diags[[index]]
my_wdm <- wdm_diags[[index]]

tseq1 <- seq(min(my_cdm[,2:3]),max(my_cdm[,2:3]),length=5000)
euler1 <- eulerChar(tseq1, my_cdm, maxdimension=max(my_cdm[,1]), threshold=0)

tseq2 <- seq(min(my_wdm[,2:3]),max(my_wdm[,2:3]),length=5000)
euler2 <- eulerChar(tseq2, my_wdm, maxdimension=max(my_wdm[,1]), threshold=0)

pdf('figure_13_max_margin_2euler.pdf')
plot(euler2, 
     type="l", 
     xlab="Grid Sequence", 
     ylab="Euler Char.", 
     col="coral1",
     lwd=3,
     cex.lab=2.0, 
     cex.axis=2.0, 
     cex.main=2.0, 
     cex.sub=2.0, 
     cex=2.0)
lines(euler1,
      lwd=3, 
      col="cyan3")
legend("topright", 
       c("CDM", "WDM"), 
       col=c("coral", "cyan"), 
       lwd=3)
dev.off()

# min margin

index <- 24
my_cdm <- cdm_diags[[index]]
my_wdm <- wdm_diags[[index]]

tseq1 <- seq(min(my_cdm[,2:3]),max(my_cdm[,2:3]),length=5000)
euler1 <- eulerChar(tseq1, my_cdm, maxdimension=max(my_cdm[,1]), threshold=0)

tseq2 <- seq(min(my_wdm[,2:3]),max(my_wdm[,2:3]),length=5000)
euler2 <- eulerChar(tseq2, my_wdm, maxdimension=max(my_wdm[,1]), threshold=0)

pdf('figure_13_min_margin_2euler.pdf')
plot(euler2, 
     type="l", 
     xlab="Grid Sequence", 
     ylab="Euler Char.", 
     col="coral1",
     lwd=3,
     cex.lab=2.0, 
     cex.axis=2.0, 
     cex.main=2.0, 
     cex.sub=2.0, 
     cex=2.0)
lines(euler1,
      lwd=3, 
      col="cyan3")
legend("topright", 
       c("CDM", "WDM"), 
       col=c("coral", "cyan"), 
       lwd=3,
       cex.lab=2.0, 
       cex.axis=2.0, 
       cex.main=2.0, 
       cex.sub=2.0, 
       cex=2.0)
dev.off()
