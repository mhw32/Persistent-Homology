
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

max_index <- 41
max_cdm <- cdm_diags[[max_index]]
max_wdm <- wdm_diags[[max_index]]

min_index <- 24
min_cdm <- cdm_diags[[min_index]]
min_wdm <- wdm_diags[[min_index]]

max_tseq1 <- seq(min(max_cdm[,2:3]),max(max_cdm[,2:3]),length=5000)
max_euler1 <- eulerChar(max_tseq1, max_cdm, maxdimension=max(max_cdm[,1]), threshold=0)

max_tseq2 <- seq(min(max_wdm[,2:3]),max(max_wdm[,2:3]),length=5000)
max_euler2 <- eulerChar(max_tseq2, max_wdm, maxdimension=max(max_wdm[,1]), threshold=0)

min_tseq1 <- seq(min(min_cdm[,2:3]),max(min_cdm[,2:3]),length=5000)
min_euler1 <- eulerChar(min_tseq1, min_cdm, maxdimension=max(min_cdm[,1]), threshold=0)

min_tseq2 <- seq(min(min_wdm[,2:3]),max(min_wdm[,2:3]),length=5000)
min_euler2 <- eulerChar(min_tseq2, min_wdm, maxdimension=max(min_wdm[,1]), threshold=0)

pdf('figure_13_maxmin_margin_euler.pdf')
par(mar=c(5,6,4,2))
plot(max_euler2, 
     type="l", 
     xlab="t", 
     ylab="EC", 
     col=rgb(63, 223, 218, 255, maxColorValue=255),
     lwd=1.5,
     cex.lab=2.5, 
     cex.axis=2.5, 
     cex.main=2.5, 
     cex.sub=2.5,
     xlim=c(0, 2000))
lines(max_euler1, lwd=1.5, col=rgb(248, 137, 113 ,255, maxColorValue=255))
# lines(min_euler2, lwd=1.5, lty="solid", col=rgb(202, 109, 89 ,255, maxColorValue=255))
# lines(min_euler1, lwd=1.5, lty="solid", col=rgb(46, 165, 162, 255, maxColorValue=255))
legend("topright", 
       c("WDM", "CDM"), # , "CDM (Low)", "WDM (Low)"), 
       col=c(rgb(63, 223, 218, 255, maxColorValue=255), rgb(248, 137, 113 ,255, maxColorValue=255)), # , #rgb(202, 109, 89 ,255, maxColorValue=255), rgb(46, 165, 162, 255, maxColorValue=255)), 
       lty=c("solid", "solid"), # "solid", "solid"),
       lwd=1.5,
       cex=2.5)
dev.off()


# plot WDM and CDM plots
library(scatterplot3d)
cdm_41 <- readRDS('intermediate/fig_13_wdm_41_slice.rds')
wdm_41 <- readRDS('intermediate/fig_13_cdm_41_slice.rds')

png("figure_13_cdm_slice.png")
scatterplot3d(cdm_41, 
              xlab='', 
              ylab='', 
              zlab='', 
              pch='.',
              color=rgb(0, 0, 0, 0.5),
              tick.marks=FALSE,
              label.tick.marks=FALSE)
dev.off()
png("figure_13_wdm_slice.png")
scatterplot3d(wdm_41,
              xlab='', 
              ylab='', 
              zlab='',
              pch='.',
              color=rgb(0, 0, 0, 0.5),
              tick.marks=FALSE,
              label.tick.marks=FALSE) 
dev.off()


