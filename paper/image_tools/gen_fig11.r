# -- Figure 11 : cdm/wdm eagle

library('TDA')
library(scatterplot3d)
source('tools.r')
source("process_eagle.r")

whole_cdm <- load_CDM()
whole_wdm <- load_WDM()

data <- c(whole_wdm, whole_cdm)
names <- c('wdm', 'cdm')

for (i in c(1, 2)) {
    pdf(paste("figure_11_", names[i], "_plot.pdf"))
    scatterplot3d(data[i], 
                  xlab='X Axis', 
                  ylab='Y Axis', 
                  zlab='Z Axis', 
                  pch='.',
                  color=rgb(0, 0, 0, 0.01),
                  cex.axis=1.5,
                  cex.lab=2)
    dev.off()

    Xlim <- c(min(data[i][,1]), max(data[i][,1]))
    Ylim <- c(min(data[i][,2]), max(data[i][,2]))
    Zlim <- c(min(data[i][,3]), max(data[i][,3]))

    res <- 1
    diag <- gridDiag(data[i], 
                     dtm, 
                     lim=cbind(Xlim,Ylim,Zlim), 
                     by=res, 
                     sublevel=T, 
                     printProgress=T, 
                     m0=0.001)

    diag <- cleanDiag(diag)
    pdf(paste('figure_11_', names[i], '_pd.pdf'))
    X <- diag$diagram
    # X is the persistence diagram
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 1, 0, 0)) 
    plot(X[,2], 
         X[,3], 
         pch = c(X[,1]+1), 
         col = c(X[,1]+1), 
         xlab = 'Death', 
         ylab = 'Birth', 
         main = '', 
         cex.lab=2.0, 
         cex.axis=2.0, 
         cex.main=2.0, 
         cex.sub=2.0, 
         cex=2.0, 
         xlim=c(0,0.30), 
         ylim=c(0,0.30))
    abline(a = 0, b = 1)
    legend('bottomright', 
           c('0','1','2'), 
           pch = c(1,2,3), 
           col = c(1,2,3), 
           cex=2.0)
    dev.off()
}