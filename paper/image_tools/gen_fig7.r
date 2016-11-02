# -- Figure 7 : different voronoi sims and pd diags --

library(TDA)
source('../../tools.r')
library(scatterplot3d)
source('../../Voronoi3Dfct.r')

boxlim <- c(0, 20)
res <- 0.5
perturb <- 1
Xlim <- boxlim
Ylim <- boxlim
Zlim <- boxlim
N <- 10000

pfs <- c(0.1, 0.5, 0.9)

for (pf in pfs) {
    vf <- voronoi3d(boxlim,
                    res,
                    perturb,
                    Ncells=64,
                    N,
                    percClutter=0,
                    percWall=1-0.02-pf,
                    percFil=pf,
                    percClust=0.02)

    pdf(paste('figure_7_plot_pf_', pf, '.pdf', sep=""))
    scatterplot3d(vf,
                  xlab='',
                  ylab='',
                  zlab='',
                  pch='.',
                  color=rgb(0, 0, 0, 0.25),
                  tick.marks=FALSE,
                  label.tick.marks=FALSE)
    dev.off()

    diag <- gridDiag(vf,
                     dtm,
                     lim=cbind(Xlim,Ylim,Zlim),
                     by=res,
                     sublevel=T,
                     printProgress=T,
                     m0=0.001)

    diag <- cleanDiag(diag$diagram)

    pdf(paste('figure_7_pd_', pf, '.pdf', sep=""))
    X <- diag
    # X is the persistence diagram
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 1, 0, 0))
    plot(X[,2],
         X[,3],
         pch = c(X[,1]+1),
         col = c(X[,1]+1),
         xlab = 'Birth',
         ylab = 'Death',
         main = '',
         cex.lab=2.0,
         cex.axis=2.0,
         cex.main=2.0,
         cex.sub=2.0,
         cex=2.0)
    abline(a = 0, b = 1)
    legend('bottomright',
           c('0','1','2'),
           pch = c(1,2,3),
           col = c(1,2,3),
           cex=2.0)
    dev.off()
}
