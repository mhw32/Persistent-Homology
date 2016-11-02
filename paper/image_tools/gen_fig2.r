# -- Figure 2 : persistence diagram --
library(TDA)
source('../../tools.r')

infile<- 'intermediate/fig_2_data.bin'
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
Mat <- matrix(readBin(con, "numeric", prod(dim)), 
              dim[1], 
              dim[2])
close(con)

res <- 0.05
boxlim <- cbind(c(-1.5, 1.5), 
                c(-1.5, 1.5))
diag <- gridDiag(Mat, 
                 kde, 
                 lim=boxlim, 
                 by=res, 
                 sublevel=T, 
                 printProgress=T,
                 h=.3)

diag <- cleanDiag(diag$diagram)

pdf('figure_2_pd.pdf')
X <- diag
# X is the persistence diagram
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0)) 
plot(X[,2], 
     X[,3], 
     pch = c(X[,1]+1), 
     col = c(X[,1]+1), 
     xlab = "Death", 
     ylab = "Birth", 
     main = "", 
     cex.lab=2.0, 
     cex.axis=2.0, 
     cex.main=2.0, 
     cex.sub=2.0, 
     cex=2.0)
abline(a = 0, b = 1)
legend("bottomright", 
       c("0","1","2"), 
       pch = c(1,2,3), 
       col = c(1,2,3), 
       cex=2.0)
dev.off()
