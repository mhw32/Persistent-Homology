rm(list=ls(all=TRUE))  # clear all variables


library(TDA)

### make a landscape function

Xdiag <- rbind(c(.1, .35), c(.3, .6), c(.2, .7), c(.7, .75), c(.2, .6), c(.5, .75), c(.25, .3), c(.1, .6), c(.05, .2))
plot(Xdiag)

lx <- landscape(cbind(rep(0,nrow(Xdiag)),Xdiag), dimension = 0, KK = 1:10, tseq = seq(min(Xdiag), max(Xdiag), length=500))
landx <- cbind((Xdiag[,1] + Xdiag[,2])/2, (Xdiag[,2] - Xdiag[,1])/2)
#
sx <- silhouette(cbind(rep(0,nrow(Xdiag)),Xdiag), p = 1, dimension = 0, tseq = seq(min(Xdiag), max(Xdiag), length=500))
#



par(mfrow = c(1,2), mar = c(4,5,2,2))
plot(Xdiag[,1], Xdiag[,2], pch = 1, col = "black", lwd = 3, cex = 1.5, xlab = "Death", ylab = "Birth", cex.axis = 2, cex.lab = 2, main = "(a) Persistence Diagram", xlim = range(Xdiag), ylim = range(Xdiag), cex.main = 2)
abline(a = 0, b = 1, col = "black", lwd = 2)
legend("bottomright", expression(H[0]), col = "black", pch = 1,  cex = 2, lwd = 3, lty = NA)
#
plot(seq(min(Xdiag), max(Xdiag), length=500), lx[,1], "l", lwd = 1, xlab = "(Birth+Death)/2", ylab = "(Birth-Death)/2", cex.axis = 2, cex.lab = 2, cex = 1.5, main = "(b) Silhouette and Landscape functions", cex.main = 2)
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,2], col = "black", lwd = 1)
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,3], col = "black", lwd = 1)
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,4], col = "black", lwd = 1)
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,5], col = "black", lwd = 1)
#
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,1], col = "magenta", lwd = 5, lty = 2)
lines(seq(min(Xdiag), max(Xdiag), length=500), lx[,2], col = "blue", lwd = 5, lty = 3)
lines(seq(min(Xdiag), max(Xdiag), length=500), sx, col = "cyan", lwd = 5, lty = 1)
#
points(landx[,1], landx[,2], pch = 1, lwd = 3, cex = 1.5)
legend("topleft", c(expression(phi[0]), expression(lambda(1,)), expression(lambda(2,)), expression(lambda(3:5,))), col = c("cyan", "magenta", "blue", "black"), lty = c(1,2,3,1), lwd = c(3,3,3,1), cex = 2)
