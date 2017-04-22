library(TDA)
library(pracma)
source('tools.r')

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

eulerCharDim <- function(tseq, diagram, which.dim, threshold=0) {
	eulerFct <- function(t) {
		if (threshold > 0) {
			persistence <- diagram[,3] - diagram[,2]
			diagram <- diagram[which(persistence>=threshold),]
		}
		betti <- sum(diagram[,1] == which.dim & diagram[,2] <= t & diagram[,3] > t)
		return(betti)
	}
	out <- numeric(length(tseq))
	count <- 0
	for (t in tseq) {
		count <- count + 1
		out[count] <- eulerFct(t)
	}
	return(out)
}

# Find the area underneath a curve.
integrate <- function(xarr, yarr) {
	auc <- trapz(xarr, yarr)
	return(auc)
}

eulerKernel <- function(P, N=1000) {
	tseq <- seq(min(P[,2:3]),max(P[,2:3]),length=N)
	euler <- eulerChar(tseq, P, maxdimension=max(P[,1]), threshold=0)
	auc <- integrate(tseq, euler)
	return(auc)
}

eulerStat <- function(X, L) {
	G1 <- X[L == 0]
	G2 <- X[L == 1]
  n <- length(G1)
  m <- length(G2)
  G1arr <- sapply(seq(1:n), function(i) { eulerKernel(G1) })
  G2arr <- sapply(seq(1:m), function(i) { eulerKernel(G2) })
  proba <- t.test(G1arr, G2arr, conf.level=0.95, paired=TRUE)
  return(proba$p.value)
}

eulerVanillaFunction <- function(diagram, length=1000) {
  tseq <- seq(min(diagram[,2:3]),max(diagram[,2:3]),length=length)
  euler <- eulerChar(tseq, diagram, maxdimension=max(diagram[,1]), threshold=0)
  return euler;
}

# Wrapper fxn to combine generation and integration.
eulerIntegration <- function(diagram) {
  tseq <- seq(min(diagram[,2:3]),max(diagram[,2:3]),length=1000)
  euler <- eulerChar(tseq, diagram, maxdimension=max(diagram[,1]), threshold=0)
  # Always return absolute value of the euler characteristic.
  auc <- integrate(tseq, abs(euler))
  return(auc)
}

eulerDimIntegration <- function(diagram, dim=1) {
	tseq <- seq(min(diagram[,2:3]),max(diagram[,2:3]),length=1000)
  euler <- eulerCharDim(tseq, diagram, dim, threshold=0)
  # Always return absolute value of the euler characteristic.
  auc <- integrate(tseq, abs(euler))
  return(auc)
}

eulerPlot <- function(diagram, main='', path='') {
	tseq <- seq(min(diagram[,2:3]),max(diagram[,2:3]),length=5000)
	euler <- eulerChar(tseq, diagram, maxdimension=max(diagram[,1]), threshold=0)
	png(filename=path)
	plot(euler, type="l", cex=2.0)
	dev.off()
}

eulerDualPlot <- function(diagram1, diagram2, main='', path='') {
  tseq1 <- seq(min(diagram1[,2:3]),max(diagram1[,2:3]),length=5000)
  euler1 <- eulerChar(tseq1, diagram1, maxdimension=max(diagram1[,1]), threshold=0)

  tseq2 <- seq(min(diagram2[,2:3]),max(diagram2[,2:3]),length=5000)
  euler2 <- eulerChar(tseq2, diagram2, maxdimension=max(diagram2[,1]), threshold=0)

  png(filename=path)
  plot(euler1, type="l", xlab="Grid Sequence", ylab="Euler Char.", col="blue", main=main)
  lines(euler2, col="red")
  legend("topright", c("CDM", "WDM"), col=c("blue", "red"), lty=1)
  dev.off()
}
