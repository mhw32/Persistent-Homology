source('VoronoiFoam.r')
N <- 5000
Boxlim <- c(0, 50)
res <- 0.5
perturb <- 1
groupN <- 15
baseline <- 0.1
for (i in 1:100) {
	print(paste("Iteration : ", i))
	voronoi_compilation(N=N, Boxlim=Boxlim, res=res, perturb=perturb, groupN=groupN, nameId=i)
}

