source('VoronoiFoam.r')
N <- 2000
Boxlim <- c(0, 50)
res <- 0.5
perturb <- 1
groupN <- 10
baseline <- 0.1
for (i in 4:7) {
	voronoi_compilation(N=N, Boxlim=Boxlim, res=res, perturb=perturb, groupN=groupN, baseline=baseline, nameId=i)
}

