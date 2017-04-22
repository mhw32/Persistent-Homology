source('VoronoiFoam.r')
N <- 5000
Boxlim <- c(0, 20)
res <- 0.5
perturb <- 1
groupN <- 15
baseline <- 0.1
for (i in 80:100) {
	voronoi_compilation(N=N, Boxlim=Boxlim, res=res, perturb=perturb, groupN=groupN, baseline=baseline, nameId=i)
}

