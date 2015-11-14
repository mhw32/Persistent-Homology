source('VoronoiFoam.r')
N <- 10000
Boxlim <- c(0, 20)
for (i in 1:20) {
	voronoi_compilation(N=N, Boxlim=Boxlim, nameId=i)
}
