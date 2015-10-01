library(FNN)
library(mvtnorm)

## INPUT
##############################
## Boxlim: limit of the side of the box
## resolution: grid space for approximating the voronoi cells.
## perturb: sd of gaussian noise to perturb the grid
## Ncells: number of voronoi cells
## N: total number of particles
## percClutter: percentage of clutter noise
## percWall: percentage of particles on the walls
## percFil: percentage of particles on the filaments
## percClust: percentage of particles on the clusters

voronoi3d <- function(Boxlim, resolution, perturb, Ncells, N, percClutter=0.05, percWall=0.8, percFil=0.1, percClust=0.05) {
	
	Xlim <- Boxlim
	Ylim <- Boxlim
	Zlim <- Boxlim

	## generate random nuclei inside the box
	cells <- cbind(runif(Ncells,Xlim[1],Xlim[2]),
				runif(Ncells,Ylim[1],Ylim[2]),
				runif(Ncells,Zlim[1],Zlim[2]))


	## construct the grid over which we will evaluate voronoi cells
	by <- resolution
	Xseq <- seq(Xlim[1],Xlim[2], by=by)
	Yseq <- seq(Ylim[1],Ylim[2], by=by)
	Zseq <- seq(Zlim[1],Zlim[2], by=by)
	Grid <- expand.grid(Xseq,Yseq,Zseq)

	# we find the 3 closest nuclei for each point of the grid
	closest3 <- get.knnx(cells, Grid, k=3, algorithm=c("kd_tree"))

	# ids of voronoi cells
	id <- closest3$nn.index[,1]

	# we check the neighborhood of each point of the grid
	# to see how many different ids we find (stored in what)
	closestGrid <- get.knnx(Grid, Grid, k=8, algorithm=c("kd_tree"))$nn.index
	what <- sapply(1:nrow(closestGrid), FUN=function(i) { 
								ids=length(unique(id[closestGrid[i,]]))
								return(ids) })
	id.cluster <- which(what>3)
	id.filament <- which(what==3)
	id.wall <- which(what==2)

	# plot3d(Grid[id.wall,])
	# points3d(Grid[id.filament,], col=2, size=5)
	# points3d(Grid[id.cluster,], col=4, size=10)

	# Gaussian noise to perturb the grid
	noise <- rmvnorm(N,rep(0,3), diag(perturb^2,3))

	clusterPoints <- NULL
	filPoints <- NULL
	wallPoints <- NULL
	clutterPoints <- NULL

	# Generate clusters
	Nclust <- percClust*N
	if (Nclust > 0) {
		clusterPoints <- Grid[sample(id.cluster,Nclust, replace=T),]+noise[1:Nclust,]
	}

	# Generate filaments
	Nfil <- percFil*N
	if (Nfil > 0) {
		filPoints <- Grid[sample(id.filament,Nfil, replace=T),]+noise[(Nclust+1):(Nclust+Nfil),]
	}

	# Generate walls
	Nwall <- percWall*N
	if (Nwall > 0) {
		wallPoints <- Grid[sample(id.wall,Nwall, replace=T),]+noise[(Nclust+Nfil+1):(Nclust+Nfil+Nwall),]
	}

	# Generate clutter
	Nclutter <- percClutter*N
	if (Nclutter > 0){
		clutterPoints <- bind(runif(Nclutter, Xlim[1], Xlim[2]),
							runif(Nclutter, Ylim[1], Ylim[2]),
							runif(Nclutter, Zlim[1], Zlim[2]))
		colnames(clutterPoints)=c("Var1", "Var2", "Var3")
	}
	out <- rbind(clusterPoints, filPoints, wallPoints, clutterPoints)
	rownames(out)=NULL

	return(out)
}


