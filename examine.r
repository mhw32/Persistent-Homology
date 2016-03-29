# Automate the process of testing this stuff!
# This will save me a lot of time in the case that I have to generate lots of new images.

source('tools.r')
source('VoronoiFoam.r')
source('euler.r')
source('localtest.r')

plotDiag <- function(diag, main, path) {
  png(filename=path)
  plot(diag, main=main)
  dev.off()
}

cleanDiag <- function(diag) { diag[2:nrow(diag),] }

boxlimTest <- function() {
  # Plot the visuals
  voronoi_plot(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=10000, main="(10k,50box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-vtest.png")
  voronoi_plot(boxlim=c(0,10),percFil=0.4,res=0.5,perturb=1,N=10000, main="(10k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k05-vtest.png")
  voronoi_plot(boxlim=c(0,10),percFil=0.4,res=0.1,perturb=1,N=10000, main="(10k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k01-vtest.png")
  voronoi_plot(boxlim=c(0,10),percFil=0.4,res=0.5,perturb=1,N=2000, main="(2k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k05-vtest.png")
  voronoi_plot(boxlim=c(0,10),percFil=0.4,res=0.1,perturb=1,N=2000, main="(2k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k01-vtest.png")
  # Save the diagrams.
  diagbase <- voronoi_example(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=10000)$diagram
  diag10k05 <- voronoi_example(boxlim=c(0,10),percFil=0.4,res=0.5,perturb=1,N=10000)$diagram
  diag10k01 <- voronoi_example(boxlim=c(0,10),percFil=0.4,res=0.1,perturb=1,N=10000)$diagram
  diag2k05 <- voronoi_example(boxlim=c(0,10),percFil=0.4,res=0.5,perturb=1,N=2000)$diagram
  diag2k01 <- voronoi_example(boxlim=c(0,10),percFil=0.4,res=0.1,perturb=1,N=2000)$diagram
  # Normalize the diagrams.
  diagbase <- normalize(diagbase)
  diag10k05 <- normalize(diag10k05)
  diag10k01 <- normalize(diag10k01)
  diag2k05 <- normalize(diag2k05)
  diag2k01 <- normalize(diag2k01)
  # Plot these.
  plotDiag(diagbase, main="(10k,50box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-diag-vtest.png")
  plotDiag(diag10k05, main="(10k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k05-diag-vtest.png")
  plotDiag(diag10k01, main="(10k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k01-diag-vtest.png")
  plotDiag(diag2k05, main="(2k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k05-diag-vtest.png")
  plotDiag(diag2k01, main="(2k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k01-diag-vtest.png")
  # Clean the diagrams.
  diagbase <- cleanDiag(diagbase)
  diag10k05 <- cleanDiag(diag10k05)
  diag10k01 <- cleanDiag(diag10k01)
  diag2k05 <- cleanDiag(diag2k05)
  diag2k01 <- cleanDiag(diag2k01)
  # Calculate euler tests.
  eulerPlot(diagbase, main="(10k,50box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-euler-vtest.png")
  eulerPlot(diag10k05, main="(10k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k05-euler-vtest.png")
  eulerPlot(diag10k01, main="(10k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/10k01-euler-vtest.png")
  eulerPlot(diag2k05, main="(2k,10box,0.5res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k05-euler-vtest.png")
  eulerPlot(diag2k01, main="(2k,10box,0.1res,0.4fil)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2k01-euler-vtest.png")
  # Now do local KDE plots.
  localdiagplot(diagbase, diag10k05, dim=0, main="(10k,10box,0.5res,0.4fil,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k05d0-local-vtest.png")
  localdiagplot(diagbase, diag10k05, dim=1, main="(10k,10box,0.5res,0.4fil,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k05d1-local-vtest.png")
  localdiagplot(diagbase, diag10k05, dim=2, main="(10k,10box,0.5res,0.4fil,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k05d2-local-vtest.png")
  localdiagplot(diagbase, diag10k01, dim=0, main="(10k,10box,0.1res,0.4fil,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k01d0-local-vtest.png")
  localdiagplot(diagbase, diag10k01, dim=1, main="(10k,10box,0.1res,0.4fil,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k01d1-local-vtest.png")
  localdiagplot(diagbase, diag10k01, dim=2, main="(10k,10box,0.1res,0.4fil,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-10k01d2-local-vtest.png")
  localdiagplot(diagbase, diag2k05, dim=0, main="(2k,10box,0.5res,0.4fil,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k05d0-local-vtest.png")
  localdiagplot(diagbase, diag2k05, dim=1, main="(2k,10box,0.5res,0.4fil,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k05d1-local-vtest.png")
  localdiagplot(diagbase, diag2k05, dim=2, main="(2k,10box,0.5res,0.4fil,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k05d2-local-vtest.png")
  localdiagplot(diagbase, diag2k01, dim=0, main="(2k,10box,0.1res,0.4fil,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k01d0-local-vtest.png")
  localdiagplot(diagbase, diag2k01, dim=1, main="(2k,10box,0.1res,0.4fil,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k01d1-local-vtest.png")
  localdiagplot(diagbase, diag2k01, dim=2, main="(2k,10box,0.1res,0.4fil,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2k01d2-local-vtest.png")
}

particleTest <- function () {
  # Plot the visuals
  voronoi_plot(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=10000, main="(10k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-ptest.png")
  voronoi_plot(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=5000, main="(5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/5k-ptest.png")
  voronoi_plot(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=2500, main="(2.5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2.5k-ptest.png")
  voronoi_plot(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=1000, main="(1k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/1k-ptest.png")
  # Save the diagrams.
  diagbase <- voronoi_example(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=10000)$diagram
  diag5k <- voronoi_example(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=5000)$diagram
  diag2k <- voronoi_example(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=2500)$diagram
  diag1k <- voronoi_example(boxlim=c(0,50),percFil=0.4,res=0.5,perturb=1,N=1000)$diagram
  # Normalize the diagrams.
  diagbase <- normalize(diagbase)
  diag5k <- normalize(diag5k)
  diag2k <- normalize(diag2k)
  diag1k <- normalize(diag1k)
  # Plot these.
  plotDiag(diagbase, main="(10k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-diag-ptest.png")
  plotDiag(diag5k, main="(5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/5k-diag-ptest.png")
  plotDiag(diag2k, main="(2.5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2.5k-diag-ptest.png")
  plotDiag(diag1k, main="(1k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/1k-diag-ptest.png")
  # Clean the diagrams.
  diagbase <- cleanDiag(diagbase)
  diag5k <- cleanDiag(diag5k)
  diag2k <- cleanDiag(diag2k)
  diag1k <- cleanDiag(diag1k)
  # Calculate euler tests.
  eulerPlot(diagbase, main="(10k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-euler-ptest.png")
  eulerPlot(diag5k, main="(5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/5k-euler-ptest.png")
  eulerPlot(diag2k, main="(2.5k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/2.5k-euler-ptest.png")
  eulerPlot(diag1k, main="(1k,50box)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/1k-euler-ptest.png")
  # Now do local KDE plots.
  localdiagplot(diagbase, diag5k, dim=0, main="(5k,50box,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-5kd0-local-ptest.png")
  localdiagplot(diagbase, diag5k, dim=1, main="(5k,50box,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-5kd1-local-ptest.png")
  localdiagplot(diagbase, diag5k, dim=2, main="(5k,50box,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-5kd2-local-ptest.png")
  localdiagplot(diagbase, diag2k, dim=0, main="(2.5k,50box,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2.5kd0-local-ptest.png")
  localdiagplot(diagbase, diag2k, dim=1, main="(2.5k,50box,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2.5kd1-local-ptest.png")
  localdiagplot(diagbase, diag2k, dim=2, main="(2.5k,50box,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-2.5kd2-local-ptest.png")
  localdiagplot(diagbase, diag1k, dim=0, main="(1k,50box,0dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-1kd0-local-ptest.png")
  localdiagplot(diagbase, diag1k, dim=1, main="(1k,50box,1dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-1kd1-local-ptest.png")
  localdiagplot(diagbase, diag1k, dim=2, main="(1k,50box,2dim)", path="/Users/grub/Desktop/Cisewski-Lab/saved_states/test_samples/base-1kd2-local-ptest.png")
}
