''' Measure a statistic by two-point correlation function. 
	Use NN (function of number counts).

	For the NN and NNN correlations, the raw calculation 
	is not sufficient to produce the real correlation function. 
	You also need to account for the survey geometry (edges, 
	mask, etc.) by running the same calculation with a random 
	catalog (or several) that have a uniform density, but the 
	same geometry.

	This calculates xi = (DD-2DR+RR)/RR for each bin. This is 
	the Landy-Szalay estimator, which is the most widely used 
	estimator for count-count correlation functions. However, 
	if you want to use a simpler estimator i = (DD-RR)/RR, then 
	you can omit the dr parameter. The simpler estimator is 
	slightly biased though, so this is not recommended.
'''

from __future__ import print_function
import numpy as np
import numpy.random as npr
import treecorr
from scipy.integrate import simps


def get_corr(x, y, z, min_r=1., max_r=100., L=500., Lc=-0.5, ngal=10000, simple=False):
	'''
	Use a simple probability distribution for the galaxies:

	n(r) = (2pi s^2)^-1 exp(-r^2/2s^2)

	The Fourier transform is: n~(k) = exp(-s^2 k^2/2)
	P(k) = <|n~(k)|^2> = exp(-s^2 k^2)
	xi(r) = (1/2pi) int( dk k P(k) J0(kr) ) 
	      = 1/(4 pi s^2) exp(-r^2/4s^2)

	However, we need to correct for the uniform density background, so the real result
	is this minus 1/L^2 divided by 1/L^2.  So:

	xi(r) = 1/4pi (L/s)^2 exp(-r^2/4s^2) - 1
	'''
	nrand = ngal

	# initialize a catalog
	cat = treecorr.Catalog(x=x, y=y, z=z)
	# initialize a NN correlation
	dd = treecorr.NNCorrelation(bin_size=0.1, min_sep=min_r, max_sep=max_r, 
															sep_units='arcmin', verbose=2)
	# process the dd
	dd.process(cat)

	# initialize some random states
	rx = npr.random_sample(nrand) * L + Lc
	ry = npr.random_sample(nrand) * L + Lc
	rz = npr.random_sample(nrand) * L + Lc

	# calculate catalog and NN for random one
	rand = treecorr.Catalog(x=rx, y=ry, z=rz)
	rr = treecorr.NNCorrelation(bin_size=0.1, min_sep=min_r, max_sep=max_r, 
															sep_units='arcmin', verbose=2)

	# process the rr
	rr.process(rand)

	if simple: 
		# calculate the correlation function (simple version)
		xi, varxi = dd.calculateXi(rr)
	else:
		# relate the dd and rr
		dr = treecorr.NNCorrelation(bin_size=0.1, min_sep=min_r, max_sep=max_r, 
																sep_units='arcmin', verbose=2)
		dr.process(cat, rand)

		# calculate the correlation function
		xi, varxi = dd.calculateXi(rr, dr)
	return xi, varxi 


def get_corr_func(data, min_r=1., max_r=100., L=30., Lc=-0.5, ngal=5000, simple=False):
	x, y, z = data[:, 0], data[:, 1], data[:, 2]
	corr_fun, corr_var = get_corr(x, y, z, min_r=min_r, max_r=max_r, L=L, ngal=ngal, simple=simple)
	return corr_fun


def get_corr_stat(data, min_r=1., max_r=100., L=30., Lc=-0.5, ngal=5000, simple=False):
	x, y, z = data[:, 0], data[:, 1], data[:, 2]
	corr_fun, corr_var = get_corr(x, y, z, min_r=min_r, max_r=max_r, L=L, ngal=ngal, simple=simple)
	corr_fun = np.abs(corr_fun)
	return simps(corr_fun, dx=0.01)


