#gaussian fit

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def fitgauss(bin_edges, counts, guess):
	"""fit offset gaussian, with offset, amplitude, and sigma constrained to be positive"""
	bincen = np.zeros(np.size(bin_edges)-1)
	bincen[:] = bin_edges[0:np.size(bin_edges)-1]
	binstep = bin_edges[1]-bin_edges[2]
	bincen -= binstep/2
	#print(bincen)
	#print(counts)
	#plt.plot(bincen, counts, label = 'data')
	#plt.show()

	coeff, var_matrix = curve_fit(gauss, bincen, counts, guess, 
		bounds=([0,-np.inf, 0, 0],[np.inf,np.inf, np.inf, np.inf]))
	[A, mu, sigma, yoff] = coeff
	c_exp = gauss(bincen, A, mu, sigma, yoff)
	resid2 = sum((c_exp - counts)**2)
#print coeff
	
#	print var_matrix
	
	#display check
	#hist_fit=gauss(bincen, *coeff)
	#plt.plot(bincen, counts, label='data')
	#plt.plot(bincen, hist_fit, label='fit')
	#plt.show()
	return [A, mu, sigma, yoff, resid2]
	

def gauss(x, *p):
    A, mu, sigma, yoff = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + yoff