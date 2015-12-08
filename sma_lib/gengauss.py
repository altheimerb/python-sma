#generate gaussian peak for data weighting
import numpy as np
import math
def gen_gauss():
	g_peaks = np.zeros((2,2,7,7))
	for k in range(0,2):
		for l in range(0,2):
			offx = -0.5*float(k)
			offy = -0.5*float(l)
			for i in range(0,7):
				for j in range(0,7):
					dist = 0.4 * ((float(i)-3.0+offx)**2 + (float(j)-3.0+offy)**2)
					g_peaks[k,l,i,j] = 2.0*math.exp(-dist)
	print 'gaussian masking is not tested'				
	return g_peaks