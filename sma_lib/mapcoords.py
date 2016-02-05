#input x,y and mapping. output x',y' in the mapped region
#does not apply the large shift - eg half fov shift
#this is based on IDL mapping type. careful!
import numpy as np
#not tested
def map_coords(x,y,P,Q):
	deg = P.shape[0]
	newx = 0
	newy = 0
	for i in range(0,deg):
		for j in range (0,deg):
			newx += P(i,j) * (x**i) * (y**j)
			newy += Q(i,j) * (x**i) * (y**j)
	result = np.array([newx,newy])
	return result
	