#input x,y and mapping. output x',y' in the mapped region
#does not apply the large shift - eg half fov shift
#this is based on IDL mapping type.

#TO DO: Make sure all the mapping formatting is right!!
import numpy as np

def map_coords(x,y,P,Q):
	deg = P.shape[0] #this is actually mapdeg +1
	newx = 0
	newy = 0
	for i in range(0,deg):
		for j in range (0,deg):
			newx += P[i,j] * (x**i) * (y**j)
			newy += Q[i,j] * (x**i) * (y**j)
	result = np.array([newx,newy])
	return result
	
	
#generate mapx and map for opencv mapping. using IDL style.
#note format is (y,x) to be consistent with elsewhere - input as map[y,x]
#indices are the (y,x) values in the destination image (mapped); 
#entries are the source pixel, either x or y depending on which array
def genmapxy(par,P,Q):
	#use src to decide on size of mapx, mapy
	dimx = par.dimx/par.emchs
	dimy = par.dimy
	origx = np.linspace(0,dimx-1,dimx)
	origx = np.tile(origx,dimy)
	origx = np.reshape(origx,(dimy,dimx))
	
	origy = np.linspace(0,dimy-1,dimy)
	origy = np.tile(origy,dimx)
	origy = np.reshape(origy,(dimx,dimy))
	origy = np.transpose(origy)
	#use map_coords to get values for each
	#for i in range(0,dimx):
#		for j in range(0,dimy):#
	#		mapx[j,i], mapy[j,i] = map_coords(i,j,P,Q)
	mapx = np.zeros((dimy,dimx))
	mapy = np.zeros((dimy,dimx))
	for i in range(0,par.mapdeg+1): #this is way faster than repeatedly calling map_coords above
		for j in range (0,par.mapdeg+1):
			mapx += P[i,j] * (origx**i) * (origy**j)
			mapy += Q[i,j] * (origx**i) * (origy**j)
	#format: mapx[i,j]: value is the x coordinate in the source of the point (i,j) in the image.
	
	
	
	return mapx,mapy
	#return origx,origy #testing only. gives null mapping

#read in mapping file and format
def readmapping(par,direction):
	if direction=='l2r':
		fileptr = open(par.mapfilel2r,'r')
	elif direction=='r2l':
		fileptr = open(par.mapfiler2l,'r')
	else:
		print 'not a valid mapping option'
		
	map = np.zeros(2*((par.mapdeg+1)**2))
	i=0
	for line in fileptr:
		map[i] = float(line)
		i +=1
	P,Q = mapformatIDL(map,par.mapdeg)
	return P,Q
	
	
	
#take single row/column mapping information and format for easy transform according to IDL's poly2D algorithm.
#for now, assumes input is from IDL polywarp. 
#FIXME: determine input type based on number elts; pad with zeros if higher order crossterms missing. (new python mapping code, probably)
#Confirmed that this outputs P,Q
#P[i,j] here is the same as P[j,i] in IDL code (since IDL and python index differently)
def mapformatIDL(maplist,deg):
	Plin = maplist[0:(deg+1)**2]
	Qlin = maplist[(deg+1)**2::]
	P = np.reshape(Plin,(deg+1,deg+1))
	Q = np.reshape(Qlin,(deg+1,deg+1))
	return P,Q