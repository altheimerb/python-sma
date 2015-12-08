#calculate background. following method in the IDL code
#Ben Oct 2015
import numpy as np
import scipy as sp
import scipy.ndimage

#img is a np array. bksize is # pixels.
#works as is
def sm_bkgr(img,bksize):
	dimx = img.shape[1]
	dimy = img.shape[0]
	fr_bk = np.zeros((dimy/bksize,dimx/bksize))
	for i in range(bksize/2,dimy,bksize):
		for j in range(bksize/2,dimx,bksize):
			tempimg = img[(i-bksize/2):(i+bksize/2),(j-bksize/2):(j+bksize/2)]
			#print tempimg
			np.reshape(tempimg,(bksize**2),1)
			fr_bk[(i-bksize/2)/bksize,(j-bksize/2)/bksize] = np.median(tempimg)
			#print np.median(tempimg)
	#print fr_bk
	#print dimy/bksize
	fr_bk = np.repeat(fr_bk,bksize,axis = 0) #equiv to 'rebin' in IDL	
	fr_bk = np.repeat(fr_bk,bksize,axis=1)
	#idl: boxcar average, 20 frames, edge truncate to smooth this array
	fr_bk = scipy.ndimage.filters.uniform_filter(fr_bk, (20,20))
	#could also use a gaussian filter instead:
	#fr_bk = scipy.ndimage.filters.gaussian_filter(fr_bk,(20,20)) #smooth out background
	return fr_bk