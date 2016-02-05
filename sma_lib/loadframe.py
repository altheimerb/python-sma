#load a single frame
#Ben Oct 2015
import numpy as np

#returns frame specified by frnum in fileptr
#x value:  column. y value: row. deal with how to display later.
#IDL: (x,y).  python: (row,column) = (y,x)
def load_frame(fileptr,frnum,par):
	fileptr.seek(frnum * par.dimx*par.dimy*2) #*2 because this is looking through 8bit; data is 16bit
	frame = np.fromfile(fileptr,dtype='int16',count=par.dimx*par.dimy)
	frame = np.reshape(frame,[par.dimy,par.dimx])
	#if par.from_hal==1:
	frame = np.transpose(frame)
		#need a 270 deg rotation cw to match the results of reading array into IDL and rotate(,5) there
	frame = np.rot90(np.rot90(np.rot90(frame)))
		#frame = np.transpose(np.rot90(frame))
	#elif par.from_hal ==0:
	#	pass
	if par.endian ==1:
		frame.byteswap(True)
	
	return frame