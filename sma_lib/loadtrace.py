import numpy as np

def load_trace(filename):
	fileptr = open(filename,'rb')
	#read in header info
	[trnum]  = np.fromfile(fileptr, dtype = 'uint32',count=1)
	[xpos,ypos] = np.fromfile(fileptr,dtype ='float32',count=2)
	[cstartfr,cstopfr,estartfr,estopfr,trlen] = np.fromfile(fileptr,dtype='uint32',count=5)
	
	c0int = np.fromfile(fileptr,dtype='int32', count = trlen)
	xx = np.fromfile(fileptr,dtype = 'float32', count = trlen)
	yy = np.fromfile(fileptr,dtype = 'float32', count = trlen)
	#there's more to read in, but I don't need it anytime soon.
	
	fileptr.close
	
	#store the data as a dict
	tr = {'trnum':trnum, 'xpos': xpos, 'ypos':ypos,'trlen':trlen,'c0int':c0int,'xx':xx,'yy':yy}
	
	return tr