#input filename of pks file and parameters file. 
#return properly formatted peaks numpy array
import numpy as np
import math

def load_peaks(filename,par):
	fileptr = open(filename,'r')
	
	if par.pks_type == 1 and par.emchs ==1: #pks3d, one channel
		peaks = np.zeros((99999,5))
		peaksraw = np.zeros((99999,5)) #without frame buffering
		no_p = 0
		xpos =0
		ypos = 0
		startt = 0
		stopt = 0 
		foo=0
		
		for line in fileptr:
			#print line
			#foo, xpos,ypos,startt,stopt = line.strip().split("\t")#thisfails if importing .pks3d fromIDL
			foo,xpos,ypos,startt,stopt = line.strip().split()
			startt = float(startt)
			stopt = float(stopt)
			peaksraw[no_p,:] = [float(foo),float(xpos), float(ypos), float(startt), float(stopt)]
			startt = startt - par.buffer_fr
			stopt = stopt + par.buffer_fr
			if startt < par.apst_fr:
				startt = par.apst_fr
			if stopt > par.apmax_fr:
				stopt = par.apmax_fr
			if par.ALEX4 ==1: #event should be over full sets of 4 camera frames
				startt = math.floor(float(startt)/4)*4
				stopt = (math.floor(float(stopt)/4)*4) -1 
			#what if the event occurs entirely outside of the frames to be analyzed? then apparent event length now < 0
			if (stopt-startt)<0.0:
				stopt = startt
			peaks[no_p,:] = [float(foo),float(xpos), float(ypos), float(startt), float(stopt)]
			#print foo, xpos, ypos, startt,stopt
			no_p += 1;
		peaks = peaks[0:no_p,:]
		if no_p > 99999:
			print "too many traces to handle in one trdir!"
		
		
#Allow using only a subset of the peaks
		if par.an_sub == 1: #continuous subset
			peaks = peaks[par.subcont1:par.subcont2+1,:]
			peaksraw = peaksraw[par.subcont1:par.subcont2+1,:]

		elif par.an_sub ==2: #every subi'th peak.
			peaks = peaks[par.suboff:no_p+1:par.subi,:]
			peaksraw = peaksraw[par.suboff:no_p+1:par.subi,:]

		no_p = peaks.shape[0]
		
#allow cutting out events which start too early 
		peaks2 = np.zeros((no_p,5))
		no_p2 = 0
		for i in range(0,no_p):
			if peaksraw[i,3] >=par.startcut:
				peaks2[no_p2,:] = peaks[i,:]
				no_p2 +=1
		peaks2 = peaks2[0:no_p2,:]
		peaks = peaks2
		no_p = no_p2
#if ignore_bounds = 1, reset the start and stop times
		if par.ignore_bounds == 1:
			peaks[:,3] = 0
			peaks[:,4] = par.apmax_fr
			print 'ignoring bounds from peak picking...'
			
			

	print no_p , " events were found in file (after filtering, if applicable)"
	
	return peaks

