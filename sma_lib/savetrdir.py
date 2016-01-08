#save trace directory files

import numpy as np
import datetime
import os


def save_trdir(filename, par, peaks, time_tr, crds_tr,done):
	if done == 0:
		print 'saving PARTIAL .trdir at ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	elif done == 1:
		print 'saving complete .trdir at ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	no_peaks = peaks.shape[0]
	peaks_dim = peaks.shape[1] #useful since size of peaks depends on analysis type
	print "number of peaks: " + str(no_peaks)
	trd = par.file + 'trdir'
	trlist = trd + '\\trlist.txt'
	d=os.path.dirname(trlist)	
	if not os.path.exists(d):
		os.makedirs(d)
	#first, save a list of the analyzed traces.
	trlptr = open(trlist,'w')
	for tr in range(0,no_peaks):
		trs = str(int(peaks[tr,0]))
		trs = trs.zfill(5) #pad with zeros to get 5 digits
		trlptr.write('%s\n' %trs)
	trlptr.close()
	#save an unformatted file with some useful overal expt info
	infoptr = open(trd+'\\analysisdetails.inf','wb')
	infsave = np.zeros(3)
	infsave1 = np.zeros(1)
	if par.ALEX4 == 0: #effective number frames.
		#infoptr.write(bytearray(long(par.apmax_fr - par.apst_fr)))
		infsave[0] = par.apmax_fr - par.apst_fr
	if par.ALEX4 == 1:
		#infoptr.write(bytearray(long((par.apmax_fr - par.apst_fr)/4)))
		infsave[0] = (par.apmax_fr - par.apst_fr)/4
	#infoptr.write(bytearray(long(par.apmax_fr - par.apst_fr)))#number camera frames
	#infoptr.write(bytearray(long(no_peaks))) #number peaks
	infsave[1] = par.apmax_fr - par.apst_fr#number camera frames
	infsave[2] = no_peaks
	#infoptr.write(bytearray(int(par.ALEX4)))
	infsave1[0] = int(par.ALEX4)
	infsave = infsave.astype('int32')
	infsave1=infsave1.astype('int16')
	#infsave = infsave.byteswap(True)
	#infsave1 = infsave1.byteswap(True)
	infsave.tofile(infoptr)
	infsave1.tofile(infoptr)
	infoptr.close()
		
		
	#then, save each trace. set up to be compatible with the old IDL output / ORBITv6 igor. (assuming single emission channel)
	for tr in range(0,no_peaks):
		trs = str(int(peaks[tr,0]))
		trs = trs.zfill(5) #pad with zeros to get 5 digits
		trpt =open(filename+'trdir\\'+trs + '.tr','wb')
		#'header' info
		first = peaks[tr,peaks_dim-2]
		last = peaks[tr,peaks_dim-1]
		trlen = last - first + 1
		if par.ALEX4 == 1:
			first = first/4
			last = last/4
			trlen = last-first+1
		infsave = peaks[tr,0]
		infsave = infsave.astype('uint32')
		infsave1 = peaks[tr,1:peaks_dim-2]
		infsave1 = infsave1.astype('float32')
		infsave2 = np.array([peaks[tr,peaks_dim-2],peaks[tr,peaks_dim-1],first,last,trlen])
		infsave2 = infsave2.astype('uint32')
		infsave.tofile(trpt)
		infsave1.tofile(trpt)
		infsave2.tofile(trpt)
		#actual trace. different for alex4 = 0 or 1
		ctimetr = time_tr[tr] #extract the single trace from the list
		ccrdstr = crds_tr[tr]
		ctimetr = ctimetr.astype('int32')
		ccrdstr=ccrdstr.astype('float32')
		if par.ALEX4 == 0:
			ctimetr.tofile(trpt)
			temp = ccrdstr[:,:,0]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,1]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,2]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,3]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,4]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,5]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,6]
			temp.tofile(trpt)
			temp = ccrdstr[:,:,7]
			temp.tofile(trpt)
		elif par.ALEX4 ==1:#save color1,then color2
		#color1
			temp = timetr[:,0:-2:2]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,0]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,1]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,2]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,3]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,4]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,5]
			temp.tofile(trpt)	
			temp = ccrdstr[:,0:-2:2,6]
			temp.tofile(trpt)
			temp = ccrdstr[:,0:-2:2,7]
			temp.tofile(trpt)	
		#color2	
			temp = timetr[:,1:-1:2]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,0]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,1]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,2]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,3]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,4]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,5]
			temp.tofile(trpt)	
			temp = ccrdstr[:,1:-1:2,6]
			temp.tofile(trpt)
			temp = ccrdstr[:,1:-1:2,7]
			temp.tofile(trpt)	
			
		trpt.close()