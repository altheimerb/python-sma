#takes peaks or pks3d file and outputs .traces file or trdir 
#(eventually. for now, limited functionality.)
#seems to work; limited to storm2; VERY limited testing so far
import sys
import numpy as np
import sa_library.parameters as params
import sma_lib.fixpar as fixpar
import math
import sma_lib.loadframe as loadframe
import sma_lib.smbkgr as smbkgr
import sma_lib.writexml as writexml
import datetime
import os
import sma_lib.gengauss as gengauss
import sma_lib.loadpeaks as loadpeaks
#import sa_library.gaussfit as gaussfit #looks like this doesn't allow fitting over all parameters - just center position
from scipy import optimize
codeversion = "20151208"
		
def ap_dax(filename,xmlname):
	print "apdax started at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " on file: " + filename
	#fixme: reduce precision on time measurement. don't need subsecond precision.

	#read in the settings in the .xml file using hazen's Parameter Class
	par = params.Parameters(xmlname+'.xml')
	par = fixpar.fix_par(par,filename,'apdax') #function to 'fix' par by reading in needed stuff from setup file, 
	#and making other changes which might be needed based on the settings - 
	#eg round frame numbers to mutiple of 4 if alternating laser on STORM2

	#for now, ignoring CMOS calibration. FIXME
	if par.hcam_cal == 1:
		print "Not set up for CMOS calibration"
	#for now, only set up for .trdir output and .pks3d input
	if par.pks_type == 0:
		print "Not set up for .pks file. Only .pks3d"
	if par.outtype ==0:
		print "Not set up for .traces file. only trdir"
		
	if par.mt == 1: #gaussian mask
		g_peaks = gengauss.gen_gauss()
		
	#read in peaks info
	if par.pks_type ==0:
		print 'not set up for pks file'
	elif par.pks_type == 1:
		pkfile = filename + '.pks3d'
		peaks = loadpeaks.load_peaks(pkfile,par) #loads and processes (adds buffer) peaks file
	else:
		print "format not recognized"
	no_peaks = peaks.shape[0]
	peaks_dim = peaks.shape[1] #useful since since of peaks depends on analysis type
		
	#make list of arrays to hold data. allows unequal lengths for different traces
	time_tr = [] #for intensities. always need.
	crds_tr = [] #for fitting. sometimes need.
	for p in range(0,no_peaks):
		curlen = int(peaks[p,peaks_dim-1]-peaks[p,peaks_dim-2])+1
		if par.ALEX4 == 1:
			curlen = curlen / 2 #note this rounds down
			#note colors are interleaved. emchs is based on emission path
		
		time_tr.append(np.zeros((par.emchs,curlen)))
		crds_tr.append(np.zeros((par.emchs,curlen,8)))
		#for holding coordinates (x,y, x_stdev, y_stdev, fitting flag, quality metric, tilt angle, fit height)
	#keep track of which position in the array next frame's info will be added to.
	addfr = np.zeros(no_peaks)
	#fixme: use peak position to decide where to center gaussian mask. see idl code.



	#start going through the frames.
	incr =1
	if par.ALEX4 ==1:
		incr = 2 #go through frames 2 at a time
	fileptr = open(filename+'.dax','rb')

	#if using constant background subtraction, determine that now:
	if par.bst==0:
		#load first 10 frames to determine background
		frs = loadframe.load_frame(fileptr,0,par)
		for i in range(1,10):
			frs += loadframe.load_frame(fileptr,i,par)
		frs = frs / 10
		fr_bkgd = smbkgr.sm_bkgr(frs,par.bksize)	
			

	#fitting stuff. fixme: choose fit type in xml file? for now, don't allow tilt
	fitfunc = lambda p,x,y: p[0] + p[1]*np.exp(-((x-p[2])/p[4])**2 -((y-p[3])/p[5])**2)
	errfunc = lambda p,x,y,data: np.ravel(fitfunc(p,x,y) - data)
	fitxval = np.zeros((par.fit_box*2+1,par.fit_box*2+1)) #independent variables for fit
	fityval = np.zeros((par.fit_box*2+1,par.fit_box*2+1))
	for i in range(0,par.fit_box*2+1): #fill arrays with row or column number
		fitxval[:,i] = i
		fityval[i,:] = i
	fitxval1D = np.ravel(fitxval)
	fityval1D = np.ravel(fityval)
	#fixme: adjust error function to allow weighting

	for i in range(par.apst_fr,par.apmax_fr+1,incr): #i is always a real camera frame number
		if i%1000 == 0:
			print "working on : " + str(i) + " " + str(par.apmax_fr) + "at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		
		if par.ALEX4 == 0:
			frame = loadframe.load_frame(fileptr,i,par)
		if par.ALEX4 ==1:
			frame1 = loadframe.load_frame(fileptr,i,par)
			frame2 = loadframe.load_frame(fileptr,i+1,par)
			frame = frame1+frame2
		#either way, now work with frame array
		
		#if set to determine background on each frame, do that
		if par.bst == 1:
			fr_bkgd = smbkgr.sm_bkgr(frame,par.bksize)
		
		#analyze each trace
		for j in range(0,no_peaks):
			if peaks[j,peaks_dim-1] >= i and peaks[j,peaks_dim-2] <= i: #this event is happening now!
				for ch in range(0,par.emchs):
					#determine intensity
					curx = peaks[j,1+2*ch]
					cury = peaks[j,2+2*ch]
					if par.mt == 0:#simple integration
						local = frame[cury-par.sibs:cury+par.sibs+1,curx-par.sibs:curx+par.sibs+1] - fr_bkgd[cury-par.sibs:cury+par.sibs+1,curx-par.sibs:curx+par.sibs+1]
						time_tr[j][ch,addfr[j]] = np.sum(local)
					elif par.mt == 1:#gaussian masking
						print 'not set up for gaussian masking yet'
				
					#if appropriate, fit
					if time_tr[j][ch,addfr[j]] > par.fit_thr:
						#use least squares optimization
						#loc = np.ravel(frame[cury-par.fit_box:cury+par.fit_box+1,curx-par.fit_box:curx+par.fit_box+1])
						loc = frame[cury-par.fit_box:cury+par.fit_box+1,curx-par.fit_box:curx+par.fit_box+1]
						[yguess,xguess] = np.unravel_index(loc.argmax(),loc.shape)
						p0 = np.array([float(fr_bkgd[cury,curx]),float(loc[yguess,xguess]),xguess,yguess,1.1,1.1])#a smarter initial guess
						loc = np.ravel(loc)
						#p0 = np.array([float(fr_bkgd[cury,curx]), float(frame[cury,curx]), par.fit_box,par.fit_box,1.0,1.0]) #initial guess

						[p1, cov_x, infodict, mesg, success] = optimize.leastsq(errfunc,p0, args = (fitxval1D,fityval1D,loc), full_output = 1,xtol=par.fitxtol)		
						
						#if fit successful and center is within box, store result
						if success>0 and success<5 and p1[2]>0 and p1[2] < (2*par.fit_box+1) and p1[3]>0 and p1[3] < (2*par.fit_box+1):
							crds_tr[j][ch,addfr[j],:] = [p1[2]+curx-par.fit_box,p1[3]+cury-par.fit_box,abs(p1[4]),abs(p1[5]),1.0,0.0,0.0,p1[1]]
							#for coordinates (x,y, x_stdev, y_stdev, fitting flag, quality metric, tilt angle, fit height)
							#fixme: not using any metric of quality.
				addfr[j] = addfr[j] + 1
	fileptr.close()
				
	print "saving data at" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	#fixme: save data here.
	if par.outtype == 0:
		print 'not set up to save .traces file'
	elif par.outtype ==1:
		print 'saving .trdir'
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
	#output xml file with actual settings in par.
	outxml = filename + "apdaxOUT.xml"
	writexml.write_xml(par,outxml,'apdax',filename,[])
			


	print "apdax done at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " on file: " + filename

if __name__ == "__main__":
	#format from user: ffpdax filename xmlfile
	# check input
	if(len(sys.argv)==3):
		filename = sys.argv[1]
		xmlname = sys.argv[2]
	else:
		print "usage: <movie> <parameters.xml>"
		exit()
		
	if '.dax' in filename:
		filename = filename[:-4]
	if '.xml' in xmlname:
		xmlname = xmlname[:-4]
		
	ap_dax(filename,xmlname)