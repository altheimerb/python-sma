import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")))
import sma_lib.loadtrace as loadtrace
import math
import datetime
import circfit
import gaussfit
import os

#For each trace, fit to circle, generate histogram, fit to gaussian.
#for DP analysis
#assumes single color
#BDA 20180802
#under development

def get_angHist(trdir):
	print("Fitting circles, generating histograms, and fitting " + 
	"histograms on %s at %s" % (trdir, datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
	
	#parameters? may convert to xml file later?
	circfit_start = 1 #which frames to fit to circle?
	circfit_end = 2000
	circcen_constr = 3 #circle fitting constraints; values taken from igor code. x,y constraints, relative to initial guess
	circfit_middle1 = int((circfit_end-circfit_start+1)/2) #for fitting over two intervals to check for drift.
	circfit_middle2 = circfit_middle1#circfit_end - int((circfit_end-circfit_start+1)/4)

	
	hist_start = circfit_start #which frames to include in histogram
	hist_end = circfit_end
	hist_len = hist_end-hist_start+1
	radguess = 0.4 #initial guess for circle radius. currently not used
	binsize = 10 #bin size for histogram, in deg
	allow_yoff = 0 #set to 1 to allow non-zero y-offset during gaussian fit of histogram; 0 for 'pure' gaussian, no offset
	
	if(allow_yoff==1): print("gaussian with y offset")
	if(allow_yoff==0): print("gaussian without y offset")
	
	#open peak list
	listptr = open(trdir+'\\trlist.txt','r')
	peaklist = []
	for line in listptr:
		peaklist.append(line[0:5])
	n_tr = len(peaklist)
	
	#for testing
	#n_tr = 5
	#fixme
	print("Number of traces: "+str(n_tr))
	
	
	#other stuff that doesn't need to get re-calculated
	n_bins = 360 / binsize
	bins1 = np.zeros(n_bins+1) #need the right edge also
	for i in range(0,n_bins+1):
		bins1[i]= -180 + i*10
	#print(bins1)	
	
	#create arrays to fill in with each trace's info
	circres = np.zeros((n_tr, 4))#(x,y,radius,avg sq residual)
	circresp1 = np.zeros((n_tr, 4))#(x,y,radius,avg sq residual) #fit using first half of pts
	circresp2 = np.zeros((n_tr, 4))#(x,y,radius,avg sq residual) #fit using second half of pts
	anghists_raw = np.zeros((n_tr,n_bins)) #not re-centered
	anghists = np.zeros((n_tr,n_bins)) #recenter: max bin at angle = 0
	gaussfits = np.zeros((n_tr,5)) #gaussian fits - amplitude, mean, width, y_offset, residual
	
	#go through each trace.
	for tr in range(0,n_tr):
		if 1:# tr%100 ==0:
			print 'working on trace ' + str(tr)
		#load trace
		trace = loadtrace.load_trace(trdir+'//'+peaklist[tr]+'.tr')
		
		#replace 0's -failed fits - with nan
		trace['xx'][trace['xx']==0] = np.nan
		trace['yy'][trace['yy']==0] = np.nan
		
		#extract x, y traces separately, over the desired range
		xcur = trace['xx'][circfit_start:circfit_end+1]
		ycur = trace['yy'][circfit_start:circfit_end+1]
		#fixme: what happens if desired range is longer than trace?

		#fit to a circle
		xav = np.nanmean(xcur) #mean, ignoring Nan values
		yav=np.nanmean(ycur)
		finiteMask = np.isfinite(xcur) #to remove Nans and infs before fitting - otherwise, scipy least squares fails
		#also generate arrays with halves of the data for drift check
		xcurp1 = trace['xx'][circfit_start:circfit_middle1]
		ycurp1 = trace['yy'][circfit_start:circfit_middle1]
		xcurp2 = trace['xx'][circfit_middle2:circfit_end]
		ycurp2 = trace['yy'][circfit_middle2:circfit_end]
		finiteMaskp1 = np.isfinite(xcurp1)
		finiteMaskp2 = np.isfinite(xcurp2)
		
		
		#print(finiteMask)
		#[circfit_x, circfit_y, circfit_rad, circfit_resid] = circfit.cfit(xcur, ycur, trace['xpos'], trace['ypos'],radguess, 
		#bounds=([trace['xpos']-3,trace['ypos']-3],[trace['xpos']+3,trace['ypos']+3]))
		try: 
			[circfit_x, circfit_y, circfit_rad, circfit_resid] = circfit.cfit(xcur[finiteMask], ycur[finiteMask], xav, yav,radguess, 
				bounds=([xav-3,yav-3],[xav+3,yav+3]))
			circres[tr,:] = [circfit_x, circfit_y, circfit_rad, circfit_resid]
			#fit to circle over half intervals to check for drift
			[circfit_xp1, circfit_yp1, circfit_radp1, circfit_residp1] = circfit.cfit(xcurp1[finiteMaskp1], ycurp1[finiteMaskp1], xav, yav,radguess, 
				bounds=([xav-3,yav-3],[xav+3,yav+3]))
			circresp1[tr,:] = [circfit_xp1, circfit_yp1, circfit_radp1, circfit_residp1]
			[circfit_xp2, circfit_yp2, circfit_radp2, circfit_residp2] = circfit.cfit(xcurp2[finiteMaskp2], ycurp2[finiteMaskp2], xav, yav,radguess, 
				bounds=([xav-3,yav-3],[xav+3,yav+3]))
			circresp2[tr,:] = [circfit_xp2, circfit_yp2, circfit_radp2, circfit_residp2]
		except(ValueError):
			print("Trace %d: invalid input to circle fits; probably trace is all NaNs" %tr)
			circres[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circresp1[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circresp2[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circfit_x = 0
			circfit_y = 0 #keeps downstream code happy.
		except(ZeroDivisionError): #this occurs when trace or half-trace is all NaNs	
			print("Trace %d: divide by zero error. Likely all NaNs or half of trace NaNs" %tr)
			circres[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circresp1[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circresp2[tr,:] = [np.nan, np.nan, np.nan, np.nan]
			circfit_x = 0
			circfit_y = 0 #keeps downstream code happy
		
		#calculate angular information from xcur, ycur
		#angcur = np.zeros(hist_len)
		poscomp=np.zeros(hist_len)
		poscomp=(xcur[:]-circfit_x)-1j*(ycur[:]-circfit_y) #convert to complex value after centering at (0,0)
		#the minus sign is for consistency with igor code; convention only

		angcur = np.angle(poscomp,deg=1)

		#generate histogram from angcur
		[anghistcur, bins1] = np.histogram(angcur,bins=bins1, density=True) #normalized, assuming equal sized bins
		anghistcur *= binsize #np histogram normalizes for integration, not summation over bins
		anghists_raw[tr,:] = anghistcur
	
		#recenter around peak. careful b/c 'mod' is always positive in python
		peakpos = np.argmax(anghistcur)
		angcur_sh = ((angcur - bins1[peakpos]+180) % 360) - 180 
		[anghistcur_sh, bins1] = np.histogram(angcur_sh, bins=bins1, density=True) #normalized, assuming equal sized bins
		anghistcur_sh *= binsize #np histogram normalizes for integration, not summation over bins
		anghists[tr,:] = anghistcur_sh
		#note: bins1 values are bin edges.
		
		#fit each to an offset gaussian
		if(allow_yoff==1):
			guess = [(1.0/n_bins), 0, 90, (0.1/n_bins)]
		elif(allow_yoff==0):
			guess = [(1.0/n_bins), 0, 90, -1] #-1 as guess for y_off is instruction to fix at zero.

		try:
			[A, mu, sigma, yoff, resid] = gaussfit.fitgauss(bins1,anghistcur_sh,guess)
			
			#maybe more reliable than assuming max is the center position:
			#recenter histogram using fit result, then re-fit
			angcur_sh2 = ((angcur_sh - mu+180) % 360) - 180 
			[anghistcur_sh2, bins1] = np.histogram(angcur_sh2, bins=bins1, density=True) #normalized, assuming equal sized bins
			anghistcur_sh2 *= binsize #np histogram normalizes for integration, not summation over bins
			anghists[tr,:] = anghistcur_sh2
		
			#refit
			if(allow_yoff==1):
				guess = [(1/n_bins), 0, 90, (0.1/n_bins)]
			elif(allow_yoff==0):
				guess = [(1/n_bins), 0, 90, -1] #-1 as guess for y_off is instruction to fix at zero.
			try:
				[A, mu, sigma, yoff, resid] = gaussfit.fitgauss(bins1,anghistcur_sh2,guess)
			except(ValueError):
				print("Trace %d: gaussian fit error. Likely NaNs." %tr)
				[A, mu, sigma, yoff, resid] = [np.nan, np.nan, np.nan, np.nan, np.nan]
				#anghistcur_sh2 = np.zeros((n_bins))
				#anghistcur_sh2 = np.nan
				anghists[tr,:] = anghistcur_sh2
			except(RuntimeError):
				print("fitting failure - RuntimeError")
				[A, mu, sigma, yoff, resid] = [np.nan, np.nan, np.nan, np.nan, np.nan]
				#anghistcur_sh2 = np.zeros((n_bins))
				#anghistcur_sh2 = np.nan
				anghists[tr,:] = anghistcur_sh2
		except(ValueError):
			print("Trace %d: gaussian fit error. Likely NaNs." %tr)
			[A, mu, sigma, yoff, resid] = [np.nan, np.nan, np.nan, np.nan, np.nan]
			#anghistcur_sh2 = np.zeros((n_bins))
			anghistcur_sh2 = anghistcur_sh #might as well output less optimally shifted histogram
			#anghistcur_sh2 = np.nan
			anghists[tr,:] = anghistcur_sh2
		except(RuntimeError):
			print("fitting failure - RuntimeError")
			[A, mu, sigma, yoff, resid] = [np.nan, np.nan, np.nan, np.nan, np.nan]
			#anghistcur_sh2 = np.zeros((n_bins))
			anghistcur_sh2 = anghistcur_sh #might as well output less optimally shifted histogram
			anghists[tr,:] = anghistcur_sh2

		
		
		gaussfits[tr,:] = [A, mu, sigma, yoff, resid]
		gf = np.array([A, mu, sigma, yoff, resid])
		
		#save histograms, as individual files
		d=os.path.dirname(trdir+"\\anghists\\")	
		if not os.path.exists(d):
			os.makedirs(d)		
		
		trpt =open(trdir+'\\anghists\\'+"anghist_raw"+str(int(tr)) + '.hi','wb')
		#header: how many bins, as uint32
		infsave = np.array([n_bins])
		infsave = infsave.astype('uint32')
		infsave.tofile(trpt)
		b = np.copy(bins1)
		b=b.astype('float32')
		h = np.copy(anghistcur)
		h = h.astype('float32')
		#main text: bins and counts, as float32
		b.tofile(trpt)
		h.tofile(trpt)
		trpt.close()
		
		trpt =open(trdir+'\\anghists\\'+"anghist"+str(int(tr)) + '.hi','wb')
		#header: how many bins, as uint32
		infsave = np.array([n_bins])
		infsave = infsave.astype('uint32')
		infsave.tofile(trpt)
		b = np.copy(bins1)
		b=b.astype('float32')
		h = np.copy(anghistcur_sh2)
		h = h.astype('float32')
		#main text: bins and counts, as float32
		b.tofile(trpt)
		h.tofile(trpt)
		trpt.close()
		
		#save gaussian fits, as individual files
		trpt = open(trdir+'\\anghists\\'+"gaussfit"+str(int(tr)) + '.fit','wb')
		f=np.copy(gf)
		f=f.astype('float32')
		f.tofile(trpt)
		trpt.close()
	
	#output circle fits as a text file (similar to pks3d file), for easy reference
	c_tosave = np.zeros((n_tr,5))
	c_tosave[:,1:5] = circres[0:n_tr,:]
	for i in range(0,n_tr):
		c_tosave[i,0] =i 
	#c_tosave = np.transpose(c_tosave)
	format = ['%- i','%-.6f','%-.6f','%-.6f','%-.6f']
	np.savetxt(trdir +'\\circlefits.txt',c_tosave,fmt=format,delimiter='\t')
	
	#output as binary file, for easy import
	trpt = open(trdir +'\\circlefits.b', 'wb')
	c = np.copy(circres)
	c=c.astype('float32')
	#print c
	c0=c[:,0]
	c1=c[:,1]
	c2=c[:,2]
	c3=c[:,3]
	c0.tofile(trpt)
	c1.tofile(trpt)
	c2.tofile(trpt)
	c3.tofile(trpt)
	trpt.close()
	
	#save circle fits on partial trace - two files p1 and p2
	#p1:
	c_tosave = np.zeros((n_tr,5))
	c_tosave[:,1:5] = circresp1[0:n_tr,:]
	for i in range(0,n_tr):
		c_tosave[i,0] =i 
	#c_tosave = np.transpose(c_tosave)
	format = ['%- i','%-.6f','%-.6f','%-.6f','%-.6f']
	np.savetxt(trdir +'\\circlefitsp1.txt',c_tosave,fmt=format,delimiter='\t')
	
	#output as binary file, for easy import
	trpt = open(trdir +'\\circlefitsp1.b', 'wb')
	c = np.copy(circresp1)
	c=c.astype('float32')
	#print c
	c0=c[:,0]
	c1=c[:,1]
	c2=c[:,2]
	c3=c[:,3]
	c0.tofile(trpt)
	c1.tofile(trpt)
	c2.tofile(trpt)
	c3.tofile(trpt)
	trpt.close()
	#p2:
	c_tosave = np.zeros((n_tr,5))
	c_tosave[:,1:5] = circresp2[0:n_tr,:]
	for i in range(0,n_tr):
		c_tosave[i,0] =i 
	#c_tosave = np.transpose(c_tosave)
	format = ['%- i','%-.6f','%-.6f','%-.6f','%-.6f']
	np.savetxt(trdir +'\\circlefitsp2.txt',c_tosave,fmt=format,delimiter='\t')
	
	#output as binary file, for easy import
	trpt = open(trdir +'\\circlefitsp2.b', 'wb')
	c = np.copy(circresp2)
	c=c.astype('float32')
	#print c
	c0=c[:,0]
	c1=c[:,1]
	c2=c[:,2]
	c3=c[:,3]
	c0.tofile(trpt)
	c1.tofile(trpt)
	c2.tofile(trpt)
	c3.tofile(trpt)
	trpt.close()
	
	
	
	#output .info file on settings
	trpt = open(trdir +'\\anghistSettings.info', 'wb')
	a=np.array([circfit_start])
	a=a.astype('uint32')
	a.tofile(trpt)
	a=np.array([circfit_end])
	a=a.astype('uint32')
	a.tofile(trpt)
	a=np.array([hist_start])
	a=a.astype('uint32')
	a.tofile(trpt)
	a=np.array([hist_end])
	a=a.astype('uint32')
	a.tofile(trpt)
	a=np.array([binsize])
	a=a.astype('uint32')
	a.tofile(trpt)
	trpt.close()
	
	
	print("get_angHist done at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
	
	
if __name__ == "__main__":
	# check input
	if(len(sys.argv)==2):
		trdir = sys.argv[1]
	else:
		print "usage: <trdir>"
		exit()
	if '.trdir' in trdir:
		trdir = trdir[:-6]
	
	get_angHist(trdir)