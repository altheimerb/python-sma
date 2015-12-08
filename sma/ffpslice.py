#given an image and settings, find fluorescent peaks. 
#peak finding based on hazen's IDL code ffp_dax2.
import numpy as np
import math
import sma_lib.smbkgr as smbkgr
from PIL import Image
import os 
#import sa_library.arraytoimage
import sma_lib.np2ctimg as np2ctimg
#import matplotlib.pyplot as plt
from matplotlib import cm #colormaps
#image could be raw frame, median filtered frameset (ORBIT), mean of some frames
#may be entire frame or one channel of a multi-channel expt. 
#may also be the sum of a multi-channel expt's channels.
#circle for checking peaks below
circle = np.zeros((11,11))
circle[:,0] = [0,0,0,0,0,0,0,0,0,0,0]
circle[:,1] = [0,0,0,1,1,1,1,1,0,0,0]
circle[:,2] = [0,0,1,0,0,0,0,0,1,0,0]
circle[:,3] = [0,1,0,0,0,0,0,0,0,1,0]
circle[:,4] = [0,1,0,0,0,0,0,0,0,1,0]
circle[:,5] = [0,1,0,0,0,0,0,0,0,1,0]
circle[:,6] = [0,1,0,0,0,0,0,0,0,1,0]
circle[:,7] = [0,1,0,0,0,0,0,0,0,1,0]
circle[:,8] = [0,0,1,0,0,0,0,0,1,0,0]
circle[:,9] = [0,0,0,1,1,1,1,1,0,0,0]
circle[:,10]= [0,0,0,0,0,0,0,0,0,0,0]

#set up color table. there is probably a better way to do this.
ct = np.zeros((4,256))
for i in range(0,256):
	for j in range(0,3):
		#ct[j,i]=cm.jet(i)[j] * 255
		#ct[j,i] = cm.nipy_spectral(i)[j]*255
		ct[j,i] = cm.CMRmap(i)[j]*255 #THIS ONE LOOKS GOOD.

def ffp_slice(img,frnum,par):
	dimx = img.shape[1]
	dimy = img.shape[0] #not necessarily the same as the camera's whole frame
	#bordersize : how much of the edge to ignore
	bordersize = 8
	#first, set all values < threshold to zero
	med = np.median(np.reshape(img,(dimx*dimy,1)))
	thresh = par.std + med
	#or, for keep, threshk
	threshk = par.std / par.keepratio  + med
	img_g = np.copy(img) #for finding good peaks
	img_k = np.copy(img) #for finding keep peaks
	img_disp = np.copy(img) #for display
	img_disp = np.rint(img_disp) #to integer values
	img_disp = img_disp.astype(int)

#gaussian peak pattern for position refinement
	g_peaks = np.zeros((3,3,7,7))
	for k in range(0,2):
		for l in range(0,2):
			offx = -0.5*float(l)
			offy = -0.5*float(k)
			for i in range(0,7):
				for j in range(0,7):
					dist = 0.4 * ((float(j) - 3.0 + offx)**2 + (float(i) - 3.0 + offy)**2)
					g_peaks[k,l,i,j] = 2.0*math.exp(-dist)
	
	#remove all pixels below threshold
	for i in range(0,dimx):
		for j in range(0,dimy):
			if(img_g[j,i] < thresh): img_g[j,i] = 0
			if(img_k[j,i]<threshk): img_k[j,i] = 0
			
	good = np.zeros((2,100000))
	keep = np.zeros((2,100000))
	no_good = 0
	no_keep = 0
	
	for i in range(bordersize,dimx - bordersize+1):
		for j in range(bordersize,dimy-bordersize+1):
			if img_k[j,i] > 0:
				gp = 0 #flag for peak is good, not just keep
				if img_g[j,i] > 0 : gp = 1 #passes std test; check around peak below
				#find the nearest maximum
				foob = img_k[j-3:j+4,i-3:i+4]
				y,x = np.unravel_index(foob.argmax(),foob.shape) #make sure this works as expected
				z = foob[y,x]
				y-=3
				x-=3
				#analyze only peaks in the current column
				if (y ==0 and x ==0):
					y = y+j
					x = x + i
					
					#check if it is a good peak - surrounding points near enough zero intensity
					quality = 1
					for k in range(-5,6):
						for l in range(-5,6):
							if circle[k+5,l+5] > 0 :
								if img[j+k,i+l] > (med+0.5*float(z)): gp= 0 #peak isn't good.
								#less picky for keep:
								if img[j+k,i+l] > (med+0.5*float(z)*par.keepratio) : quality = 0
					
					if quality == 1 :
						#refine peak position to nearest half pixels
						#print "pre refine (x,y):"
						#print x,y
						curbest = 10000
						best_x = 1
						best_y = 1
						for k in range(0,3):
							for l in range(0,3):
								#how well is the data described by a gaussian at this position. faster than fitting to refine position
								diff = np.sum(np.absolute(float(z) *g_peaks[k,l,:,:] - img[y-3:y+4,x-3:x+4]))
								if diff < curbest:
									best_x = l
									best_y = k 
									curbest = diff
						flt_x = float(x) - 0.5*float(best_x-1)
						flt_y = float(y) - 0.5*float(best_y-1)
						#print 'post refine x,y'
						#print flt_x,flt_y
						#now, add peak to list if appropriate
						if gp == 1:
							good[0,no_good] = flt_x
							good[1,no_good] = flt_y
							no_good = no_good + 1
						#at this point, any peak should be added to keep
						keep[0,no_keep] = flt_x
						keep[1,no_keep] = flt_y
						no_keep = no_keep + 1
						#'draw' peak where it was found, after refinement
						if gp == 1:
							for k in range(-5,6):
								for l in range(-5,6) :
									if circle[k+5,l+5] == 1:
										img_disp[flt_y+k,flt_x+l] = 100
						elif par.keeptype == 0: #only display these results if they are being used
								for l in range(-5,6) :
									if circle[k+5,l+5] == 1:
										img_disp[flt_y+k,flt_x+l] = 60 #different color for circle
						
	
	#FIXME: while this does overwrite the existing .bmp, the date/time stamp isn't updated. (at least on my laptop)
	#even if the original file is deleted. can fix by deleting the whole folder before rerunning analysis
	#im = Image.fromarray(img_disp.astype(np.uint8))# use mode = 'RGB' after converting grayscale to rainbow values in RGB
	#im.show()
	imname = par.file + "_peakfind\\good_fr"+ str(frnum)#+'.bmp'
	d=os.path.dirname(imname)
	if not os.path.exists(d):
		os.makedirs(d)
	np2ctimg.saveimg(img_disp,imname, colortable = ct, autoscale = False)
	
	print "no_good = %d. no_keep = %d" %(no_good,no_keep)
	#put good and keep into a list to return together
	peaks = [good[:,0:no_good],keep[:,0:no_keep]]
	
	return peaks
	

#function which uses a different type of threshold to decide whether items on the active list should be kept. uses intensity only.
#prevlist: peak positions to check.
def ffp_keep(img,frnum,par,prevlist):
	no_act = prevlist.shape[1]
	#print "no active in ffp_keep: %i" %no_act
	keeplist = np.zeros(no_act)
	keep_num = 0
	dimx = img.shape[1]
	dimy = img.shape[0]
	#generate background. note img is assumed to already be background subtracted.
	#instead, want to compare signal / background to decide whether to keep
	fr_bk = smbkgr.sm_bkgr(img,par.bksize)
	
	for i in range(0,no_act):
		if (fr_bk[prevlist[1,i],prevlist[0,i]] * par.int_thr < img[prevlist[1,i],prevlist[0,i]]): keeplist[i] = 1
		
	#display result
	img_dispk = np.copy(img)
	for i in range(0,no_act):
		if keeplist[i] == 0:
			for k in range(-5,6):
				for l in range(-5,6):
					if circle[k+5,l+5] == 1:
						img_dispk[prevlist[1,i]+k,prevlist[0,i]+l] = 60
		else:
			for k in range(-5,6):
				for l in range(-5,6):
					if circle[k+5,l+5] == 1:
						img_dispk[prevlist[1,i]+k,prevlist[0,i]+l] = 100
	
	#save img_dispk
	#im = Image.fromarray(img_dispk.astype(np.uint8))
	#im.show()
	imname = par.file + "_peakfind\\keep_fr"+ str(frnum)#+'.bmp'
	d=os.path.dirname(imname)
	if not os.path.exists(d):
		os.makedirs(d)
	np2ctimg.saveimg(img_dispk,imname, colortable = ct,autoscale = False)
	
	print "based on intensity, keeping %d out of %d spots" %(np.sum(keeplist),no_act)
	return keeplist