#ffpdax: primary script for finding fluorescence peaks in .dax files
#for single molecule analyis

#mapping essentially works now (2channel); but I haven't verified I'm reading in the mapping info right
#ie - are we doing the right mapping or a different mapping?

#format from user: ffpdax filename xmlfile
import sys
import numpy as np
import sa_library.parameters as params
import sma_lib.fixpar as fixpar
import math
import sma_lib.loadframe as loadframe
import sma_lib.smbkgr as smbkgr
import ffpslice
from PIL import Image
import sma_lib.writexml as writexml
import matplotlib.pyplot as plt
import datetime
import os
import cv2
#from skimage import transform 
import sma_lib.mapcoords as mapcoords
codeversion = "20160215"

#12/8/15 - changed architecture to allow multithreading. 
#to run single, still call the script - see if statement at the very end
	
def ffp_dax(filename,xmlname):	
	print "ffpdax started at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " on file: " + filename
	#read in the settings in the .xml file using hazen's Parameter Class
	par = params.Parameters(xmlname+'.xml') #par is an object of type Parameters, defined in sa_library
	#to access parameters, use par.parameter name. eg par.start_frame
	#note these values can be manually changed: par.frameset = 200 replaces whatever was there.


	par = fixpar.fix_par(par,filename,'ffpdax') #function to 'fix' par by reading in needed stuff from setup file, 
	#and making other changes which might be needed based on the settings - 
	#eg round frame numbers to mutiple of 4 if alternating laser on STORM2

	print "x pixels: %d. y pixels: %d" %(par.dimx,par.dimy)

	if par.emchs == 2:	#read in mapping files
		Pr2l, Qr2l = mapcoords.readmapping(par,'r2l')
		Pl2r, Ql2r = mapcoords.readmapping(par,'l2r')
		#generate mapx and mapy for l-->r map (transforms right ch onto left)
		mapxl2r,mapyl2r = mapcoords.genmapxy(par,Pl2r,Ql2r) #does not include the dimx/2 offset
	
	fileptr = open(filename+'.dax','rb')
	#open the dax file
	#print "Filename: ", fileptr.name
	if par.d3peaks == 1:
		no_sets = long(math.floor(float(par.max_frame - par.start_frame + 1) / float(par.frameset)))
		#arrays for the pks data
		active = np.zeros((3,50000))
		complete = np.zeros((4,50000))
		#keep track of how many active, complete
		no_a = 0 #number active 
		no_com = 0#number complete
		
		#start reading in frame sets
		currset_st = long(par.start_frame) #start frame of the current set of interest
		frames = np.zeros((par.dimy,par.dimx,par.frameset)) #note: x value is column, y value is row.
		print "number sets: %i" %no_sets
		for i in range(0,no_sets): #goes from 0 to no_sets -1
			if i % 1000 == 0: print "working on %i"  %i #keep track of progress
			for k in range(0,par.frameset):
				frame = loadframe.load_frame(fileptr,currset_st+k,par)	
				#im =Image.fromarray(frame)
				#im.show()
				frame.astype(float) #change to float
				frames[:,:,k] = frame
				
			#if ALEX4 =1, pick out subset of frames to be used, based on pickcol. otherwise, use all frames
			if par.ALEX4 ==1:
				if par.pickcol == 0:
					rframes = np.zeros((par.dimy,par.dimx,(par.frameset)/4))
					for k in range(0,par.frameset,4):
						rframes[:,:,k/4] = (frames[:,:,k] + frames[:,:,k+1])/2
				elif par.pickcol ==1:
					rframes = np.zeros((par.dimy,par.dimx,par.frameset/4))
					for k in range(2,par.frameset,4):
						rframes[:,:,(k-2)/4] = (frames[:,:,k] + frames[:,:,k+1])/2
				elif par.pickcol ==2:
					rframes = np.zeros((par.dimy,par.dimx,par.frameset/2))
					for k in range(0,par.frameset,2):
						rframes[:,:,k/2] = (frames[:,:,k]+frames[:,:,k+1])/2
				else:
					print "That's not a pickcol option!"
					break	
			else:
				rframes = frames
				
					
				
			#next, median filter the frames in the frameset
			medimg = np.median(rframes, axis = 2)	
			
			#background for medimg
			fr_bk = smbkgr.sm_bkgr(medimg,par.bksize) #seems okay; check once showing images added
			medimg = medimg-fr_bk

			#deal with multiple emission channels -- either use one, or combine after mapping.
			if par.emchs ==1: 
				pass #one channel - no subset is picked out.
			elif par.emchs ==2:
				if par.pickchan ==0: #'left' channel
					medimg = medimg[:,0:par.dimx/2]
				elif par.pickchan ==1: #'right' channel; still find peaks in left coords so we get both.
					src = np.zeros((par.dimy,par.dimx/2,3))
					src[:,:,0] = medimg[:,par.dimx/2:par.dimx]
					mapyl2rCV = mapyl2r.astype(np.float32)
					mapxl2rCV = mapxl2r.astype(np.float32)
					rightmapped = cv2.remap(src,mapxl2rCV,mapyl2rCV,interpolation = cv2.INTER_CUBIC)
					medimg = rightmapped[:,:,0]
				elif par.pickchan ==2: #combine channels after warping. Using the 'left' channel coordinates
					print 'warping is lightly tested. Appears to work.'
					#medimg = medimg[:,0:par.dimx/2] + transform.warp(medimg[:,par.dimx/2:par.dimx],tformr2l)
					#medimg = medimg[:,0:par.dimx/2] + medimg[:,par.dimx/2:par.dimx] #FOR TESTING ONLY. NOT WARPED
					#templeft = medimg[:,0:par.dimx/2]
					#tempright = medimg[:,par.dimx/2:par.dimx]
					src = np.zeros((par.dimy,par.dimx/2,3))
					src[:,:,0] = medimg[:,par.dimx/2:par.dimx]
					mapyl2rCV = mapyl2r.astype(np.float32)
					mapxl2rCV = mapxl2r.astype(np.float32)
					rightmapped = cv2.remap(src,mapxl2rCV,mapyl2rCV,interpolation = cv2.INTER_CUBIC)
					rightmapped = cv2.remap(src,mapxl2rCV,mapyl2rCV,interpolation = cv2.INTER_CUBIC)
					#medimg = medimg[:,0:par.dimx/2] + cv2.remap(rightch,mapyl2r,mapxl2r,interpolation = cv2.INTER_CUBIC)
					medimg = medimg[:,0:par.dimx/2] + rightmapped[:,:,0]
					#medimg = rightmapped[:,:,0]
			else:
				print "Not set up for more than two emission channel"		
			
			
			#scale the med_img for display. for convenience, also find peaks in the scaled image
			medimg=255.0*(medimg+par.disp_off)/par.disp_fact 
			#requires scaling factors to be known ahead of time. tricky to pick out of the data -often nothing present at the begining
			
			
			#ready to go find peaks!
			sliceresult = ffpslice.ffp_slice(medimg,currset_st,par) #these are all 'raw' results in medimg frame. Wait until the very end to map, if needed
			current = sliceresult[0]
			
			#fixme: if number ch > 1, somewhere, output all channels with picked peaks circled
			#add current peaks to active unless they are already present
			no_cu = current.shape[1]
			#print "number current %d" %no_cu
			for c in range(0,no_cu):
				xc = current[0,c]
				yc = current[1,c]
				ID = 0 #has it been found?
				if xc > 0.1:	#real peaks aren't at zero
					for a in range(0,no_a):
						distance = ((xc - active[0,a])**2 + (yc - active[1,a])**2)**0.5
						#print distance
						if distance < par.dist_thr:
							ID = 1 #this peak is already in active.
							break #can stop looking for it
					if ID == 0:
						active[:,no_a] = [xc,yc,float(currset_st)]
						no_a +=1
			
			#move peaks from active to complete if they aren't found in keep
			n_t = 0 #number moved
			temp_active = np.zeros((3,50000))
			if par.keeptype == 0: 
				for a in range(0,no_a): #for each active peak, look through keep to decide whether to keep it active
					keep = sliceresult[1]
					xa = active[0,a]
					ya = active[1,a]
					ID = 0
					for c in range(0,no_cu): 
						distance = ((xa - keep[0,c])**2 + (ya - keep[1,c])**2)**0.5
						if distance < dist_thr :
							ID = 1
							break
					if ID ==1:
						temp_active[:,n_t] = active[:,a]
						n_t = n_t + 1
					else: #event over. move to complete, but only if long enough but not too long
						if((currset_st - active[2,a]) > par.length_thr) and ((currset_st - active[2,a]) < par.max_len):
							complete[0:3,no_com] = active[:,a]
							complete[3,no_com] = float(currset_st-1)
							no_com +=1
			else: #keeptype = 1
				keep = ffpslice.ffp_keep(medimg,currset_st,par,active[:,0:no_a])
				for a in range(0,no_a):
					if keep[a] == 0: #event a is done
						if((currset_st - active[2,a]) > par.length_thr) and ((currset_st - active[2,a]) < par.max_len):
							complete[0:3,no_com] = active[:,a]
							complete[3,no_com] = float(currset_st-1)
							no_com +=1
					else: #keep on active list
						temp_active[:,n_t] = active[:,a]
						n_t +=1
						
			active = temp_active
			no_a = n_t
			
			currset_st += par.frameset
			#end of analysis for this frameset
		#once we're done flipping through sets, move everything from active to complete, assuming long enough
		#print 'num active: %d' %no_a
		for a in range(0,no_a):
			if (((par.max_frame - active[2,a]) > par.length_thr) and ((par.max_frame - active[2,a]) < par.max_len)) :
				complete[0:3,no_com] = active[:,a]
				complete[3,no_com] = float(par.max_frame)
				no_com += 1
		
		
		if no_com > 0:
			print "there were %d events" %no_com
			times = complete[3,0:no_com] - complete[2,0:no_com] + 1
			print 'average event length: %f' %float(np.mean(times))
			print 'median event length: %f' %float(np.median(times))
		else:
			times = np.zeros((1,1))
		#now, if appropriate, map the values
		if par.emchs ==1:
			pass
		elif par.emchs == 2:
			#mapping - have 'left' coordinates; need 'right' channel coordinates too (**this is true even if picking on the right side, since it is mapped)
			#first, expand complete to have space for more coordinates.
			bigcomplete = np.zeros((6,no_com))
			bigcomplete[0:2,:] = complete[0:2,0:no_com]
			bigcomplete[4:6,:] = complete[2:4,0:no_com]
			for a in range(0,no_com):
				bigcomplete[2:4,a]=mapcoords.map_coords(bigcomplete[0,a],bigcomplete[1,a],Pl2r,Ql2r)
				#bigcomplete[2:4,a] = complete[0:2,a] #testing only. for null mapping
				bigcomplete[2,a] += par.dimx/2
				bigcomplete[2,a] = 0.5*round(bigcomplete[2,a]*2,0)
				bigcomplete[3,a] =0.5*round(bigcomplete[3,a]*2,0)
			complete = bigcomplete
			
			#then, need to make sure all the mapped positions are also acceptable . if not, remove that peak!
			#compare to par.frameborder
			filtcomplete = np.zeros((6,no_com))
			no_filt = 0
			for a in range(0,no_com):
				okay =1 
				if(complete[0,a]-par.frameborder < 0 or complete[0,a]+par.frameborder >par.dimx/2):
					okay = 0
				if(complete[1,a] -par.frameborder<0 or complete[1,a] +par.frameborder>par.dimy):
					okay = 0
				if(complete[2,a]-par.frameborder<par.dimx/2 or complete[2,a]+par.frameborder>par.dimx):
					okay = 0
				if(complete[3,a]-par.frameborder<0 or complete[3,a]+par.frameborder>par.dimy):
					okay = 0
				if(okay==1):
					filtcomplete[:,no_filt] = complete[:,a]
					no_filt += 1
			complete = filtcomplete
			no_com = no_filt
		
			
		#save output as text file
		c_tosave = np.zeros((complete.shape[0]+1,no_com))
		c_tosave[1:complete.shape[0]+1,:] = complete[:,0:no_com]
		for i in range(0,no_com):
			c_tosave[0,i] =i 
		c_tosave = np.transpose(c_tosave)
		format = ['%- i','%-.1f','%-.1f','%-.1f','%-.1f']
		if par.emchs ==2:
			format = ['%- i','%-.1f','%-.1f','%-.1f','%-.1f','%-.1f','%-.1f']
		np.savetxt(filename+'.pks3d',c_tosave,fmt=format,delimiter='\t')
		
		#make and save a histogram of the times
		#hbinnum = np.round((times.max()-times.min())/par.frameset)
		p95 = np.percentile(times,95)#make histogram look reasonable - don't display the very long tail
		hbinnum = np.round((p95 - times.min())/par.frameset)
		if hbinnum == 0: hbinnum =1
		hist,binedges = np.histogram(times,bins=hbinnum,range = (times.min(),p95))
		plt.plot(binedges[0:-1],hist)
		plt.xlabel('duration,frames')
		#plt.show()
		plt.savefig(filename+'durationhist.jpeg')

		
	else:
		print "not set up for this yet. use IDL code or add here"

	#save par object as an xml file. mostly the same as the input file, but some things are changed by fixpar.
	outxml = filename + "ffpdaxOUT.xml"
	writexml.write_xml(par,outxml,'ffpdax',filename,[no_com,np.mean(times),np.median(times)])

	#close the dax file
	fileptr.close()

	print "done at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

	if par.autocont==1:
		print 'automatically calling apdax'
		#print 'apdax.py'
		from apdax import ap_dax
		ap_dax(filename,xmlname)

if __name__ == "__main__":
	# check input
	if(len(sys.argv)==3):
		filename = sys.argv[1]
		xmlname = sys.argv[2]
	else:
		print "usage: <movie> <parameters>"
		exit()
	if '.dax' in filename:
		filename = filename[:-4]
	if '.xml' in xmlname:
		xmlname = xmlname[:-4]
		
	ffp_dax(filename,xmlname)