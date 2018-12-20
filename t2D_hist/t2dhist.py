import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")))
import sma_lib.loadtrace as loadtrace
import math
import datetime

#output 2d histogram of trace localizations. 
#Each histogram is up to frnum frames long - if trace is longer, split into multiple histograms
#output is an integer count (uint32). Can normalize later if desired.
#saved as a stack of histograms, each resolutionxresolution in size
#assumes single color
def get_2d_hist(trdir):
	print "Generating 2D histograms starting at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + 'on ' + trdir
	resolution = 30
	fov = 1 #in pixels; extends +/- fov in each direction from center
	#du = fov / resolution / 2 #size of each bin, in pixels.
	frnum = 1000.0
	
	curhist = np.zeros((resolution,resolution)) #only need one in memory at a time.
	
	#open peak list
	listptr = open(trdir+'\\trlist.txt','r')
	peaklist = []
	for line in listptr:
		peaklist.append(line[0:5])
	n_tr = len(peaklist)
	
	#open file for 2d histograms.
	histptr = open(trdir+'\\hist2d.stack','wb')
	#open file for hist_sets info
	setsptr = open(trdir+'\\histsets.list','wb')
	hist_sets_cur = np.zeros((1,3)) #hist set info for a single trace
	
	#save parameters
	parptr = open(trdir+'\\histpar.info','wb')
	inf = np.array([resolution, fov, frnum])
	inf = inf.astype('uint32')
	inf.tofile(parptr)
	
	count = 0
	#go through each trace.
	for tr in range(0,n_tr):
		if tr%1000 ==0:
			print 'working on trace ' + str(tr)
		#load trace
		trace = loadtrace.load_trace(trdir+'//'+peaklist[tr]+'.tr')
		#get tr length --> number of sets, n_sets
		n_sets = int(math.ceil(float(trace['trlen'])/frnum))
		
		#replace 0's -failed fits - with nan
		trace['xx'][trace['xx']==0] = np.nan
		trace['yy'][trace['yy']==0] = np.nan

		for s in range(0,n_sets):
			xcur = trace['xx'][int(s*frnum):int((s+1)*frnum-1)]
			ycur = trace['yy'][int(s*frnum):int((s+1)*frnum-1)]
			
			
			#define mean as center of localizations
			xcen = np.nanmean(xcur)
			ycen = np.nanmean(ycur)
			
			#get 2d histogram
			if(np.isnan(xcen) or np.isnan(ycen))
				curhist = nan
			else
				curhist,x_edges,y_edges = np.histogram2d(xcur[~np.isnan(xcur)],ycur[~np.isnan(ycur)],bins=resolution,range=[[xcen-fov,xcen+fov],[ycen-fov,ycen+fov]])
				
			#save 2d histogram to open file
			curhist = curhist.astype('uint32')
			curhist.tofile(histptr)
			#save hist_sets info
			hist_sets_cur[0:3]=[tr,s*frnum,(s+1)*frnum]
			hist_sets_cur = hist_sets_cur.astype('uint32')
			hist_sets_cur.tofile(setsptr)
			count +=1
	c = np.array([count])
	c = c.astype('uint32')
	c.tofile(parptr)
	
	setsptr.close()
	histptr.close()	
	parptr.close()

	print "Done generating 2D histograms at " + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + 'on ' + trdir		
	
	
if __name__ == "__main__":
	# check input
	if(len(sys.argv)==2):
		trdir = sys.argv[1]
	else:
		print "usage: <trdir>"
		exit()
	if '.trdir' in trdir:
		trdir = trdir[:-6]
	
	get_2d_hist(trdir)