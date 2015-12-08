import re 
import math

#finish putting values in parameter object par; adjust as needed based on the settings
#works as is; will need to adjust for additional types of analysis
#antype = analysis type string. ignore when it doesn't matter - let it make the extra changes with no effect
def fix_par(par,filename,antype):
	par.file = filename
	
	#parameters.py only sets data type correctly for some types of data. fix longs here
	par.frameset = long(par.frameset)
	par.max_len = long(par.max_len)
	par.length_thr = long(par.length_thr)
	
	if par.hal_info : #read in values from .setup file.
		#define what we're looking for in the setup file
		#based on https://github.com/ZhuangLab/storm-analysis/blob/master/sa_library/datareader.py
		length_re = re.compile(r'number of frames = ([\d]+)')
		xs_re = re.compile(r'x_start = ([\d]+)')
		ys_re = re.compile(r'y_start = ([\d]+)')
		xstop_re = re.compile(r'x_end = ([\d]+)')
		ystop_re = re.compile(r'y_end = ([\d]+)')
		endian_re = re.compile(r' (big|little) endian')
		
		filename = filename + '.inf'
		inf_file = open(filename,"r")
		while 1 : 
			line = inf_file.readline()
			if not line: break
			
			if par.max_frame == -1: #only use length info from file if it isn't already specified!
				m = length_re.match(line) #check if the line has length info
				if m :
					par.max_frame = long(m.group(1)) -1
			if par.apmax_fr == -1:
				m = length_re.match(line)
				if m:
					par.apmax_fr = long(m.group(1)) - 1
					
			m = xs_re.match(line) 
			if m :
				par.xpix_start = long(m.group(1))
				
			m = ys_re.match(line) 
			if m :
				par.ypix_start = long(m.group(1))
				
			m = xstop_re.match(line) 
			if m :
				par.xpix_stop = long(m.group(1))
				
			m = ystop_re.match(line) 
			if m :
				par.ypix_stop = long(m.group(1))
				
			m = endian_re.search(line)
			if m:
				if m.group(1) =="big":
					par.endian = 1
				else:
					par.endian = 0
		inf_file.close()
		
	
	if par.start_frame < 0: #in case anyone uses -1
		par.start_frame = 0
				
	if par.ALEX4 == 1 : #need all frame numbers to be a multiple of four
		par.start_frame = long(math.floor(float(par.start_frame)/4) * 4)
		par.max_frame = long(math.floor((float(par.max_frame)+1)/4)*4)-1
		par.frameset = long(math.floor(float(par.frameset)/4)*4)
		par.apst_fr = long(math.floor(float(par.apst_fr)/4) * 4)
		par.apmax_fr = long(math.floor((float(par.apmax_fr)+1)/4)*4)-1
		
		
	if par.start_frame !=0 and antype == "ffpdax":
		print "Caution: Initial %d frames are being skipped" %par.start_frame
	if par.apst_fr !=0 and antype == 'apdax':
		print "Caution: Initial %d frames are being skipped" %par.apst_fr

	par.dimx = par.xpix_stop - par.xpix_start + 1
	par.dimy = par.ypix_stop - par.ypix_start + 1

	if par.dimx % par.bksize !=0 or par.dimy % par.bksize !=0:
		print "bksize needs to be a divisor of both dimensions!"
		
	print par.det_disp_set
	if(par.det_disp_set==1):
		print "Not set up for det_disp_set == True yet!!!"
	
	
	return par