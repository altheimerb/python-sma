#pragma rtGlobals=1		
#include <AxisSlider>
#include <All IP Procedures>
#include <Image Saver>
#include <GraphMagnifier>


//*********************************for capturing input events**************************************//
static constant hook_mousemoved=4
static constant hook_mousedown=3
static constant hook_mouseup=5
static constant hook_keyboard=11
//***************************************************************************************************//

//*********************************buffer (padding used in peak finding) for center finding algorithm**************************************//
static constant frame_padding=150
//***************************************************************************************************//


////ORBIT analysis v7. 
//For interactive analysis of rotational tracking traces
//Input: trace directory files from associated python code SMA
//function begin() defines variables, contains set up. function browselocs() sets up data windows.
// function get_2D_hist_all() sets up 2d histogram window for quickly looking through traces. function doit() runs those three in order to get analysis started.


//ORBIT analysis v7. BDA 9/30/15.
//alter to allow multiple traces to be considered in a single igor experiment.
//note: the acceptedtraces analysis part hasn't been adjusted for separating out the different input files
//Known issues
//circle fitting using user specified bounds - only set up for single color right now
//much of the analysis is only set up for single color right now.

////////////////////////////////CORE FUNCTIONS

//change path to data. useful when changing computers.
function changepath()
	wave/t data_path, find_data, trpath, dir_path
	data_path[0] = "DATAPATH:"
	dir_path[0]= "movie_0001"
	find_data[*]= data_path[p]+ dir_path[p]+"trdir:"
	print find_data
end

//run 'begin' to get started. specify what you want to analyze here.
//set settings here also.
function begin()
	variable/g num_expt= 1
	variable/g fov=1	//fov size in pixels for loc window and 2D histograms. 
	variable/g scr = 0 //=0 for data in order; =1 for random scrambling. =2 for interleaving (req. same number of traces and only 2 files)
	variable/g trackfilt = 3//=1 for no extra filtering. take a median filter of xx,yy of this many points for constraining 360deg jumps. ttc will still be unfiltered, but with 360 jumps determined using filtered data.
	
	variable/g std_window = 20 //size of sliding window for calculating std for stick detection. used by AccStick()
	variable/g std_thr = 1.3 // threshold for ratio of angular to radial fluctuations for stick detection
	variable/g std_freq = 0.05 //what fraction of frames during translocation can be below std_thr without calling a sticky trace? // set > 1 to turn off.
	
	
	//fit quality filtering
	variable/g fitqual_thr = 0.1 // radial standard deviation cutoff - filter for localization quality. set high to effectively turn off.
	
	variable/g load_2dhist = 1 //=0 to generate 2dhistograms here, =1 to load from python output.
	//THIS ONLY WORKS IF num_expt =1 and scr = 0, at least for now.
	if(load_2dhist ==1 && (num_expt != 1 || scr != 0))
		print "Incompatible settings. Changing load_2dhist to 0"
		load_2dhist = 0 
	endif
	
	///pause finding parameters
	variable/g smth_window = 111
	variable/g rate_thr = 2.5
	variable/g rev_rate_thr = -rate_thr
	variable/g start_thr = 35 	//extend pauses forward/backwards until the values are >start_thr different than the mean
	variable/g p_thr =35//merge pauses which are within p_thr of each other.
	variable/g fr_back =2   //allow fr_back - 1 frames out of 34 deg range while extending pauses forward/backward
	variable/g min_break_dur = 30
	variable/g run_dur_thr = 20
	variable/g smth_window2 = 600 //set to -1 to turn off second smoothing / pause finding step
	variable/g slip_thr = 100 //how many deg back during slip to accept as real
	
	
	variable/g maxcircfr = 3000 //don't let user try to do a fit with more than this many frames - can crash igor. this limit has not been finely tuned.
	make/o/n=(num_expt)/t data_path
	make/o/n=(num_expt)/t dir_path
	make/o/n=(num_expt)/t find_data
	make/o/n=(num_expt)/t movie_number //for naming accepted traces
	make/o/n=(num_expt)/t info_path //to trdir info file
	make/o/n=(num_expt)/t trpath //to traces list
	data_path[0]="DATAPATH:"
	//data_path[1] ="DATAPATH2" //can load more than one set of traces if desired. Analysis sections not fully set up for this though.
	dir_path[0]= "movie_0001"
	//dir_path[1]="movie_0002"
	find_data[*]= data_path[p]+ dir_path[p]+"trdir:"
	info_path[*] = data_path[p]+ dir_path[p] + "trdir:analysisdetails.inf" //don't change this line
	trpath[*] = data_path[p]+dir_path[p]+"trdir:trlist.txt"

	print find_data
	print Secs2Time(Datetime,1)
	movie_number[*] = dir_path[p]
	
	//illumination
	//options: 0= single color, 1=two color offset in time (alternating, with or without dual view), 2=simultaneous two color (dual view)
	//not all functions are set up for two color.
	variable/g ill_set = 0
	variable/g histch1 = 0 //basically, should we think of ch1 instead of ch0 as the primary channel.
	string/g ch0col = "g"
	string/g ch1col = "r"
	
	//various settings
	variable/g pc = 1	//allows different viewing settings, etc for different computers. 
	//on windows, igor graphs, etc all confined within the igor window (pc=1 better). on Mac, separate windows can be moved more freely (pc = 0 better)
	variable/g fc = 0 //set to 1 to fit data to circle. 
	//note: fc = 1 is currently required for two color circle-based realignment. using manual circle position doesn't work. change?
	//with single color, generally more efficient to set fc = 0 and manually set frames to fit over for traces worth fitting.
	if(ill_set !=0)
		fc=1 //angle data not meaningful without fit and re-alignment with alternating color on CMOS
		print "resetting fc to 1 since dual color in use"
	endif
	
	
	//read in some useful quantities from the .info file.
	make/o/n=(num_expt) expt_frCam, expt_frEff, expt_ntr
	variable expt_frCam0
	variable expt_frEff0
	variable expt_ntr0
	variable/g all_ntr=0
	variable infref
	variable i = 0
	string temp
	for(;i<num_expt;i+=1)
		Open/Z/R infref as info_path[i]
		FBinRead/F=3 infref, expt_frEff0 //after accounting for alternative laser reducing 'effective' number of frames.
		FBinRead/F=3 infref, expt_frCam0 //total number of camera frames
		FBinRead/F=3 infref, expt_ntr0
		Close infref
		expt_frCam[i] = expt_frCam0
		expt_frEff[i]=expt_frEff0
		expt_ntr[i]=expt_ntr0
		all_ntr +=expt_ntr0
		temp = dir_path[i]
		movie_number[i] = num2str(i) //redundant,but should keep old code happier. and easier to change later.
	endfor
	
	print "Number of traces: " + num2str(all_ntr)
	//read in the trace list, trlist.txt
	variable trlref
	string trtemp
	make/o/n=(all_ntr,num_expt)/t trnums = ""
	variable j = 0
	for(;j<num_expt;j+=1)
		Open/Z/R trlref as trpath[j]
		for(i=0;i<expt_ntr[j];i+=1)
			FReadLine trlref, trtemp
			trnums[i][j] = trtemp
		endfor
	endfor
	Close trlref
//note this is a string wave. use only for (a) accessing data and (b) crossreferencing to pks3d
	
	//make a full list of all traces: (trdir,filenumber)
	make/o/n=(all_ntr) dirlist
	make/o/n=(all_ntr)/t trlist
	variable added = 0
	for(i=0;i<num_expt;i+=1)
		dirlist[added,added+expt_ntr[i]]=i
		trlist[added,added+expt_ntr[i]]= trnums[p-added][i]
		added+=expt_ntr[i]
	endfor
	
	if(scr ==1)//random scrambling permutation. works on dirlist, trlist.
		make/o/n=(added) order=enoise(added)
		sort order, dirlist, trlist //sorts dirlist and trlist according to values in order (random)
	endif 
	
	if(scr==2)//interleave two expts.this overwrites drlist,trlist. useful for comparing two analyses (eg with different fit settings) on the same data.
		if(num_expt != 2)
			print "settings are inconsistent!"
			Abort
		endif
		make/o/n=(all_ntr) dirlist1,trlist1
		dirlist[0,2*all_ntr-1;2] = 0
		dirlist[1,2*all_ntr;2] =1
		trlist[0,2*all_ntr-1;2] =trnums[p/2]
		trlist[1,2*all_ntr;2] =trnums[(p-1)/2]
	endif
	//now, use dirlist,trlist to access traces.
	
	
	if(ill_set == 0)
		print "Single Color Illumination"
	endif
	if(ill_set ==1)
		print "Two Color Illumination, offset in time. With or without DualView"
	endif
	if(ill_set ==2)
		print "Simultaneous Two Color Illumination."
	endif
	
	
	//other stuff that needs to get initialized
	//waves for circle fitting
	make/o/n=(all_ntr*2,3) circpos0,circpos1 //hold all circles
	//make/o/n=1 ch1xcent,ch1ycent,ch2xcent,ch2ycent //waves for current trace's circle
	make/o/n=100 theta,c0x,c0y,c1x,c1y //waves for current traces circle
	theta = 2*pi * p /dimsize(theta,0) //used to draw circle around fit center position
	make/o/n=(all_ntr) effStart //effective start frame for each trace.
	
	//create waves and variables for angle measurements based on clicking only.
	make/n=(all_ntr)/o delta_angle1=nan,delta_angle2=nan,delta_angle3=nan
	variable/g TraceLevel=nan,TraceStepSize=nan
	
	//create trace classification wave
	make/n=(all_ntr)/o trace_eval=nan
	
	variable/G peakn=0 // peakn is the trace number within the trlist.txt file. not necessarily the same as the actual trace number 
	//because that leads to some challenges with data subsets. trnums[peakn] gives the actual trace number
	variable/g startfr = 0
	variable/g startfr_flag = 0
	variable/g view_acc = 0 //view accepted traces only. change manually after accepting traces to '1' to activate.
	variable/g view_ind = 0 //for view_acc - keep track of acc_traces index
	
	variable/g recalc = 1 //recalculate trace on refreshing movie only when something has changed...
	variable/g loadedTr = -1 //what trace is in memory right now.
	
	//		//frame number for color map
	make/N=(expt_frEff0)/O zz=p

END

//load a single trace into memory. tr number for is an entry in dirlist and trlist.
function loadtrace(tr)
	variable tr
	wave/t find_data
	nvar peakn, all_ntr, recalc, loadedTr
	
	//wave/t trnums
	wave/t trlist
	wave dirlist
	wave effStart
	peakn = tr
	
	if(peakn > all_ntr-1 || peakn < 0)
		peakn = 0
		tr = 0
		print "Attempted to access non-existent trace. peakn reset to 0"
	endif
	
	if (loadedTr != peakn)
		recalc = 1
	endif
	
	
	variable dataref
	string trs = trlist[tr]
	trs = trs[0,4] //trim away trailing space.
	//print "loading trace: " +find_data + trs + ".tr"
	Open/Z/R dataref as find_data[dirlist[tr]] + trs + ".tr"
	//print dataref
	
	//read in the header info
	variable trnum
	variable/g xpos,ypos,cstartfr,cstopfr, estartfr, estopfr
	variable/g trlen //number of frames saved for this trace
	FBinRead/F=3 dataref, trnum
	FBinRead/F=4 dataref, xpos
	FBinRead/F=4 dataref, ypos
	FBinRead/F=3 dataref, cstartfr
	FBinread/F=3 dataref, cstopfr
	FBinRead/F=3 dataref,estartfr
	FBinRead/F=3 dataref,estopfr
	FBinRead/F=3 dataref,trlen //effective frames, not camera frames.
	
	effStart[tr] = estartfr
	
	//read in the trace
	make/o/n=(trlen) c0int = NaN //color zero intensity
	FBinread/F=3 dataref, c0int
	make/o/n=(trlen) xx=NaN,yy=NaN, x_sd=NaN,y_sd=NaN,locFlag=NaN,qual=NaN,tilt=NaN,fitheight=NaN
	FBinRead/F=4 dataref, xx
	FBinRead/F=4 dataref, yy
	FBinRead/F=4 dataref, x_sd
	FBinRead/F=4 dataref, y_sd
	FBinRead/F=4 dataref, locFlag
	FBinRead/F=4 dataref, qual
	FBinRead/F=4 dataref, tilt
	FBinRead/F=4 dataref, fitheight
	
	MatrixOp /o xx = Replace(xx,0,NaN) //failed fits stored as '0' but should be NaN.
	MatrixOp/o yy = Replace(yy,0,NaN)
	MatrixOp/o c0int = Replace(c0int,0,NaN)
	MatrixOp/o x_sd = Replace(x_sd,0,NaN)
	MatrixOp/o y_sd = Replace(y_sd,0,NaN)
	MatrixOp/o qual = Replace(qual,0,NaN)
	MatrixOp/o tilt = Replace(tilt,0,NaN)
	MatrixOp/o fitheight= Replace(fitheight,0,NaN)
	
	//2 color?
	nvar ill_set
	make/o/n=(trlen) c1int = NaN
	make/o/n=(trlen) xx1=NaN,yy1=NaN, x_sd1=NaN,y_sd1=NaN,locFlag1=NaN,qual1=NaN,tilt1=NaN,fitheight1=NaN

	if(ill_set !=0)
		FBinread/F=3 dataref, c1int
		FBinRead/F=4 dataref, xx1
		FBinRead/F=4 dataref, yy1
		FBinRead/F=4 dataref, x_sd1
		FBinRead/F=4 dataref, y_sd1
		FBinRead/F=4 dataref, locFlag1
		FBinRead/F=4 dataref, qual1
		FBinRead/F=4 dataref, tilt1
		FBinRead/F=4 dataref, fitheight1
		
		MatrixOp /o xx1 = Replace(xx1,0,NaN)
		MatrixOp/o yy1 = Replace(yy1,0,NaN)
		MatrixOp/o c1int = Replace(c1int,0,NaN)
		MatrixOp/o x_sd1 = Replace(x_sd1,0,NaN)
		MatrixOp/o y_sd1 = Replace(y_sd1,0,NaN)
		MatrixOp/o qual1 = Replace(qual1,0,NaN)
		MatrixOp/o tilt1 = Replace(tilt1,0,NaN)
		MatrixOp/o fitheight1= Replace(fitheight1,0,NaN)
	endif
	
	loadedTr = peakn //update loadedTr
	
	//set scaling for all the waves from the trace. use the start frame. use effective frames, not camera frames for alternating two color
	//setscale/p x estartfr,1, c0int,xx,yy,x_sd,y_sd,locFlag,qual,tilt,fitheight
	//setscale/p x estartfr,1, c1int,xx1,yy1,x_sd1,y_sd1,locFlag1,qual1,tilt1,fitheight1
	//decided to use local frame numbers instead. less confusing; can check estartfr if I want to know where it is
	setscale/p x 0,1, c0int,xx,yy,x_sd,y_sd,locFlag,qual,tilt,fitheight
	setscale/p x 0,1, c1int,xx1,yy1,x_sd1,y_sd1,locFlag1,qual1,tilt1,fitheight1
	
	//make some other waves we'll need
	make/n=(trlen)/o ttc=NaN,ddc=NaN,asst=NaN,rsst=NaN, ttc_temp=NaN, ttc_raw=NaN, ttc_unfilt = NaN, ttc_tempfilt = NaN
	//angle,distance, angular step size, radial step size, ttc_temp for calculations, ttc_raw for no-tracking version
	//setscale/p x estartfr,1, ttc,ddc,asst,rsst, ttc_temp, ttc_raw
	setscale/p x 0,1, ttc,ddc,asst,rsst, ttc_temp, ttc_raw
	
	close dataref
END


//displays windows. set up for browsing
function browse_locs()
	//start by loading trace 0
	loadtrace(0)
	
	wave xx, yy
	wave circpos0,circpos1,theta,c0x,c0y,c1x,c1y
	nvar ill_set
	nvar fc
	nvar fov, peakn, provrad, recalc
	//variable/G fov=1	//1.5
	variable/g slower=0 //movie refresh speed
	variable/g framestep = 1 //refresh how many frames at a time in the movie?

		//current localizations
	make/N=1/O xc0=nan,yc0=nan
	make/N=1/O xc1=nan,yc1=nan
	
	
	duplicate/o xx,xxc //for movie display
	duplicate/o yy,yyc
	duplicate/o xx1,xxc1
	duplicate/o yy1,yyc1
	make/o yyd,xxd,yyd1,xxd1
	
	//define center as mean of localizations
	nvar histch1
	if(histch1 ==0)
		wavestats/m=1/q xxd
		variable/g xc=V_avg
		wavestats/m=1/q yyd
		variable/g yc=V_avg
	elseif(histch1==1)
		wavestats/m=1/q xxd1
		variable/g xc=V_avg
		wavestats/m=1/q yyd1
		variable/g yc=V_avg
	endif
	
	
		//display graphs
	execute "LocalizationsGraph()"	
	execute "DistanceGraph()"
	execute "AngleGraph()"
	execute "StepSizeTrace()"
	execute "RadStepSizeTrace()"
	execute "IntTrace()"
	
	SetAxis/W=locWindow left (yc+fov),(yc-fov)
	SetAxis/W=locWindow bottom (xc-fov),(xc+fov)
	
	//fit and adjust positions to overlap circles (2color)
	if(fc ==1 && (circpos0[peakn][0] == 0 || abs(circpos0[peakn][0] - provrad) < 10^-6)) //adjusted positions. //don't refit if you've already fit this trace
		//calculate center on one or two colors
		if(ill_set !=0) //for all other settings right now there should be two colors to fit
		 	FindCenters(1)
		else
			FindCenters(0)
		endif
		nvar recalc
		recalc = 1 //new fit --> need to recalculate trace
	endif
	
	//movie ch1 localizations such that the circles are centered at the same place (re-alignment of channels)
			xxc1 += circpos0[peakn][1] - circpos1[peakn][1] // returns NaN if ch1 doesn't exist (single color)
			yyc1 += circpos0[peakn][2] - circpos1[peakn][2]
			
			c1x = circpos1[peakn][1] + circpos1[peakn][0]*cos(theta[p])
			c1y = circpos1[peakn][2] + circpos1[peakn][0]*sin(theta[p])
			c0x = circpos0[peakn][1] + circpos0[peakn][0]*cos(theta[p])
			c0y = circpos0[peakn][2] + circpos0[peakn][0]*sin(theta[p])
	
	
	//calculate separation distance and angle
	if(recalc==1)
		update_distance_angle() //will give meaningless results if it is single color and lacks a center position.
		recalc =0
	endif
	
END

//animate movie
function update_locWindow()
	nvar peakn
	nvar loadedTr, recalc
	
	if(loadedTr != peakn) //for efficiency, skip if we're already looking at that trace.
		loadtrace(peakn) //bring peakn into memory
	endif
	DoWindow/T locWindow,"Peak: " + num2str(peakn)

	nvar xc,yc,fov
	nvar slower, framestep, startfr, startfr_flag
	nvar ill_set, fc, provrad
	nvar trlen
	wave xx,yy,xxc,yyc,xxc1,yyc1,xx1,yy1,xc0,yc0,xc1,yc1,zz
	wave circpos0,circpos1
	wave c0x,c0y,c1x,c1y, theta
	wave c0int,c1int

	variable showframes = 50 //how many trailing frames to show in the movie?
	
	duplicate/o xx,xxc //for movie display
	duplicate/o yy,yyc
	duplicate/o xx1,xxc1
	duplicate/o yy1,yyc1
	
	if(fc ==1 && (circpos0[peakn][0] == 0 || abs(circpos0[peakn][0] - provrad) < 10^-6)) //adjusted positions. //don't refit if you've already fit this trace
		//calculate center on one or two colors
		if(ill_set !=0) //for all other settings right now there should be two colors to fit
		 	FindCenters(1)
		else
			FindCenters(0)
		endif
		nvar recalc
		recalc = 1 //new fit --> need to recalculate trace
	endif
	
//move ch1 positions such that centers are overlapping - realign channels.
	if(fc ==1 || ill_set !=0)
			xxc1 += circpos0[peakn][1] - circpos1[peakn][1] // returns NaN if ch1 doesn't exist (single color)
			yyc1 += circpos0[peakn][2] - circpos1[peakn][2]
	endif	
			c1x = circpos1[peakn][1] + circpos1[peakn][0]*cos(theta[p])
			c1y = circpos1[peakn][2] + circpos1[peakn][0]*sin(theta[p])
			c0x = circpos0[peakn][1] + circpos0[peakn][0]*cos(theta[p])
			c0y = circpos0[peakn][2] + circpos0[peakn][0]*sin(theta[p])
	
	//calculate separation distance and angle
	if(recalc ==1)
		update_distance_angle() //will give meaningless results if it is single color and lacks a center position.
		recalc =0
	endif
	
	//		//define center as mean of localizations
	nvar histch1
	if(histch1 ==0)
	wavestats/m=1/q xxc
	xc=V_avg
	wavestats/m=1/q yyc
	yc=V_avg
	elseif(histch1==1)
	wavestats/m=1/q xxc1
	xc=V_avg
	wavestats/m=1/q yyc1
	yc=V_avg
	endif
	SetAxis/W=locWindow left (yc+fov),(yc-fov)
	SetAxis/W=locWindow bottom (xc-fov),(xc+fov)
	
	variable i=0
	variable zstart=i
	
	if(startfr_flag ==1)
	 	i = startfr
	 	startfr_flag = 0
	else
	 	i = 0
	endif
	
	//reset waves for animation
	//xxc = nan
	//yyc = nan
	//xxc1 = nan
	//yyc1 = nan
	//use small waves to improve speed.
	make/o/n=(showframes) xxd, yyd, xxd1, yyd1 = nan
	nvar trlen
	for(;i<trlen;i+=framestep)
		if((GetKeyState(0) & 1)!=1)	//igor 7
		//if(str2num(KeyboardState("")[0])!=1)
			if(i<showframes)
				make/o/n=(i+1) xxd,xxd1,yyd,yyd1 = nan
			endif
			xxd = nan
			xxd1 = nan
			xxd = xx[max((i-showframes+1),0)+p]
			yyd = yy[max((i-showframes+1),0)+p]
			xxd1 = xx1[max((i-showframes+1),0)+p]
			yyd1= yy1[max((i-showframes+1),0)+p]
			
			xc0 = xx[i] //showing current frame with bigger point
			yc0 = yy[i]
			xc1 = xx1[i]
			yc1 = yy1[i]
			
			//				//update colormaps
			svar ch0col,ch1col
			make/o/n=(showframes+10) zzd
			zzd = p
			if(stringmatch(ch0col, "g"))
				ModifyGraph/W=locWindow zColor(yyd)={zzd,0,50,green,0}
			else
				ModifyGraph/W=locWindow zColor(yyd)={zzd,0,50,red,0}
			endif
			if(stringmatch(ch1col,"g"))
				ModifyGraph/W=locWindow zColor(yyd1)={zzd,0,50,green,0}
			else
				ModifyGraph/W=locWindow zColor(yyd1)={zzd,0,50,red,0}
			endif
			TextBox/W=locWindow/C/N=text0 "\Z36" + num2str(i)		//frame number textbox

			doupdate
			
			if(slower)
				sleep/t slower
			endif
		else
			i=trlen//skip ahead if key is pressed
		endif
	endfor
	
	//end by plotting entire trace.
	 duplicate/o xx,xxd
	 duplicate/o yy, yyd
	 duplicate/o xx1,xxd1
	 duplicate/o yy1,yyd1
	 
	//remove markers
	yc0 = nan
	xc0 = nan
	xc1 = nan
	yc1 = nan
	
	if(stringmatch(ch0col, "g"))
		ModifyGraph/W=locWindow zColor(yyd)={zz,zstart-500,trlen+300,green,0}
	else
		ModifyGraph/W=locWindow zColor(yyd)={zz,zstart-500,trlen+300,red,0}
	endif
	if(stringmatch(ch1col,"g"))
		ModifyGraph/W=locWindow zColor(yyd1)={zz,zstart-500,trlen+300,green,0}
	else
		ModifyGraph/W=locWindow zColor(yyd1)={zz,zstart-500,trlen+300,red,0}
	endif

END

//*********************************for capturing input events**************************************//

Function KeyboardPanelHook(s)	// Capture keyboard events in main graph
	STRUCT WMWinHookStruct &s
	
//	if(s.eventCode==hook_keyboard)
	Variable MousePoint, RightClick, NewFCLevel
	variable newpeak
	nvar pc
	svar offs=root:WinGlobals:locWindow:S_TraceOffsetInfo
	variable x_off,y_off
	
	if(stringmatch(winname(0,3),"locWindow"))
		switch(s.eventCode)
			case hook_mousedown:
				RightClick= (s.eventMod & (2^4))			 // bit 4, detect right click
				if(!RightClick)
					MousePoint = str2num(stringbykey("HITPOINT",Tracefrompixel(s.mouseLoc.h,s.mouseLoc.v,""),":",";"))
					//NewFCLevel = mean(current_trace, pnt2x(current_trace, MousePoint-15), pnt2x(current_trace, MousePoint+15)) 		//average over 20 points
					
					variable/g xcenter=AxisValFromPixel("locWindow","bottom",s.mouseLoc.h)
					variable/g ycenter=AxisValFromPixel("locWindow","left",s.mouseLoc.v)
				endif
				break	
		
//			case hook_mousemoved:
//				MousePoint = str2num(stringbykey("HITPOINT",Tracefrompixel(s.mouseLoc.h,s.mouseLoc.v,""),":",";"))
//				if(numtype(MousePoint)!=2)
//					MousePoint = str2num(stringbykey("HITPOINT",Tracefrompixel(s.mouseLoc.h,s.mouseLoc.v,""),":",";"))
//					//sprintf TagText, "%.2f",cw[MousePoint]												//create tag
//					//Tag/C/N=text0/F=0/L=1 $CurrentWaveName, pnt2x(cw,MousePoint),TagText			//display tag
//				endif
//				rval= 0									// we have not taken over this event
//				break
				
			case hook_mouseup:
					//align offsets for ch0 trace, used after quickdrag. also update sep dist and angle traces with new coordinates
				x_off=str2num(stringbykey("XOFFSET",offs))
				y_off=str2num(stringbykey("YOFFSET",offs))
				ModifyGraph/W=locWindow offset(yc0)={x_off,y_off}
				break	
				print hook_keyboard
			case hook_keyboard:
				wave counter
				wave trace_eval
				nvar peakn
				switch(s.keycode)
				
					case 31:	//key arrow down =  slow down, same peak
						nvar slower
						
						if (pc == 0)
							slower+=2
						elseif (pc == 1)
							slower+=2
						endif
						update_locWindow()
						break
					case 30:	//key arrow up = refresh; same peak, reset speed
						nvar slower
						slower=0
						update_locWindow()
						break
					case 29:	//key arrow right = next peak
						nvar view_acc 
						wave acc_traces
						
						if(view_acc ==0)
							peakn +=1
						else
							nvar view_ind, num_acc
							if(view_ind < num_acc)
								view_ind +=1
								peakn = acc_traces[view_ind]
							else
								view_ind = 0			
								peakn = acc_traces[view_ind]
							endif
						endif
						
						ModifyGraph/W=locWindow offset={0,0}	//remove any offset
						update_locWindow()
						break
					case 28:	//key arrow left = previous peak
						nvar view_acc
						wave acc_traces
						if(view_acc==0)
							peakn -= 1
						else
							nvar view_ind		
							if(view_ind > 0)
								view_ind-=1
								peakn = acc_traces[view_ind]
							else
								view_ind = 0
								peakn = acc_traces[view_ind]
							endif
						endif
						ModifyGraph/W=locWindow offset={0,0}	//remove any offset
						update_locWindow()
						break
						
					case 102: //key = f - fit circle with user specified bounds. only set up for SINGLE color.
						nvar maxcircfr
						if(exists("startfr"))
							nvar startfr
							nvar endfr
						else
							variable/g startfr
							variable/g endfr
						endif
						variable sf_temp = startfr
						variable ef_temp = endfr	
						do		
						Prompt sf_temp, "Enter startfr: "
						Prompt ef_temp, "Enter endfr: '"
						DoPrompt "Enter startfr and endfr", sf_temp, ef_temp
						while(!(ef_temp-sf_temp > 50) || !(ef_temp-sf_temp<maxcircfr) && V_flag ==0)
						
						if(V_flag==0) //if 'continue' pressed		
							make/o/n=(ef_temp-sf_temp) xxf, yyf
							wave xx,yy
							xxf= xx(sf_temp+p)
							yyf= yy(sf_temp+p)		
							//Display yyp vs xxp		
							fitcircle(xxf,yyf)
							wave circleguess, circpos0
							nvar peakn, recalc
							circpos0[peakn][0,2] = circleguess[q]
							recalc = 1 //new fit --> need to recalculate trace
							update_locWindow()
						endif
						break
						
					case 114:	//key r=reset graph
						nvar xc,yc,fov
						wave xxc,yyc
						
							//define center as mean of localizations
						wavestats/m=1/q xxc
						xc=V_avg
						wavestats/m=1/q yyc
						yc=V_avg
						SetAxis/W=locWindow left (yc+fov),(yc-fov)
						SetAxis/W=locWindow bottom (xc-fov),(xc+fov)
						ModifyGraph/W=locWindow offset={0,0}
						ModifyGraph/W=locWindow msize(yc0)=10,msize(yc0)=10
						break
					case 115: //key: s = startfr control
						if(exists("startfr"))
							nvar startfr
						else
							variable/g startfr
						endif
						
						if(exists("startfr_flag"))
							nvar startfr_flag
						else
							variable/g startfr_flag
						endif
						//prompt for startfr value	
						variable sf_t = startfr			
						Prompt sf_t, "Enter startfr: "
						DoPrompt "Enter startfr", sf_t
						if(V_flag==0) //if 'continue' pressed
						startfr = sf_t
						startfr_flag = 1
						
						nvar peakn
						update_locWindow()
						endif
					break
						
					case 97: ///key a= accept trace
						if(!exists("copied_ttc"))
							variable/g num_acc = 0
							make/n=(0)/o acc_traces=NaN
							make/t/n=(0)/o copied_ttc
							make/n=(0)/o acc_bound1 = NaN
							make/n=0/o acc_bound2 = NaN
							make/n=(0)/o acc_deltaAngle= NaN //trimmed delta angle
							make/n=(0)/o acc_avgrate = NaN // average rate over entire picked window. deg/point
							make/n=0/o acc_locstart = NaN //when does the trace start having localizations on ~every frame?
							make/n=0/o acc_effStart = NaN // effective start frame, from .tr file
							make/n=0/o acc_init=NaN 
							make/n=0/o acc_cp1 = NaN //optional: specify additional time points
							make/n=0/o acc_cp2 = NaN
							make/n=0/o acc_trdir=NaN
							make/n=0/o/t acc_trs = ""
							make/n=0/o acc_circx = NaN
							make/n=0/o acc_circy = NaN
							make/n=0/o acc_cstart = NaN // trace start frame from cstartfr
						
							
							//for sticking detection
							make/n=0/o stick_cat = 0 //=0 for no stick, 1 for stick
							make/n=0/o/t copied_ddc
							make/n=0/o stickquant // hold metric being used to decide on sticking
							
							make/n=0/o/t copied_c0int
							
							edit acc_locstart
							appendtotable acc_bound1
							appendtotable acc_bound2
							appendtotable acc_cp1
							appendtotable acc_cp2
						else
							wave acc_bound1,acc_bound2,acc_avgrate, acc_deltaAngle, acc_traces, acc_locstart,acc_init, acc_effStart
							wave/t copied_ttc
							nvar num_acc
						endif
						
						wave circpos0
						
						if (circpos0[peakn][0] <10^-6)
							print "cannot accept this trace until a center position is defined"
						else
							FindValue/v=(peakn)/T=0.5 acc_traces
							if(V_value == -1)
								print "accepting trace " + num2str(peakn)
								wave acc_bound1,acc_bound2,acc_avgrate, acc_deltaAngle, acc_traces, acc_locstart,acc_init, acc_effStart
								wave/t copied_ttc, copied_ddc, copied_c0int
								wave dirlist,acc_trdir, acc_circx, acc_circy, acc_cstart
								wave/t trlist, acc_trs, movie_number
								nvar num_acc
								redimension/n=(num_acc+1,-1) acc_bound1
								redimension/n=(num_acc+1,-1) acc_bound2
								redimension/n=(num_acc+1) acc_traces
								redimension/n=(num_acc+1,-1) acc_deltaAngle
								redimension/n=(num_acc+1) acc_avgrate 
								redimension/n=(num_acc+1) acc_deltaAngle
								redimension/n=(num_acc+1) copied_ttc
								redimension/n=(num_acc+1), acc_locstart
								redimension/n=(num_acc+1),acc_init
								redimension/n=(num_acc+1), acc_effStart
								redimension/n=(num_acc+1),acc_cp1
								redimension/n=(num_acc+1),acc_cp2
								redimension/n=(num_acc+1),acc_trdir,acc_trs
								redimension/n=(num_acc+1),stick_cat
								redimension/n=(num_acc+1), copied_ddc
								redimension/n=(num_acc+1), stickquant
								redimension/n=(num_acc+1) acc_circx, acc_circy, acc_cstart
								redimension/n=(num_acc+1) copied_c0int
								acc_traces[num_acc] = peakn
								string acc_name = "ttc_m"+movie_number[dirlist[peakn]] +"_"+ num2str(peakn)
								string ddc_name = "ddc_m"+movie_number[dirlist[peakn]] +"_"+ num2str(peakn)
								string int_name = "c0int_m"+movie_number[dirlist[peakn]]+"_"+num2str(peakn)
								acc_trdir[num_acc] = dirlist[peakn]
								acc_trs[num_acc]=trlist[peakn]
								duplicate ttc, $acc_name
								duplicate ddc, $ddc_name
								duplicate c0int, $int_name
								copied_ttc[num_acc] = acc_name
								copied_ddc[num_acc] = ddc_name
								copied_c0int[num_acc] = int_name
								
								acc_circx[num_acc] = circpos0[peakn][1]
								acc_circy[num_acc] = circpos0[peakn][2]
								nvar cstartfr
								acc_cstart[num_acc] = cstartfr
								
								num_acc +=1

							else 
								print "this trace has already been accepted"
							endif
						endif
					
						break
						
						
					case 99:	//key c=clear
						wave yyc,yyc1,yc0,yc1
						yyc=nan
						yyc1=nan
						yc0=nan
						yc1=nan
						break
					case 122:	//key z=zoom
						SetAxis/W=locWindow left 0,-2
						SetAxis/W=locWindow bottom -1,1
						ModifyGraph/W=locWindow msize(yc0)=20,msize(yc1)=20
						break
					case 117:	//key u=update angle graph etc, with current offset
						wave xxc1,yyc1
						variable/g update_center_flag=1	//signals that the center position should be updated
						x_off=str2num(stringbykey("XOFFSET",offs))
						y_off=str2num(stringbykey("YOFFSET",offs))
						xxc1+=x_off
						yyc1+=y_off
						update_distance_angle()
						xxc1-=x_off
						yyc1-=y_off
						break
							
						//trace classification
					case 49:	//1
						trace_eval[peakn]=1
						break
					case 50:	//2
						trace_eval[peakn]=2
						break
					case 51:	//3
						trace_eval[peakn]=3
						break
					case 52:	//4
						trace_eval[peakn]=4
						break
					case 53:	//5
						trace_eval[peakn]=5
						break
					case 54:	//6
						trace_eval[peakn]=6
						break
				endswitch
			//	histogram/b=2 trace_eval, counter		//update counter //to annotate traces
			endswitch
	endif
	return 0
End



Function AngleHook(s)	// Capture mouse/keyboard events in angle window
	STRUCT WMWinHookStruct &s
	
	NVAR TraceLevel = root:TraceLevel
	NVAR TraceStepSize =  root:TraceStepSize
	Variable MousePoint, RightClick, NewTraceLevel
	
		if(stringmatch(winname(0,3),"Angle_window"))
		switch(s.eventCode)
			case hook_mousedown:
				RightClick= (s.eventMod & (2^4))			 // bit 4, detect right click
				if(!RightClick)
					
//					wave current_trace = ttc
					
					MousePoint = AxisValFromPixel("Angle_window","left",s.mouseLoc.v)
//					MousePoint = str2num(stringbykey("HITPOINT",Tracefrompixel(s.mouseLoc.h,s.mouseLoc.v,""),":",";"))
//					print MousePoint
					NewTraceLevel = MousePoint
//					NewTraceLevel = mean(current_trace, pnt2x(current_trace, MousePoint-10), pnt2x(current_trace, MousePoint+10)) 		//average over 20 points
					TraceStepSize = NewTraceLevel - TraceLevel
					TraceLevel = NewTraceLevel					
				endif
				break	
		
			case hook_mouseup:
				break	
				
			case hook_keyboard:
				switch(s.keycode)
					case 49:								//key 1=store delta_angle1 value
						wave delta_angle1
						nvar peakn
						delta_angle1[peakn]=TraceStepSize
						print delta_angle1[peakn]
						break
					case 50:								//key 2=store delta_angle1 value
						wave delta_angle2
						nvar peakn
						delta_angle2[peakn]=TraceStepSize
						print delta_angle2[peakn]
						break
					case 51:								//key 3=store delta_angle1 value
						wave delta_angle3
						nvar peakn
						delta_angle3[peakn]=TraceStepSize
						print delta_angle3[peakn]
						break
					case 45:								//key "-'=minus 360 deg
						wave delta_angle
						nvar peakn
						TraceStepSize-=360
						print TraceStepSize				
						break
					case 61:								//key '='=plus 360 deg
						wave delta_angle
						nvar peakn
						TraceStepSize+=360
						print TraceStepSize				
						break
				endswitch
			endswitch
		
		endif		
	return 0
End



////*********************************Graphs**************************************//
// 
// 
Window LocalizationsGraph() : Graph
	PauseUpdate; Silent 1		// building window...
		//build window
	Display/K=1/W=(1,45,746,734) yyd vs xxd		//ch0 trace
	Appendtograph yyd1 vs xxd1							//ch1 trace
	Appendtograph yc0 vs xc0						//current ch0 loc
	Appendtograph yc1 vs xc1							//current ch1 loc
	DoWindow/C locWindow
	SetWindow locWindow,hook(testhook)=KeyboardPanelHook
	ModifyGraph gfSize=20,height={Aspect,1},gbRGB=(5000,5000,5000)
	ModifyGraph lSize(yyd)=2,lSize(yyd1)=2
	if(stringmatch(ch0col,"g"))
	ModifyGraph zColor(yyd)={zz,*,*,green,0}
	ModifyGraph mode(yc0)=3,marker(yc0)=19,msize(yc0)=20,rgb(yc0)=(0,65535,0)
	else
	ModifyGraph zColor(yyd)={zz,*,*,red,0}
	ModifyGraph mode(yc0)=3,marker(yc0)=19,msize(yc0)=20,rgb(yc0)=(65535,0,0)
	endif
	if(stringmatch(ch1col,"g"))
	ModifyGraph zColor(yyd1) = {zz,*,*,green,0}
	ModifyGraph mode(yc1)=3,marker(yc1)=19,msize(yc1)=10,rgb(yc1)=(0,65535,0)
	else
	ModifyGraph zColor(yyd1)={zz,*,*,red,0}	
	ModifyGraph mode(yc1)=3,marker(yc1)=19,msize(yc1)=10,rgb(yc1)=(65535,0,0)
	endif
	ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={0,0}
	ModifyGraph manTick(bottom)={0,1,0,0},manMinor(bottom)={0,0}
	ModifyGraph mirror=1
	ModifyGraph grid=1,gridRGB=(8738,8738,8738)
	//ModifyGraph lblMargin(left)=36,lblMargin(bottom)=23
	ModifyGraph standoff=0
	Label left "Y"
	Label bottom "X"
	SetAxis left (yc+fov),(yc-fov)
	SetAxis bottom (xc-fov),(xc+fov)
	TextBox/N=text0/A=RT/F=0/B=1/G=(65535,65535,65535) "0"
		
//		//add beam
//	appendtograph beam_y vs beam_x
//	ModifyGraph mode(beam_y)=0,lsize(beam_y)=10,rgb(beam_y)=(65535,65535,65535)
		
		//enable quickdrag for red trace
	newdatafolder/o root:WinGlobals
	newdatafolder/o root:WinGlobals:locWindow
	string/G root:WinGlobals:locWindow:S_TraceOffsetInfo
	ModifyGraph quickdrag(yyd1)=1
	

	appendtograph/w=locWindow c0y vs c0x
	ModifyGraph/w=locWindow mode(c0y) = 0, marker(c0y) = 19, msize(c0y) = 10,rgb(c0y) =(24576,24576,65535)
		appendtograph/w=locWindow c1y vs c1x
		ModifyGraph/w=locWindow mode(c1y) = 0, marker(c1y) = 19, msize(c1y) = 10,rgb(c1y) =(25500,12800,0)
EndMacro

 Window DistanceGraph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(758,45,1439,381) ddc
	DoWindow/C Distance_window
	if (pc == 0)
	ModifyGraph gfSize=16,width=550,height=250,gbRGB=(4369,4369,4369)
	endif
	if (pc == 1)
	ModifyGraph gbRGB=(4369,4369,4369)
	endif
	ModifyGraph mode=0
	ModifyGraph marker=19
	ModifyGraph lsize(ddc)=1
	ModifyGraph rgb(ddc)=(32768,54615,65535)
	ModifyGraph msize(ddc)=2
	Label left "separation distance (px)"
	SetAxis left 0,1.1
EndMacro


Window AngleGraph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(750,407,1441,743) ttc
	DoWindow/C Angle_window
	if (pc == 0)
		ModifyGraph gfSize=16,width=550,height=250,gbRGB=(4369,4369,4369)
	endif
	if (pc == 1)
		ModifyGraph gbRGB=(4369,4369,4369)
	endif
	ModifyGraph mode(ttc)=0
	ModifyGraph lsize(ttc)=1
	ModifyGraph rgb(ttc)=(40000,40000,40000)
	Label left "angle"
	//SetAxis left -180,180
	SetWindow Angle_window,hook(testhook2)=AngleHook	//capture mouse&keybord events in graph
EndMacro

Window StepSizeTrace() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(250,407,941,543) asst
	DoWindow/C Stepsize_window
	if (pc == 0)
		ModifyGraph gfSize=16,width=550,height=250,gbRGB=(4369,4369,4369)
	endif
	if (pc == 1)
		ModifyGraph gbRGB=(4369,4369,4369)
	endif
	ModifyGraph mode(asst)=2
	ModifyGraph lsize(asst)=4
	ModifyGraph rgb(asst)=(65535,32768,0)
	Label left "angular step size"
	SetAxis left -200,200
EndMacro

Window RadStepSizeTrace() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(250,407,941,543) rsst
	DoWindow/C RadStepsize_window
	if (pc == 0)
		ModifyGraph gfSize=16,width=550,height=250,gbRGB=(4369,4369,4369)
	endif
	if (pc == 1)
		ModifyGraph gbRGB=(4369,4369,4369)
	endif
	ModifyGraph mode(rsst)=2
	ModifyGraph lsize(rsst)=4
	ModifyGraph rgb(rsst)=(0,32768,65535)
	Label left "radial step size"
	SetAxis left -0.5,0.5
EndMacro

Window IntTrace() : Graph
	PauseUpdate; Silent 1		// building window...
	Display/K=1 /W=(300,450,941,543) c0int,c1int
	DoWindow/C Intensity_window
	if (pc == 0)
		ModifyGraph gfSize=16,width=550,height=250,gbRGB=(4369,4369,4369)
	endif
	if (pc == 1)
		ModifyGraph gbRGB=(4369,4369,4369)
	endif
	ModifyGraph mode(c0int)=0,mode(c1int) = 0
	ModifyGraph lsize(c0int)=1, lsize(c1int) = 1
	//ModifyGraph rgb(vertline2)=(34952,34952,34952),rgb(vertline)=(34952,34952,34952)
	if(stringmatch(ch0col, "g"))
		ModifyGraph rgb(c0int)=(0,65535,0)
		ModifyGraph rgb(c1int) = (65535,0,0)
	else
		ModifyGraph rgb(c1int)=(0,65535,0)
		ModifyGraph rgb(c0int) = (65535,0,0)
	endif
	
	
	
	Label left "Intensity"
EndMacro


////////////////////////////////////////
/////////////////////////////////////// "SHORTCUTS"
//
////close all windows
function CloseWindows()
	DoWindow/K RadStepsize_window
	DoWindow/K Stepsize_window
	DoWindow/K Angle_window
	DoWindow/K Distance_window
	DoWindow/K locWindow
	DoWindow/K Intensity_window
end
//
//
function doit()
	begin()
	browse_locs()
	get_2D_hist_all()
end
////////////////////////////////////
////////////////////////////////////
//

////////////////////////////////////////////////////
////////////////////////////////////////////////////GENERATE PROCCESSED DATA WAVES - ANGLE, DISTANCE, INTENSITY,....
//
	//calculate separation distance and angle for current trace
function update_distance_angle()
	wave xxc,yyc,xxc1,yyc1,ddc,ttc,asst,rsst, ttc_temp, ttc_tempfilt
	nvar fc, trlen
	nvar ill_set
	nvar peakn
	variable/g provrad = 0.3
	wave circpos0
	nvar update_center_flag
	if(fc==0 && ill_set ==0 && update_center_flag)//note this only happens if update_center_flag manually set to 1 by pressing 'u'
		nvar xcenter,ycenter
		wave circpos0
		//xxc1=xcenter
		//yyc1=ycenter
		circpos0[peakn][0]=provrad	//provisional radius
		circpos0[peakn][1]=xcenter
		circpos0[peakn][2]=ycenter
		variable/g update_center_flag=0
	endif
	
	//create median filtered positions
	nvar trackfilt
	duplicate/o xxc,xxcfilt
	duplicate/o yyc, yycfilt
	duplicate/o xxc1,xxc1filt
	duplicate/o yyc1, yyc1filt
	smooth/m=0 trackfilt, xxcfilt, yycfilt, xxc1filt, yyc1filt
	
	
	if(ill_set != 0)  //two color
		ddc  = sqrt((xxc[p]-xxc1[p])^2+(yyc[p]-yyc1[p])^2)
		ttc_temp	=(180/pi)*(imag(r2polar(cmplx(xxc[p]-xxc1[p],-(yyc[p]-yyc1[p])))))
		ttc_tempfilt=(180/pi)*(imag(r2polar(cmplx(xxcfilt[p]-xxc1filt[p],-(yycfilt[p]-yyc1filt[p])))))

		duplicate/o ttc_temp,ttc_raw
		clean_theta()
	
	else //single color:
		wave circpos0
		ddc   = sqrt((xxc[p]-circpos0[peakn][1])^2+(yyc[p]-circpos0[peakn][2])^2)
		ttc_temp	=(180/pi)*(imag(r2polar(cmplx(xxc[p]-circpos0[peakn][1],-(yyc[p]-circpos0[peakn][2])))))	
		ttc_tempfilt	=(180/pi)*(imag(r2polar(cmplx(xxcfilt[p]-circpos0[peakn][1],-(yycfilt[p]-circpos0[peakn][2])))))	

		duplicate/o ttc_temp,ttc_raw
		clean_theta()

	endif
	
		//rescale angle graph
	wavestats/q/m=1 ttc_temp
	SetAxis/w=Angle_window left V_min-10,V_max+10
	duplicate/o ttc_temp,ttc
	
	//rescale frame axis for all
	variable xstart = 0
	variable xstop = trlen
	SetAxis/W=Distance_window bottom xstart,xstop
	SetAxis/W=Angle_window bottom xstart,xstop
	SetAxis/W=Stepsize_window bottom xstart,xstop
	SetAxis/W=RadStepsize_window bottom xstart,xstop
	SetAxis/W=Intensity_window bottom xstart,xstop
	
		//angle&radial step size traces
	asst = ttc[p] - ttc[p-1]
	asst[0] = NaN
	rsst = ddc[p] - ddc[p-1]
	rsst[0] = NaN
	
		//update filtered angle trace
	Duplicate/O ttc,ttc_smth
	Smooth/M=0 8, ttc_smth
END

	//cleans theta trace to remove >180deg jumps
function clean_theta()

	wave ttc_temp, ttc_tempfilt
	nvar trlen, trackfilt
	variable prev_theta
	variable i
	
		//skip to first non-nan value
	i=-1
	do
		i+=1
	while(numtype(ttc_tempfilt[i])==2 && i<trlen)
	variable st = i
	prev_theta=ttc_tempfilt[i]
	variable adj = 0
	for(;i<trlen;i+=1)
		if(numtype(ttc_tempfilt[i])!=2)		//if not nan
			if(abs(ttc_tempfilt[i]+adj-prev_theta)>180)
				if((ttc_tempfilt[i]+adj)>prev_theta)
					adj -=360
					//ttc_temp[i,*]-=360
				else
					adj +=360
					//ttc_temp[i,*]+=360
				endif
			endif
			ttc_tempfilt[i] += adj
			prev_theta=ttc_tempfilt[i]
		endif
	endfor
	
	if(trackfilt ==1) //no median filter was actually applied. be efficient.
		ttc_temp = ttc_tempfilt 
	else //constrain ttc_temp to be within 180 deg of ttc_tempfilt
		i=st
		variable diff
		for(;i<trlen;i+=1)	
			do
				if((ttc_temp[i] -ttc_tempfilt[i]) > 0)
					ttc_temp[i] -= 360
				else
					ttc_temp[i] += 360
				endif
			while(abs(ttc_tempfilt[i] - ttc_temp[i]) >180)
		endfor
	endif
END


////////////////////////
///////////////////////


///////////////////////
//////////////////////CIRCLE FITTING FUNCTIONS
	//find center and size of circle (single color at a time)
	//11jun15 bda
	//see http://www.wavemetrics.net/doc/igorman/III-08%20Curve%20Fitting.pdf
function fitcircle(xxf,yyf) //takes x and y pos waves as inputs.
	wave xxf, yyf
	variable xcenter, ycenter, radius
	variable xguess, yguess, rguess
	
	//initial guesses
	rguess = 0.5
	Wavestats/q xxf
	xguess = V_avg
	wavestats/q yyf
	yguess = V_avg
	
	make/o/d circleguess = {rguess, xguess, yguess}
	//print circleguess
	Duplicate/o xxf, xxfit, yyfit
	variable V_FitError = 0
	//add constraints
	Make/o/T ConstText = {"K0 > 0","K0 < 5","K1>"+num2str(xguess-3),"K1<"+num2str(xguess+3),"K2>"+num2str(yguess-3),"K2<"+num2str(yguess+3)}
	FuncFit/ODR=3/q/n/w=2 circlefn, circleguess /X ={xxf, yyf} /XD={xxfit,yyfit} /C=ConstText //note result ends up in circleguess
END

	//fit function for fitcircle
Function circlefn(w,x,y) : FitFunc
	wave w
	variable x, y
	//CurveFitDialog/
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = r
	//CurveFitDialog/ w[1] = x0
	//CurveFitDialog/ w[2] = y0
	return ((x-w[1])/w[0])^2 + ((y-w[2])/w[0])^2 - 1  //=0, for implicit fitting.
END

//puts centers of circles on the two channels into circpos0 and circpos1
//based on preset frames, not user specified for that trace.
Function FindCenters(mode)
//mode = 0 for single color, mode = 1 for two color
variable mode
nvar peakn, trlen
wave xxc,yyc,xxc1,yyc1,circpos0,circpos1,circleguess
//wave ch1xcent,ch1ycent,ch2xcent,ch2ycent,
//wave c0x,c0y,c1x,c1y,theta, 

//for mode = 0 or 1, find the first circle

variable startfr = frame_padding
	
	variable endfr = min(frame_padding + 500, trlen)
	//don't let length for fit be >1000 to make faster.
	//fit this subset of frames to a circle.
	make/o/n=(endfr-startfr) xxp, yyp
	xxp= xxc[startfr+p]
	yyp = yyc[startfr+p]
	duplicate/o xxp,xxf
	duplicate/o yyp,yyf
	 fitcircle(xxp,yyp)
	 wave circleguess //need to 'recall' after fitcircle.
	circpos0[peakn][] = circleguess[q]

if(mode==1)
	//if mode =1 also find second channel circle fit
	make/o/n=(endfr-startfr) xxp1, yyp1
	xxp1= xxc1[startfr+p]
	yyp1 = yyc1[startfr+p]
	 fitcircle(xxp1,yyp1)
	 wave circleguess
	 circpos1[peakn]= circleguess[q]
endif	
END

///Can add circle_batch from v5 if desired. But, no longer in use since 2d histogram wall added

/////////////////////////////
/////////////////////////////




///////////////////////////
//////////////////////////FUNCTIONS FOR QUICKLY GOING THROUGH TRACES
//	//generate 2D histograms of all traces and displays them as a scrollable wall
function get_2d_hist_all()
	nvar load_2dhist
	
	variable normalize=1	// 0 - no normalization; 1 - normalizes each trace histogram by number of localizations 
	variable/g nrows=20		//histograms per column in wall view
	variable/g zsize
	variable hist_frnum
	nvar pc
	
	if(load_2dhist ==0) // generate 2d histograms here. 
		variable/g resolution=30
		hist_frnum = 1000 //approximate number of frames to use for each 2d histogram
		get_hist_sets(hist_frnum)
		wave hist_sets
		zsize = dimsize(hist_sets,0)
		
		wave xx,yy, xx1,yy1
		wave xxc,yyc
		nvar fov, peakn
		variable xc,yc
	
			//create and scale 2d hist wave
		make/n=(resolution,resolution,zsize)/o t2d_hist_all=0
		setscale/i x (-fov),fov,t2d_hist_all
		setscale/i y (-fov),fov,t2d_hist_all
	
		//reference waves for mapping
		make/n=(dimsize(t2d_hist_all,0))/o du_x
		setscale/p x dimoffset(t2d_hist_all,0),dimdelta(t2d_hist_all,0),du_x
		make/n=(dimsize(t2d_hist_all,1))/o du_y
		setscale/p x dimoffset(t2d_hist_all,1),dimdelta(t2d_hist_all,1),du_y
	
	//waves for partial traces
		make/o/n=(hist_frnum) xxcp, yycp
		
		variable trace,i, xn,yn
		nvar histch1
		
	
		for(trace=0;trace<zsize;trace+=1) //'trace' is now a partial trace, duration set by hist_frnum
				//load current trace
			if(!mod(trace,1000))
				print trace
			endif
			nvar peakn
			if(hist_sets[trace][0] != peakn)
				loadtrace(hist_sets[trace][0])
			endif
					
			if(histch1 ==0) //make histogram with ch0
			xxcp=xx[hist_sets[trace][1]+p][hist_sets[trace]*2] //partial trace
			yycp=yy[hist_sets[trace][1]+p][hist_sets[trace]*2]
		elseif(histch1==1)
			xxcp=xx1[hist_sets[trace][1]+p][hist_sets[trace]*2] //partial trace
			yycp=yy1[hist_sets[trace][1]+p][hist_sets[trace]*2]
			endif
			//define center as mean of localizations; offset traces to be centered at origin
			wavestats/m=1/q xxcp
			xxcp-=V_avg
			wavestats/m=1/q yycp
			yycp-=V_avg
			nvar trlen
			for(i=0;i<trlen;i+=1)
				if(numtype(xxcp[i])!=2)
					xn=x2pnt(du_x,xxcp[i])
					yn=x2pnt(du_y,yycp[i])
				t2d_hist_all[limit(xn,0,resolution-1)][limit(yn,0,resolution-1)][trace]+=1
				endif
			endfor
			if(normalize)
					//normalize by number of localizations
				t2d_hist_all[][][trace]/=V_npnts
			endif
		endfor
	else //load 2d histograms.
		wave/t find_data
		
		//load  info filef
		variable dataref
		variable/g resolution
		variable hfov
		Open/Z/R dataref as find_data[0] + "histpar.info"
		FBinRead/F=3/U dataref, resolution
		FBinRead/F=3/U dataref, hfov
		FBinRead/F=3/U dataref, hist_frnum
		FBinRead/F=3/U dataref, zsize
		Close dataref
	
		print "resolution in loaded 2d histogram: " +num2str(resolution)
		print "fov in loaded 2d histogram: " + num2str(hfov)
		
		//load hist_sets
		make/o/n=(3,zsize) hist_sets
		Open/Z/R dataref as find_data[0] + "histsets.list"
		FBinRead/F=3/U dataref, hist_sets
		matrixtranspose hist_sets
		Close dataref
		
		//load stack of histograms, one at a time, normalizing each if appropriate
		make/n=(resolution,resolution,zsize)/o t2d_hist_all=0
		make/n=(resolution,resolution)/o t2d_hist = 0
		Open/Z/R dataref as find_data[0] + "hist2d.stack"
		i=0
		for(;i<zsize;i+=1)
			FBinRead/F=3/U dataref,  t2d_hist
			matrixtranspose t2d_hist
			if(normalize)
				wavestats/q t2d_hist
				t2d_hist /= V_sum
			endif
			t2d_hist_all[][][i] = t2d_hist[p][q]
		endfor
	endif
	
	
	variable/g ncols=ceil(zsize/nrows)
	make/n=(nrows*resolution,ncols*resolution)/o t2d_hist_wall=0
	
	variable/g axis_left = -0.5
	variable/g axis_right = 1.5*dimsize(t2d_hist_wall,0)-0.5
	if(pc==1)
		axis_right = axis_right * 0.5
	endif
	variable/g axis_step = 0.9*axis_right
	
		//note: rows and columns are flipped in images 
	t2d_hist_wall=t2d_hist_all[mod(p,resolution)][mod(q,resolution)][min(zsize-1,floor(q/resolution)*nrows+floor(p/resolution))]
	
		//build graph
	execute "t2d_histGraph()"
	variable satcount = 15 /hist_frnum
	if(normalize)
		ModifyImage/W=t2D_hist_wall_window t2d_hist_wall ctab= {*,satcount,YellowHot,0}
	else
		ModifyImage/W=t2D_hist_wall_window t2d_hist_wall ctab= {*,30,YellowHot,0}
	endif
END


//////split each trace into pieces for 2d histogram.
function get_hist_sets(frnum)
	variable frnum
	nvar all_ntr
	wave expt_frEff
	wavestats/q expt_frEff
	make/o/n=(max(all_ntr * V_max/ frnum,all_ntr),3) hist_sets
	hist_sets = NaN
	variable curr = 0
	variable trace, dur, numpieces, endfr, startfr
	//split up each trace appropriately.
	
	for(trace=0;trace<all_ntr;trace+=1)
		//load current trace
		loadtrace(trace)
		wave xx
		nvar trlen
		//nvar estartfr //effective start frame
		numpieces = ceil(trlen/frnum)
		hist_sets[curr,curr+numpieces-1;1][0] = trace
		hist_sets[curr,curr+numpieces-1;1][1] =  frnum*(p-curr) //index. not true frame.
		hist_sets[curr,curr+numpieces-1;1][2] = frnum*(p-curr+1)
		curr +=numpieces
	endfor
	//get rid of empty spaces at the end
	redimension/n=(curr,3) hist_sets
END

//*********************************Graphs**************************************//

Window t2d_histGraph() : Graph
	PauseUpdate; Silent 1		// building window...
	if(pc==0)
	Display/K=1/W=(0,45,1433,893)
	else
	Display/K=1/W=(716,45,1433,800)
	endif
	DoWindow/C t2D_hist_wall_window
	AppendImage/T t2d_hist_wall
	ModifyImage t2d_hist_wall ctab= {*,*,YellowHot,0}
	ModifyGraph margin(left)=30,margin(bottom)=14,margin(top)=14,margin(right)=14,gfSize=16
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=10,nticks(top)=2
	ModifyGraph minor=1
	ModifyGraph fSize=9
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	//SetAxis/R left -0.5,1.5*dimsize(t2d_hist_wall,0)-0.5//299.5
	SetAxis/R left axis_left,axis_right//299.5
	//print axis_left
	//print axis_right
	SetAxis top -0.5,dimsize(t2d_hist_wall,0)-0.5//199.5
	ModifyGraph swapXY=1
	TextBox/C/N=textTraceNo/F=0/B=1/G=(0,0,0)/A=LT/X=-1.00/Y=-1.99 ""
	SetWindow t2d_hist_wall_window,hook(testhook3)=t2d_hist_Hook	//capture mouse&keybord events in graph
	
end

//*********************************Capturing input events**************************************//



Function t2d_hist_Hook(s)	// Capture mouse/keyboard events in histogram wall window
	STRUCT WMWinHookStruct &s

	Variable MousePoint_x,MousePoint_y, RightClick
	
	wave t2d_hist_all,t2d_hist_wall
	wave hist_sets
	variable resolution=dimsize(t2d_hist_all,0)
	variable nrows=dimsize(t2d_hist_wall,0)/resolution
	variable/g selectedTrace
	nvar axis_left,axis_right,axis_step
	
	if(stringmatch(winname(0,3),"t2D_hist_wall_window"))
	switch(s.eventCode)
		case hook_mousemoved:
				//set cursor icon so it is visible on dark background
			s.cursorCode=20
			s.doSetCursor=1
	
				//display trace number		
			MousePoint_y = AxisValFromPixel("t2D_hist_wall_window","left",s.mouseLoc.v)
			MousePoint_x = AxisValFromPixel("t2D_hist_wall_window","bottom",s.mouseLoc.h)
			//selectedTrace=2*(floor(MousePoint_x/resolution)*nrows+floor(MousePoint_y/resolution))
			selectedTrace=(floor(MousePoint_x/resolution)*nrows+floor(MousePoint_y/resolution))
			TextBox/W=t2d_hist_wall_window/C/N=textTraceNo "\Z10" + "peakn: " +  num2str(hist_sets[selectedTrace][0])	 + " ::: " + num2str(hist_sets[selectedTrace][1]) 	//trace number textbox
			break
			
		case hook_mousedown:
			RightClick= (s.eventMod & (2^4))			 // bit 4, detect right click
			if(!RightClick)
				MousePoint_y = AxisValFromPixel("t2D_hist_wall_window","left",s.mouseLoc.v)
				MousePoint_x = AxisValFromPixel("t2D_hist_wall_window","bottom",s.mouseLoc.h)
				selectedTrace= (floor(MousePoint_x/resolution)*nrows+floor(MousePoint_y/resolution))
			endif
			break	
			
		case hook_keyboard:
				nvar peakn,selectedTrace, pc
				switch(s.keycode)

					case 103:	//key g=go to trace
						peakn=hist_sets[selectedTrace][0]
						print "Jumped to trace " + num2str(peakn)
						update_locWindow()
						break
						
					case 28:  //left arrow key - shift left
						axis_left -= axis_step
						axis_right -= axis_step
						
						SetAxis bottom, axis_left, axis_right
					break
					case 29: //right arrow key - shift right
						axis_left += axis_step
						axis_right += axis_step
						
						SetAxis bottom, axis_left, axis_right

					break
				endswitch
		endswitch
	
	endif		
	return 0
END

////////////////////////////////////////////////////////
///////////////////////////////////////////////////////

/////////////////////////////////////////////////
///////////////////////////////////////////////// MISC. ANALYSIS FUNCTIONS
//
//display that isn't reset automatically all the time.
function dttc()
	display ttc
	ModifyGraph grid(left)=1,manTick(left)={0,360,0,1},manMinor(left)={0,0};DelayUpdate
	ModifyGraph gridRGB(left)=(0,0,0)
end

	//locates traces classified as 'val' in the trace_eval wave, and creates a list of these trace numbers
function get_traces(val)
	 variable val
	
	wave trace_eval
	make/n=0/o trace_list
	
	variable i=0
	
	do
		FindValue/v=(val)/s=(i) trace_eval
		redimension/n=(numpnts(trace_list)+1) trace_list
		trace_list[numpnts(trace_list)-1]=V_value
		i=V_value+1
	while(V_value!=-1 && i<numpnts(trace_eval))
	
		//remove last entry if it is -1
	if(trace_list[numpnts(trace_list)-1]==-1)
		redimension/n=(numpnts(trace_list)-1) trace_list
	endif
	
 end


/////////////////////////
/////////////////////////


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////ANALYSIS FUNCTIONS - TRIM AND DISPLAY ACCEPTED
////display all accepted traces, without trimming the ends
function dispFullAcc()
	nvar num_acc
	wave acc_traces
	wave/t copied_ttc
	string windowname = "Full accepted traces"
	DoWindow/K FullAcc

	Display/K=1/n=windowname $copied_ttc[0]
	DoWindow/C FullAcc

	variable i
	for(i=1;i<num_acc;i+=1)
	AppendToGraph $copied_ttc[i]
	endfor

	Label left "angle (degree)"
	Label bottom "frames"

END


//set bounds for part of a trace to exclude the beginning and end for display
//change_flag should usually be zero
//also, save the change in angle from initial to final time
//uses PCSR()
//and, determines the average rate, in degrees/point. if laser is not alternating, this is degrees/camera frame. 
//otherwise, be careful
//this works, but is inconvenient since you can only see the ttc window while picking bounds. better to manually enter numbers
//then, call this with change_flag = 0 to quickly calculate other stats
function AccTrim(change_flag)
	variable change_flag //=0: don't change existing bounds. only look at those with zero values.
	//=1 go back through all the traces and allow changing points. but use the previous values as starting points
	nvar num_acc, expt_frEff
	wave acc_traces, acc_bound1, acc_bound2, acc_deltaAngle, acc_avgrate, acc_locstart
	wave acc_init,acc_effStart, effStart
	wave/t copied_ttc
	
	make/o/n=(num_acc) acc_radstd = NaN // uses metric for radial std which excludes points > 2*mean in AccStick(i)
		
	DoWindow/K setbounds
	duplicate/o ttc,ttc_temp
	Display/K=1 ttc_temp //just to give it something to build the window with
	//TextBox/N=text0/A=RB " "
	Dowindow/C setbounds
	ShowInfo
	//variable dummy = 0
	//variable startx, endx
	
	acc_effStart = effStart[acc_traces[p]]
	variable i
	for(i=0;i<num_acc;i+=1)
		if(acc_bound2[i] == 0 || change_flag == 1)
		duplicate/o $copied_ttc[i],ttc_temp
//have to manually find bounds b/c don't want to reload the .tr
		//FindValue/v=0/T=(10^6) ttc_temp //find first not NaN entry. INDEX, not x
		//startx = max(V_value-100,0)
		duplicate/o ttc_temp, ttc_tempr
		//Reverse ttc_tempr
		//FindValue/v=0/T=(10^6) ttc_tempr //find last not NaN entry. INDEX.
	 	//endx = min(expt_frEff-1 - V_value + 100,expt_frEff)
	
		//SetAxis bottom, startx, endx
		TextBox/C/N=text0/A=RB "Trace: " + copied_ttc[i]
		cursor/a=1 a ttc_temp acc_bound1[i]
		cursor/a=1 b ttc_temp acc_bound2[i]
		variable pick_flag = 0
		do
			NewPanel/K=2 as "Pause for Cursor"
			Dowindow/c tmp_pauseforcursor
			AutopositionWindow/E/m=1/R=setbounds
			DrawText 21,20, "Adjust the cursors, then"
			Drawtext 21,40,"Click to continue"
			Button button0,pos={80,58},size={92,20},title = "continue"
			button button0, proc=UserCursorAdjust_ContButtonProc
			PauseForUser tmp_Pauseforcursor, setbounds
			if(xcsr(b) !=0)
			pick_flag = 1
			else 
				 print "Please try again"
			endif
		while(pick_flag == 0)
			acc_bound1[i] = pcsr(a)
			acc_bound2[i] = pcsr(b)
		endif
	endfor
		variable startx, endx
		//update stats for all accepted traces.
		i = 0
		for(;i<num_acc;i+=1)
			svar movie_number
			string acc_name = copied_ttc[i]
				duplicate/o $copied_ttc[i], ttc_tempr
				startx = acc_bound1[i]
				endx = acc_bound2[i]
				variable initx = acc_locstart[i]
				acc_deltaAngle[i] = ttc_tempr(endx) - ttc_tempr(startx)
				acc_avgrate[i] = acc_deltaAngle[i] / (endx- startx)
				acc_init[i] = startx - initx		
		
			//determine whether each trace has surface sticking
			AccStick(i)

		
		endfor
		DoWindow/K setbounds
		
		//determine if fits are acceptable
		fitqual_filt()
		//set fitfilt_cat based on accstick and fitqual_filt results
		filter_accepted()
		
		//generate stats for traces which were accepted by filtering
		make/o/n=(num_acc) fa_avgrate, fa_deltaAngle, fa_init, fa_traces //filtered accepted data traces
		make/o/n=(num_acc) fa_bound1, fa_bound2, fa_locstart
		variable/g num_filt = 0
		wave filt_val
		for(i=0; i<num_acc;i+=1)
			if(filt_val[i] ==0) //good trace
				fa_traces[num_filt] = acc_traces[i]
				fa_bound1[num_filt] = acc_bound1[i]
				fa_bound2[num_filt] = acc_bound2[i]
				fa_locstart[num_filt] = acc_locstart[i]
				
				fa_avgrate[num_filt] = acc_avgrate[i]
				fa_deltaAngle[num_filt] = acc_deltaAngle[i]
				fa_init[num_filt] = acc_init[i]
				
				num_filt +=1
			endif
		endfor
		redimension/n=(num_filt) fa_avgrate, fa_deltaAngle, fa_init, fa_traces, fa_bound1, fa_bound2, fa_locstart
		
		
		//pull out the intensity at the first and last frame of translocation.
		make/o/n=(num_acc) trstart_int, trstop_int, travg_int
		wave/t copied_c0int
		for(i=0;i<num_acc;i+=1)
			duplicate/o $copied_c0int[i], int_temp
			trstart_int[i] = int_temp[acc_bound1[i]]
			trstop_int[i]  = int_temp[acc_bound2[i]]
			wavestats/q/r=(acc_bound1[i],acc_bound2[i]) int_temp
			travg_int[i] = V_avg
		endfor
		
end

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor // Kill panel
End



//for indicated trace, determine whether it shows signs of surface sticking
//set up for single color only
//i is accepted trace number, not peakn.
Function AccStick(i)
	variable i	
	nvar ill_set
	if(ill_set !=0)
		print "Warning - not set up for dual color"
	endif
	wave/t copied_ttc, copied_ddc
	wave stick_cat, circpos0, acc_traces,acc_bound1,acc_bound2, stickquant, acc_radstd
	nvar num_acc
	nvar std_window, std_thr, std_freq
	
	variable translocfr, x1,x2, ddc_std_tr, trl, j, counter
		translocfr = acc_bound2[i] - acc_bound1[i]
		
		
		duplicate/o $copied_ttc[i], ttc_tempa
		duplicate/o $copied_ddc[i], ddc_temp
		wavestats/q/m=1 ttc_tempa
		trl = V_npnts
		
		duplicate/o ttc_tempa,pttc
		pttc = ttc_tempa * ((circpos0[acc_traces[i]][0] * 2 * pi) / 360) // arc motion, converted to pixels
		duplicate/o ttc_tempa,pttc_std,ttc_std,ddc_std
		pttc_std = NaN
		ttc_std = NaN
		ddc_std = NaN
		
		x1 = pnt2x(ddc_temp, acc_bound1[i])
		x2 = pnt2x(ddc_temp,acc_bound2[i])
		//improvement 20160530 BDA - remove outliers from ddc before getting ddc_std_tr --
		//generally are due to other stuff floating by. prevents false positives on sticky calling due to couple frames 
		//with fits way off the circle
		duplicate/o ddc_temp, ddc_temp2
		variable fr = x1
		variable count = 0
		wavestats/q/r=(x1,x2) ddc_temp
		for(;fr<x2+1;fr+=1)
			if(ddc_temp2[fr]>(V_avg*2))
				ddc_temp2[fr] = NaN
				count +=1
			endif
		endfor
		//print i
		//print count
		if(count/(x2-x1+1)> 0.05)
			print "For accepted trace " + num2str(i) + " >5% of ddc values were removed during AccStick"
		endif
		ddc_std_tr = variance(ddc_temp2,x1,x2)^0.5
		acc_radstd[i] = ddc_std_tr
		
		//calculate std for each sliding window; trim down after.
		//position i in the std trace is for the window starting at i, not centered at i.
		j = 0
		for(;j<trl-std_window+1;j+=1)
			x1 = pnt2x(ddc_temp, j)
			x2 = pnt2x(ddc_temp,j+std_window-1)
			pttc_std[j] = variance(pttc,x1,x2)^0.5
			ttc_std[j] = variance(ttc_tempa,x1,x2)^0.5
			ddc_std[j] = variance(ddc_temp,x1,x2)^0.5
		endfor
		
		duplicate/o pttc_std, std_ratio
		std_ratio = pttc_std / ddc_std_tr
		//print ddc_std_tr
		
		
		//how often is std_ratio below the threshold?
		j = acc_bound1[i]
		counter = 0 
		for(;j<acc_bound2[i] -std_window+1;j+=1)
			if(std_ratio[j] < std_thr)
				counter +=1
			endif
		endfor
		counter = counter / (acc_bound2[i] - acc_bound1[i]-std_window+1)
		stickquant[i] = counter
		if(counter < std_freq) //not sticky
			stick_cat[i] = 0
		else
			stick_cat[i] = 1			
		endif
		doupdate
END

function fitqual_filt()
	nvar num_acc, fitqual_thr
	wave acc_radstd, stick_cat, stickquant
	make/n=(num_acc)/o fitfilt_cat = NaN

	variable i = 0
	for(;i<num_acc;i+=1)
		if(acc_radstd[i] < fitqual_thr)//good trace
		fitfilt_cat[i] = 0
		else //bad trace
		fitfilt_cat[i] = 1
		stick_cat[i] = NaN
		stickquant[i] = NaN
		endif
	endfor
end

function filter_accepted()
	nvar num_acc
	make/n=(num_acc)/o filt_val = NaN // this is the final value determining what traces are included
	wave stick_cat, fitfilt_cat
	variable i = 0
	for(;i<num_acc;i+=1)
		if(stick_cat[i] == 0 & fitfilt_cat[i] ==0) //good trace, passing all filters
			filt_val[i] = 0
		else //bad trace, for one or more reasons
			filt_val[i] = 1
		endif
	endfor
	
end

//display trimmed version of the accepted traces
function dispTrimAcc([filter])
	variable filter
	
	if(paramisdefault(filter))
		filter = 0
	endif	
	
	nvar num_acc
	wave acc_traces,acc_bound1, acc_bound2
	wave/t copied_ttc
	string windowname
	if(filter==0)
		 windowname = "Trimmed Accepted Traces"
	elseif(filter==1)
		 windowname = "Filtered Trimmed Accepted Traces"
		wave filt_val
	endif
	DoWindow/K TrimAcc
	string tempname
	variable buildindex
	if(filter==0)
		buildindex = 0
	elseif(filter==1)
		//find first acceptable trace
		findlevel/q filt_val, 0
		buildindex = V_LevelX
	endif
	tempname = copied_ttc[buildindex]
	Display/K=1/n=windowname $tempname[acc_bound1[buildindex],acc_bound2[buildindex]]
	DoWindow/C TrimAcc
	ModifyGraph zero(left)=1
	duplicate/o $tempname, temptr
	variable yoff = temptr[acc_bound1[buildindex]]
	variable xoff = -acc_bound1[buildindex]
	ModifyGraph offset($tempname) = {xoff,-yoff}
	//print "xoff original:" + num2str(xoff)
	//variable xoff =acc_bound2[0]
	if(filter==0)
		xoff = find_xoffset(0,xoff,1)
	elseif(filter==1)
		findlevel/q/r=[buildindex+1,num_acc] filt_val,0
		xoff = find_xoffset(0,xoff,V_LevelX)
	endif
	//print "xoff after:" +num2str(xoff)
	variable i
	
	for(i=1;i<num_acc;i+=1)
		if(filter==0 || filt_val[i] ==0) //filtering off or good trace
			tempname = copied_ttc[i]
			AppendToGraph $tempname[acc_bound1[i],acc_bound2[i]]
			duplicate/o $tempname, temptr
			yoff = temptr[acc_bound1[i]]
			//modifygraph offset($tempname) = {xoff-acc_bound1[i],-yoff}
			modifygraph offset($tempname) = {xoff,-yoff}
			//xoff += acc_bound2[i] - acc_bound1[i] //older version: simple offset
			if(filter==0)
				xoff = find_xoffset(i,xoff, i+1)
			elseif(filter==1)
				findlevel/q/r=[i+1,num_acc] filt_val,0
				xoff = find_xoffset(i,xoff,V_LevelX) //doesn't really work for the last trace. position might be a little odd. can manually adjust one.
			endif
		endif
	endfor
//to do: it would be nice to update xoff in a more intelligent way, to get tighter packing without overlap. Here, no overlap but loosely packed traces.
	Label left "angle (degree)"
	Label bottom "frames"
	
	Killwaves temptr
end
//
//determine min x_off value that gives no intersection in displayed trimmed traces
//tr1,tr2 are accepted trace numbers. not raw. --> tr2 = tr1+1
//intersect: if trace 2 is ever higher than trace 1 at the same adjusted x value, then they must have intersected
//return acceptable xoff for trace 2. seems to work fine. generally not as pretty as manually adjusting
//6/22/16have to edit for filtered displaying only filtered traces. explicitly pass tr2
Function find_xoffset(tr1,xoff1, tr2)
	variable tr1,xoff1, tr2
	variable yoff1,yoff2
	wave acc_bound1,acc_bound2
	wave/t copied_ttc
	//variable tr2 = tr1+1
	
	variable stpt1 = xoff1 + acc_bound1[tr1] //plotted starting point of previous trace.
	variable off_step = 200 //how many frames to jump the offset at a time
	variable stpt2 = stpt1 + off_step  // proposed starting point for next trace
	variable xoff2 = stpt2 - acc_bound1[tr2] //offset to make that happen.
	make/n=(acc_bound2[	tr1] - acc_bound1[tr1]+1)/o tr1temp
	make/n=(acc_bound2[	tr2] - acc_bound1[tr2]+1)/o tr2temp
	
	//print tr1
	
	duplicate/o $copied_ttc[tr1], ttemp
	tr1temp =ttemp[p+acc_bound1[tr1]]
	duplicate/o $copied_ttc[tr2], ttemp
	tr2temp = ttemp[p+acc_bound1[tr2]]
	
	yoff1 = tr1temp[0]
	yoff2 = tr2temp[0]
	
	//apply y offsets.
	tr1temp -=yoff1
	tr2temp -=yoff2
	//set scale so that the trace starts at the 'right' place
	setscale/p x stpt1,1, tr1temp
	setscale/p x stpt2,1,tr2temp
	//what range do we need to examine? from x = 0 to
	//variable xmax = min(dimsize(tr1temp,0),dimsize(tr2temp,0)) //shorter trace.
	variable xmax,xmin
	//note both temp traces start at zero now.
	variable done = 0 //have we found the required offset?
	variable i
	variable fail = 0
	do
		//print xoff2
		xmax = min(stpt2 + dimsize(tr2temp,0), stpt1 + dimsize(tr1temp,0))
		xmin = max(stpt1,stpt2)
		for(i=xmin;i<xmax;i+=1)
			if(tr1temp(i) < tr2temp(i)+360) //compare using x values, not index.
			//print xoff2
				//print i
				//print tr1temp[(i+xoff2
				fail =1
				break //this xoff is no good.
			endif
		endfor
		
		if(fail==1)
			stpt2 +=off_step
			setscale/p x stpt2,1,tr2temp
			fail =0
		else
			done = 1
		endif
	while(done ==0)
	//xoff2 += off_step //add a bit of extra spacing.
	//this can actually cause conflict if there is backtracking!! better to just use a bigger jump (off_step)
	xoff2 = stpt2 - acc_bound1[tr2]
	return xoff2
END
/////////////////////////////////
/////////////////////////////////



////////////////
///////////////PAUSE FINDER FUNCTIONS
//new and improved pause finder. with more careful consideration of filtering. BDA 4/4/16-4/5/16
Function findPause3(tr)
	string tr //note tr is string with trace to be analyzed, not accepted number. (eg ttc or ttc_m#_####)
	//because this can be run on traces before they are accepted -- just doesn't permanently keep the results
	
	///parameters are now set in begin(). can be reset manually anytime.
//	variable smth_window = 111//should be odd
//	variable rate_thr = 2.5
//	variable rev_rate_thr = -rate_thr
//	variable start_thr = 35 	//extend pauses forward/backwards until the values are >start_thr different than the mean
//	variable p_thr =35//merge pauses which are within p_thr of each other.
//	variable fr_back =2   //allow fr_back - 1 frames out of 34 deg range while extending pauses forward/backward
//	variable min_break_dur = 30
	///////
	nvar smth_window, rate_thr, rev_rate_thr, start_thr, p_thr, fr_back, min_break_dur, smth_window2
	
	wave acc_bound1,acc_bound2
	wave/t copied_ttc
	duplicate/o $tr,ttc_fp
	smooth/M=(NaN) 3, ttc_fp //replace NaN with 3 point median. don't change non-NaN. fails if too many NaNs.
	//(most traces have no NaN's during translocation)
	
	variable trstart,trstop,dist,rate, n_fr
	wave ttc
	
	//where is the translocation phase?
	if(stringmatch(tr,"ttc")) //use cursors from Angle_window
		trstart =pcsr(A,"Angle_window")
		trstop=pcsr(B,"Angle_window")
		 n_fr = dimsize(ttc,0)
	else //trace is from accepted traces. use acc_bound1 and acc_bound2
		FindValue/TEXT=tr/txop=4 copied_ttc
		variable acc_pos = V_value
		trstart = acc_bound1[acc_pos]
		trstop = acc_bound2[acc_pos]
		 n_fr = dimsize($copied_ttc[acc_pos],0)
		// print n_fr
	endif
	
	//dist = ttc_fp[trstop]-ttc_fp[trstart]
	//rate = dist / (trstop-trstart)
	
	make/n=(n_fr)/o pause_ID
	make/n=(n_fr)/o fit_rate
	pause_ID[] = 0
	pause_ID[0,trstart-1] = -1
	pause_ID[trstop+1,n_fr-1] = -1
	
	//filter and differentiate
	duplicate/o ttc_fp,ttc_sm
	//Smooth/S=2 smth_window, ttc_sm//SG 2nd order filter
	//Smooth/b smth_window, ttc_sm //boxcar average
	Smooth smth_window, ttc_sm //binomial weighted average -- see notes 20160404 -- USE THIS
	Differentiate ttc_sm/D=ttc_sm_DIF //using central differences
	
	//assign initial value to each frame
	//-1 : outside translocation window (already set); 0: run, 1: pause, 2: backtrack
	variable i = 0
	for(i=trstart;i<=trstop; i+=1)
		if(ttc_sm_DIF[i] < rev_rate_thr) //backtrack
			pause_ID[i] = 2
		elseif(ttc_sm_DIF[i] < rate_thr) //pause
			pause_ID[i] = 1
		else //run //note run is the default if DIF is NAN because the filtering gave 
		//NAN at the beginning and end because there weren't enough points outside the translocation window
			pause_ID[i] = 0
		endif
	endfor

	
	//apply a second round of pause finding with more smoothing, if turned on.
	//avoids breaking up long pauses due to occassional bigger fluctuations
	if(smth_window2 > 0)
		duplicate/o ttc_fp, ttc_sm2
		Smooth smth_window2, ttc_sm2
		Differentiate ttc_sm2/D=ttc_sm2_DIF
		//leave all old pause values alone. but if previously assigned as a run or backtrack and now a pause, reassign as a pause.
		//that is, only change values if they are becoming a pause.
		for(i=trstart;i<=trstop; i+=1)
			if(abs(ttc_sm2_DIF[i]) < rate_thr) //pause
				pause_ID[i] = 1
			endif
		endfor
	endif
	
	/////////clean up results	
	//merge pauses which are essentially at the same angular position - avoid letting brief fluctuations split up a pause
	//for each pause, compare its average position to the previous pause's average position. if they are less than
	//a threshold value, fill in the values between them as a pause. 
	variable in_pause =0 //0 while in a run, 1 while in a pause
	variable pause_start = 0
	variable pause_end = 0
	variable pause_pos
	//variable pause_start,pause_end, pause_pos
	variable prev_pause_pos, prev_pause_end,prev_pause_start
	variable run_pos
	prev_pause_pos = -100000 //initialize to a value that will never match anything.
	for(i=trstart;i<trstop+1;i+=1)
		//wait for the pause to start
		if(pause_ID[i]==1 && in_pause==0)
			in_pause =1 //now in a pause
			//print "pause started at " + num2str(i)
			pause_start = i
		endif
		if((pause_ID[i] != 1  || i == trstop) && in_pause ==1)//end of a pause
			in_pause = 0
			pause_end = i-1
			//print "pause ended at " +num2str(i-1)
			pause_pos = mean(ttc_fp,pause_start,pause_end)
			run_pos = mean(ttc_fp,prev_pause_end+1,pause_start-1)
			//compare to previous pause position and run position
			if((abs(pause_pos-prev_pause_pos) < p_thr) && ((abs(run_pos - pause_pos)<p_thr)  || (abs(run_pos-prev_pause_pos) <p_thr)))//merge pauses
				//print "merging pauses: " +num2str(prev_pause_start) +":" +num2str(pause_end)
				pause_ID[prev_pause_end,pause_start] = 1
				//reset pause variables for the new extended pause
				prev_pause_end = pause_end
				//don't change prev_pause_start
				prev_pause_pos = mean(ttc_fp,prev_pause_start,prev_pause_end)
			else //set current pause to prev pause and continue on. no merge.
				prev_pause_pos = pause_pos
				prev_pause_end = pause_end
				prev_pause_start = pause_start
			endif
		endif
	endfor

//extend pauses backwards until the values are >start_thr different than the mean
	//extend pauses forward in the same way
	//allow 'fr_back'-1 number of outliers while still continuing the pause back or forward
	//this is the same as findPauses2(). corrected a minor problem 20160609.
	in_pause = 0
	pause_start = 0
	pause_end = 0
	pause_pos = 0
	variable j
	variable cont
	variable p_end
	pause_ID[0,trstart -1] = -1
	pause_ID[trstop+1,n_fr] = -1
	for(i=trstart;i<trstop+1;i+=1)
		//wait for the pause to start
		if(pause_ID[i]==1 && in_pause==0)
			in_pause =1 //now in a pause
			//print "pause started at " + num2str(i)
			pause_start = i
		endif
		//wait for pause to end
		if((pause_ID[i] != 1 || i ==trstop) && in_pause ==1)
			pause_end = i-1
			p_end = pause_end
			//print "pause ended at " + num2str(i-1)
			pause_pos = mean(ttc_fp,pause_start,pause_end)
			//go backwards until we've left the pause. pause_pos is not updated.
			j = pause_start -1
			cont = fr_back
			do
				//	print "pause position: " +num2str(pause_pos)
				if(abs(ttc_fp[j]-pause_pos)<=start_thr) //continue the pause backwards
					//print j
					if(pause_ID[j-1] !=1 || pause_ID[j-2] !=1)//don't continue the pause backwards if that will merge pauses. unless that pause is only a single point.
						pause_ID[j,j+fr_back-1]  = 1	
						j -= 1//go back another step.
						cont = fr_back
					else //if we would be merging pauses, stop here. will check if these pauses should be merged later.
						cont =0
					endif
				else 
					cont -= 1 //if appropriate, check the next frame back but note that we've skipped one.
					j -= 1
				endif
				//			if(pause_ID[j-1] ==1) //reached the previous pause
				//				cont = 0
				//			endif
			while(cont>0)
			
			//go forward until we've left the pause
			j = pause_end+1
			cont =fr_back
			do
				//	print j
				if(abs(ttc_fp[j]-pause_pos)<=start_thr) //continue the pause forward
					//	print j
					if(pause_ID[j+1] != 1 || pause_ID[j+2] !=1) //don't continue the pause forward if that will merge pauses. unless that pause is really short
						//		print j
						pause_ID[j-fr_back+cont,j]  =1
						p_end = j // update pause end.
						j += 1//go forward another step.
						cont = fr_back
					else //if we would be merging pauses, stop here. will check if these pauses should be merged later.
						cont =0
					endif
				else
					cont -=1
					j+=1
				endif
				//		if(pause_ID[j+1] ==1) //reached the next pause
				//			cont = 0
				//		endif
			while(cont>0)		
			in_pause = 0
			i = p_end+1 //jump ahead to avoid re-extending this pause
			
		endif
	endfor
	
	//remove short 'breaks' from forward motion, 
	//breaks = pauses and backtracking, combined together if continuous.
	variable in_break = 0
	variable break_start = 0
	variable break_end = 0
	for(i=trstart;i<trstop+1; i+=1)
		//wait for break to start
		if(pause_ID[i] !=0 && in_break ==0)
			in_break =1
			break_start = i
		//	print "break starts : " + num2str(break_start)
		endif	
		//wait for break to end
		if((pause_ID[i] ==0 || i==trstop) && in_break ==1)
			break_end = i-1
			//print "break ends : " +num2str(break_end)
			if((break_end-break_start+1) < min_break_dur)
				pause_ID[break_start,break_end] = 0
				//print "set to run"
			endif
			in_break =0
		endif
	endfor

	//note if this is run through accPauses(), pauseID is saved by that function right after this runs.
END




///find pauses for all accepted traces between acc_bound1 and acc_bound2
//just need to save a copy of pause_ID generated by each iteration of findPauses.
//copy into pauseID_m###_####
function accPauses()
	wave/t copied_ttc
	wave acc_traces, acc_bound1,acc_bound2
	nvar num_acc

	make/n=(num_acc)/o/t copied_pauseID
	variable i
	string cpause, name_temp


	//for storing information about pauses and runs. want this to refresh each time this is run.
	//these hold results for all accepted traces, regardless of whether they pass the sticking test
	variable/g num_pause = 0
	variable/g num_run = 0
	variable/g num_break = 0 //breaks from runs - combined pauses and backtracking. redundant, perhaps but useful
	variable/g num_bktr = 0 //backtracking
	variable/g num_slip = 0 //slips
	variable/g num_brkns = 0 //breaks without slips  
	make/o/n=10000 pause_source=NaN, run_source =NaN //save which trace the pause or run comes from. use the acc_traces index, not true pause number
	make/o/n=10000 pause_dur=NaN, pause_disp=NaN, pause_loc =NaN //pause stats
	make/o/n=10000 run_dur=NaN, run_disp=NaN, run_rate=NaN, run_loc=NaN //run_loc is angle where the run starts
	make/o/n=10000 break_source = NaN, break_dur = NaN, break_disp = NaN, break_loc = NaN
	make/o/n=10000 bktr_source = NaN, bktr_dur = NaN, bktr_disp = NaN, bktr_loc = NaN, bktr_rate = NaN	
	make/o/n=10000 pause_angSTD = NaN, pause_radSTD = NaN
	make/o/n=(num_acc) brkfree_rate = NaN // = (acc_deltaAngle / duration - break frames)
	///'slips' = continguous pause+backtracking events (breaks which include at least one backtracking frame)
	make/o/n=10000 slip_source = NaN, slip_dur = NaN, slip_recdur = NaN, slip_disp = NaN, slip_loc = NaN, slip_recdur2 = Nan
	//slip_dur is entire even. slip_recdur is from first backtracking frame until next 'run' -- recovery time
	//slip_disp = how far back does it go. slip_loc = location of initial pause
	make/o/n=10000 brkns_source = NaN, brkns_dur = NaN, brkns_disp = NaN, brkns_loc = NaN //breaks which aren't slips

  
	if(cmpstr(getwavesdatafolder(pause_source,1),"root:")==0) // don't generate lots of windows if this is part of a combined expt.
	edit/k=1 pause_source
	appendtotable pause_loc,pause_dur,pause_disp
	appendtotable run_source,run_loc,run_disp,run_dur,run_rate
	appendtotable bktr_source, bktr_dur, bktr_disp, bktr_loc, bktr_rate
	appendtotable break_source,break_dur, break_disp, break_loc, brkfree_rate
	appendtotable pause_angSTD, pause_radSTD
	endif

	//AccTrim(0) //need this later, if it hasn't already been run

	for(i=0;i<num_acc;i+=1)
		findPause3(copied_ttc[i])
		wave pause_ID
		name_temp = copied_ttc[i]
		cpause = "pauseID_" + name_temp[4,strlen(name_temp)-1]
		duplicate/o pause_ID,$cpause
		duplicate/o/R=[acc_bound1[i],acc_bound2[i]] pause_ID,$(cpause+"p")
		copied_pauseID[i] = cpause
		analyzePause2(i)
	endfor

	redimension/n=(num_pause) pause_source,pause_dur,pause_disp,pause_loc, pause_angSTD, pause_radSTD
	redimension/n=(num_run) run_source, run_dur,run_disp,run_rate,run_loc
	redimension/n=(num_bktr) bktr_source, bktr_dur,bktr_disp,bktr_loc,bktr_rate
	redimension/n=(num_break) break_source, break_dur, break_disp, break_loc
	redimension/n=(num_slip) slip_source, slip_dur, slip_recdur,slip_disp, slip_loc, slip_recdur2
	redimension/n=(num_brkns) brkns_source, brkns_dur, brkns_disp, brkns_loc



	//display trimmed traces, color coded. unless this window already exists
	DoWindow/F PauseAcc
	if(V_flag==0) //window does not exist
	string windowname = "Pauses in Accepted Traces"
	DoWindow/K PauseAcc
	
	//find first acceptable trace
	findlevel/q filt_val, 0
	variable buildindex = V_LevelX
	
	string tempname = copied_ttc[buildindex]
	Display/K=1/n=windowname $tempname[acc_bound1[buildindex],acc_bound2[buildindex]]
	ModifyGraph zColor($tempname)={$(copied_pauseID[buildindex]+"p"),0,2,rainbow,0}
	DoWindow/C PauseAcc
	ModifyGraph zero(left)=1
	duplicate/o $copied_ttc[buildindex], temptr
	variable yoff = temptr[acc_bound1[buildindex]]
	ModifyGraph offset($tempname) = {-acc_bound1[buildindex],-yoff}

	variable xoff = acc_bound2[buildindex] - acc_bound1[buildindex]
	
	wave filt_val
	
	for(i=buildindex+1;i<num_acc;i+=1)
		if(filt_val[i] == 0)
		tempname = copied_ttc[i]
		AppendToGraph $tempname[acc_bound1[i],acc_bound2[i]]
		ModifyGraph zColor($tempname)={$(copied_pauseID[i]+"p"),0,2,rainbow,0}
		duplicate/o $copied_ttc[i], temptr
		yoff = temptr[acc_bound1[i]]
		modifygraph offset($tempname) = {xoff-acc_bound1[i],-yoff}
		xoff += acc_bound2[i] - acc_bound1[i]
		//print xoff
		//print i
		endif
	endfor

	Label left "angle (degree)"
	Label bottom "frames"
	endif
	Killwaves temptr
	
	filt_pause_data()
	proc_rpdata()
END


//filter the pause data based on trace quality filters -- 
//create a set of waves with pause data for good traces
function filt_pause_data()
	//generate stats for traces which were accepted by AccStick algorithm
	nvar num_filt, num_acc
	nvar num_pause, num_run, num_break, num_bktr, num_slip, num_brkns
	wave pause_source, run_source, pause_dur, pause_disp, pause_loc, run_dur, run_disp, run_rate, run_loc
	wave break_source,break_dur, break_disp, break_loc, brkfree_rate, bktr_source, bktr_dur, bktr_disp, bktr_loc, bktr_rate
	wave slip_source, slip_dur, slip_recdur, slip_loc, slip_disp, slip_recdur2
	wave brkns_source, brkns_dur, brkns_loc, brkns_disp
	wave pause_angSTD, pause_radSTD
	//for storing information about pauses and runs. 
	//holds results for traces which pass sticking test only 
	variable/g fnum_pause = 0
	variable/g fnum_run = 0
	variable/g fnum_break = 0 //breaks from runs - combined pauses and backtracking. redundant, perhaps but useful
	variable/g fnum_bktr = 0 //backtracking
	variable/g fnum_slip = 0 //slips
	variable/g fnum_brkns = 0 //breaks that aren't slips
	make/o/n=10000 fpause_source=NaN, frun_source =NaN //save which trace the pause or run comes from. use the acc_traces index, not true pause number
	make/o/n=10000 fpause_dur=NaN, fpause_disp=NaN, fpause_loc =NaN //pause stats
	make/o/n=10000 frun_dur=NaN, frun_disp=NaN, frun_rate=NaN, frun_loc=NaN //run_loc is angle where the run starts
	make/o/n=10000 fbreak_source = NaN, fbreak_dur = NaN, fbreak_disp = NaN, fbreak_loc = NaN
	make/o/n=10000 fbktr_source = NaN, fbktr_dur = NaN, fbktr_disp = NaN, fbktr_loc = NaN, fbktr_rate = NaN
	make/o/n=10000 fpause_angSTD = NaN, fpause_radSTD = NaN
	
	make/o/n=(num_filt) fbrkfree_rate = NaN // = (acc_deltaAngle / duration - break frames)

	make/o/n=10000 fslip_source = NaN, fslip_dur = NaN, fslip_recdur = NaN, fslip_disp = NaN, fslip_loc = NaN, fslip_recdur2
	make/o/n=10000 fbrkns_source = NaN, fbrkns_dur = NaN, fbrkns_disp = NaN, fbrkns_loc = NaN


	//separately go through pauses, runs, break, bktr and add to filtered versions if appropriate...
		wave filt_val
		variable i
		//runs
		for(i=0; i<num_run;i+=1)
			if(filt_val[run_source[i]] ==0) //good trace
				frun_source[fnum_run] = run_source[i]
				frun_dur[fnum_run] = run_dur[i]
				frun_disp[fnum_run] = run_disp[i]
				frun_rate[fnum_run] =  run_rate[i]
				frun_loc[fnum_run] = run_loc[i]
				fnum_run +=1
			endif
		endfor
		//pauses
		for(i=0; i<num_pause;i+=1)
			if(filt_val[pause_source[i]] ==0) //good trace
				fpause_source[fnum_pause] = pause_source[i]
				fpause_dur[fnum_pause] =pause_dur[i]
				fpause_disp[fnum_pause] = pause_disp[i]
				fpause_loc[fnum_pause] = pause_loc[i]
				fpause_angSTD[fnum_pause] = pause_angSTD[i]
				fpause_radSTD[fnum_pause] = pause_radSTD[i]
				fnum_pause +=1
			endif
		endfor
		//backtracking
		for(i=0; i<num_bktr;i+=1)
			if(filt_val[bktr_source[i]] ==0) //good trace
				fbktr_source[fnum_bktr] = bktr_source[i]
				fbktr_dur[fnum_bktr] = bktr_dur[i]
				fbktr_disp[fnum_bktr] = bktr_disp[i]
				fbktr_loc[fnum_bktr] = bktr_loc[i]
				fbktr_rate[fnum_bktr] = bktr_rate[i]
				fnum_bktr +=1
			endif
		endfor
		//'breaks'
		for(i=0; i<num_break;i+=1)
			if(filt_val[break_source[i]] ==0) //good trace
				fbreak_source[fnum_break] = break_source[i]
				fbreak_dur[fnum_break] = break_dur[i]
				fbreak_disp[fnum_break] = break_disp[i]
				fbreak_loc[fnum_break] = break_loc[i]
				fnum_break +=1
			endif
		endfor
		
		//'slips'
		for(i=0;i<num_slip;i+=1)
			if(filt_val[slip_source[i]] ==0) //good trace
				fslip_source[fnum_slip] = slip_source[i]
				fslip_dur[fnum_slip] = slip_dur[i]
				fslip_recdur[fnum_slip] = slip_recdur[i]
				fslip_disp[fnum_slip] = slip_disp[i]
				fslip_loc[fnum_slip] = slip_loc[i]
				fslip_recdur2[fnum_slip] = slip_recdur2[i]
				fnum_slip +=1
			endif
		endfor
		
		//breaks, no slip
		for(i=0;i<num_brkns;i+=1)
			if(filt_val[brkns_source[i]] ==0) //good trace
				fbrkns_source[fnum_brkns] = brkns_source[i]
				fbrkns_dur[fnum_brkns] = brkns_dur[i]
				fbrkns_disp[fnum_brkns] = brkns_disp[i]
				fbrkns_loc[fnum_brkns] = brkns_loc[i]
				fnum_brkns +=1
			endif
		endfor
		
	redimension/n=(fnum_pause) fpause_source,fpause_dur,fpause_disp,fpause_loc, fpause_angSTD, fpause_radSTD
	redimension/n=(fnum_run) frun_source, frun_dur,frun_disp,frun_rate,frun_loc
	redimension/n=(fnum_bktr) fbktr_source, fbktr_dur,fbktr_disp,fbktr_loc,fbktr_rate
	redimension/n=(fnum_break) fbreak_source, fbreak_dur, fbreak_disp, fbreak_loc
	redimension/n=(fnum_slip) fslip_source, fslip_dur, fslip_recdur, fslip_disp, fslip_loc, fslip_recdur2
	redimension/n=(fnum_brkns) fbrkns_source, fbrkns_dur, fbrkns_disp, fbrkns_loc

	
	
	//finally, deal with break free rate
	variable foo = 0
	for(i=0; i<num_acc; i+=1)
		if(filt_val[i]==0)
			fbrkfree_rate[foo] = brkfree_rate[i]
			foo+=1
		endif
	endfor	
	
	if(cmpstr(getwavesdatafolder(pause_source,1),"root:")==0)
	edit/k=1 fpause_source
	appendtotable fpause_loc,fpause_dur,fpause_disp
	appendtotable frun_source,frun_loc,frun_disp,frun_dur,frun_rate
	appendtotable fbktr_source, fbktr_dur, fbktr_disp, fbktr_loc, fbktr_rate
	appendtotable fbreak_source,fbreak_dur, fbreak_disp, fbreak_loc, fbrkfree_rate
	appendtotable fpause_angSTD, fpause_radSTD
	endif
	
END


//analyze pauses, runs, backtracking, breaks
function analyzePause2(trn)
	variable trn
	wave/t copied_ttc,copied_pauseID, copied_ddc
	wave acc_bound1,acc_bound2
	
	duplicate/o $copied_ttc[trn], ttc_apr
	duplicate/o $copied_pauseID[trn], pID_apr
	duplicate/o $copied_ddc[trn], ddc_apr
	
	smooth/M=(NAN) 3, ttc_apr //replace NaN with 3 pt median. don't change non-NaN
	
	nvar num_run, num_pause,num_bktr,num_break
	wave pause_source,pause_dur,pause_disp,pause_loc
	wave run_source,run_dur,run_disp,run_rate,run_loc
	wave bktr_source,bktr_dur,bktr_disp,bktr_loc,bktr_rate
	wave break_source,break_dur, break_disp, break_loc
	wave brkfree_rate, acc_deltaAngle
	wave pause_angSTD, pause_radSTD
	wave slip_source,slip_dur,slip_recdur, slip_disp, slip_loc, slip_recdur2
	
	variable i 
	variable currphase = pID_apr[acc_bound1[trn]]
	variable phase_counter = 0
		
	//first, go through and pick out boundaries and locations (type, start,end, mean location)
	make/o/n=(1000,4) pick_bnd = NaN
	variable picked = 0
	variable phase_st = acc_bound1[trn]
	for(i=acc_bound1[trn];i<acc_bound2[trn]+2;i+=1)
		if(currphase ==pID_apr[i]) //phase continues. do nothing
			phase_counter +=1
		else //phase over.	
			pick_bnd[picked][0] = currphase
			pick_bnd[picked][1] =phase_st
			pick_bnd[picked][2] = i-1
			pick_bnd[picked][3] = mean(ttc_apr,phase_st,i-1)  - ttc_apr[acc_bound1[trn]]
			
			picked +=1
			currphase = pID_apr[i]
			phase_st = i
		endif
	endfor
	//print pick_bnd[2][0]
	//print picked
	
	variable stloc,endloc
	variable x1,x2
	//then, go through and do calculations for each phase. need the locations of each phase do this
	for(i=0;i<picked;i+=1)
		if(pick_bnd[i][0] ==0) //run
			run_source[num_run] = trn
			run_dur[num_run] = pick_bnd[i][2] - pick_bnd[i][1] +1 //number of frames in run phase
			run_loc[num_run] = ttc_apr[pick_bnd[i][1]] -ttc_apr[acc_bound1[trn]] //start of the run relative to start of trace.
			//displacement: if between two pauses, take the difference of those pauses. 
			//for either end which isn't a pause, use the next frame out from this phase
			if(i==0) //start
				stloc = ttc_apr[pick_bnd[i][1]]  - ttc_apr[acc_bound1[trn]]
			elseif(pick_bnd[i-1][0] ==1) //from a pause
				stloc = pick_bnd[i-1][3] 
			else //not from a pause
				stloc=ttc_apr[pick_bnd[i][1]] - ttc_apr[acc_bound1[trn]]
			endif
			
			if(i==picked) //end
				endloc = ttc_apr[pick_bnd[i][2]] - ttc_apr[acc_bound1[trn]]
			elseif(pick_bnd[i+1][0] ==1) //to a pause
				endloc=pick_bnd[i+1][3]
			else //not to a pause
				endloc=ttc_apr[pick_bnd[i][2]] - ttc_apr[acc_bound1[trn]]
			endif
			run_disp[num_run] = endloc-stloc
			run_rate[num_run] = run_disp[num_run] / run_dur[num_run]
			
			num_run += 1
		elseif(pick_bnd[i][0] == 1) //pause
			pause_source[num_pause] = trn
			pause_dur[num_pause] = pick_bnd[i][2] - pick_bnd[i][1] + 1
			pause_loc[num_pause] = pick_bnd[i][3]
			pause_disp[num_pause] = ttc_apr[pick_bnd[i][2]] - ttc_apr[pick_bnd[i][1]] //hopefully this is nearly 0.
			
			x1 = pnt2x(ttc_apr,pick_bnd[i][1])
			x2 = pnt2x(ttc_apr,pick_bnd[i][2])
			pause_angSTD[num_pause] = variance(ttc_apr,x1,x2)^0.5
			pause_radSTD[num_pause] = variance(ddc_apr,x1,x2)^0.5
			
			num_pause +=1
		else //backtrack
			bktr_source[num_bktr] = trn
			bktr_dur[num_bktr] = pick_bnd[i][2] - pick_bnd[i][1] +1//not especially reliable - short --> depends strongly on where exactly the bounds are
			bktr_loc[num_bktr] = ttc_apr[pick_bnd[i][1]]  - ttc_apr[acc_bound1[trn]]//first frame of backtracking
			//displacement: usually, backtracking surrounded by pauses. use their locations where possible. like with runs above.
			if(i==0) //start
				stloc = ttc_apr[pick_bnd[i][1]] - ttc_apr[acc_bound1[trn]]
			elseif(pick_bnd[i-1][0] ==1) //from a pause
				stloc = pick_bnd[i-1][3] 
			else //not from a pause
				stloc=ttc_apr[pick_bnd[i][1]] - ttc_apr[acc_bound1[trn]]
			endif
			
			if(i==picked) //end
				endloc = ttc_apr[pick_bnd[i][2]] - ttc_apr[acc_bound1[trn]]
			elseif(pick_bnd[i+1][0] ==1) //to a pause
				endloc=pick_bnd[i+1][3]
			else //not to a pause
				endloc=ttc_apr[pick_bnd[i][2]] - ttc_apr[acc_bound1[trn]]
			endif
			
			bktr_disp[num_bktr] = endloc-stloc // should be negative!
			bktr_rate[num_bktr] = bktr_disp[num_bktr] / bktr_dur[num_bktr] // in general, not very reliable - bounds typically somewhat uncertain 
			num_bktr +=1
		endif
	endfor
	
	//finally, have to deal with 'break' category. merge pauses and backtracks.
		//first, go through and pick out boundaries and locations (type, start,end, mean location)
	make/o/n=(1000,4) pick_bnd = NaN
	picked = 0
	phase_st = acc_bound1[trn]
	currphase = pID_apr[acc_bound1[trn]]
	variable is_slip = 0
	nvar num_brkns
	wave brkns_source, brkns_dur, brkns_disp, brkns_loc
	
	for(i=acc_bound1[trn];i<acc_bound2[trn]+2;i+=1)
		if((currphase == 0 && pID_apr[i] == 0) || (currphase !=0 && pID_apr[i] !=0)) //phase continues. do nothing
			phase_counter +=1
		else //phase over.	
			if(currphase ==0)
			pick_bnd[picked][0] = 0 //run
			else
			pick_bnd[picked][0] = 1 //pause or backtrack
			endif
			
			pick_bnd[picked][1] =phase_st
			pick_bnd[picked][2] = i-1
			pick_bnd[picked][3] = mean(ttc_apr,phase_st,i-1)  - ttc_apr[acc_bound1[trn]]
			
			//is this not a run?
			if(currphase == 1)
				//are any of these frames backtracking?
				wavestats/q/r=(phase_st,i-1) pID_apr
				if(V_max == 2) //there is a backtracking frame
					is_slip = analyzeslip(trn, phase_st,i-1) // 0 if it was not accepted as a slip, 1 if it was.
				else
					is_slip = 0 //definitely not a slip if no backtracking frames.
				endif
				
				if(is_slip ==0)
					brkns_source[num_brkns] = trn
					brkns_dur[num_brkns] = (i-1) - phase_st + 1
					brkns_disp[num_brkns] = ttc_apr[i-1] - ttc_apr[phase_st]
					brkns_loc[num_brkns] = ttc_apr[phase_st] - ttc_apr[acc_bound1[trn]]
					num_brkns +=1
				endif
			endif
			
			picked +=1
			currphase = pID_apr[i]
			phase_st = i
		endif
	endfor
	//print picked
	//calculations for each break phase. don't want to redo run phase calculations.
	variable tot_break = 0
	for(i=0;i<picked;i+=1)
		if(pick_bnd[i][0] ==1)
			break_source[num_break] = trn
			break_dur[num_break] = pick_bnd[i][2] - pick_bnd[i][1] + 1 //number frames in the break
			break_disp[num_break] = ttc_apr[pick_bnd[i][2]] - ttc_apr[pick_bnd[i][1]] //probably not very useful
			break_loc[num_break] = ttc_apr[pick_bnd[i][1]] - ttc_apr[acc_bound1[trn]]
			tot_break += break_dur[num_break]
			num_break+=1
		endif
	endfor

	//rate with break frames removed
	brkfree_rate[trn]  = acc_deltaAngle[trn] / (acc_bound2[trn] - acc_bound1[trn] + 1 - tot_break)
END

//analyze a slip for the trace currently being analyzed by analyzePause2
function analyzeslip(trn, startfr, stopfr)
	variable trn, startfr, stopfr
	nvar num_slip
	wave ttc_apr, pID_apr
	variable initial_pos
	variable final_pos
	variable slipstart //frame of first backtracking
	variable slipend
	wave slip_source, slip_dur, slip_recdur, slip_disp, slip_loc, acc_bound1, slip_recdur2
	
	//find initial backtrack frame
	findlevel/p/r=[startfr,stopfr]/q pID_apr, 2
	slipstart  = V_LevelX
	if(V_flag ==1)
		print "SLIP DETECTION ERROR IN ANALYZESLIP()" //this never happens.
	endif
	//find final backtrack frame
	findlevel/p/r=[stopfr,startfr]/q pID_apr,2
	slipend = V_LevelX
	
	//find initial_pos. initial pause is everything before slipstart since there can't be any runs or backtracks there
	initial_pos = mean(ttc_apr, startfr, slipstart-1) - ttc_apr[acc_bound1[trn]]
	//find final_pos. this is everything after the last backtrack frame.
	findlevel/p/r=[stopfr,startfr]/q pID_apr, 2
	if(V_LevelX < stopfr)
		final_pos = mean(ttc_apr, V_LevelX+1, stopfr) - ttc_apr[acc_bound1[trn]]
	else //break ends with backtracking frame. very rare. just use last frame.
		final_pos = ttc_apr(stopfr)
	endif
	
	nvar slip_thr //threshold, in degrees, to call it a real slip
	if((final_pos-initial_pos) < -slip_thr)
		slip_source[num_slip] = trn
		slip_loc[num_slip] = initial_pos
		slip_dur[num_slip] = stopfr-startfr+1
		slip_recdur[num_slip] = stopfr - slipstart
		slip_disp[num_slip] = final_pos - initial_pos
		slip_recdur2[num_slip] = stopfr - slipend
	
		num_slip +=1
		return 1
	else	
		return 0 //not a slip.
	endif
end


//additional processing of run/pause data
function proc_rpdata()
	//allow considering run_rate only if run is sufficiently long - short runs have unreliable rates. (pause boundary calling becomes very important!)
	wave run_dur, run_rate
	nvar num_run
	//variable/g run_dur_thr = 20  //moved to begin()
	nvar run_dur_thr
	variable i
	//print "Long run duration threshold: " +num2str(run_dur_thr)
	
	make/n=(num_run)/o longrun_rate
	make/n=(num_run)/o shortrun_rate
	variable num_longrun =0
	variable num_shortrun = 0
	
	for(i=0;i<num_run+1;i+=1)
		if(run_dur[i] >= run_dur_thr)
			longrun_rate[num_longrun] = run_rate[i]
			num_longrun +=1
		else
			shortrun_rate[num_shortrun] = run_rate[i]
			num_shortrun += 1
		endif
	endfor
	
	redimension/n=(num_longrun) longrun_rate
	redimension/n=(num_shortrun) shortrun_rate
	
	//repeat for filtered data
	nvar fnum_run
	wave frun_dur, frun_rate
	make/n=(num_run)/o flongrun_rate
	make/n=(num_run)/o fshortrun_rate
	variable fnum_longrun =0
	variable fnum_shortrun = 0
	
	for(i=0;i<fnum_run+1;i+=1)
		if(frun_dur[i] >= run_dur_thr)
			flongrun_rate[fnum_longrun] = frun_rate[i]
			fnum_longrun +=1
		else
			fshortrun_rate[fnum_shortrun] = frun_rate[i]
			fnum_shortrun += 1
		endif
	endfor
	
	redimension/n=(fnum_longrun) flongrun_rate
	redimension/n=(fnum_shortrun) fshortrun_rate



end


/////////////////////////
////////////////////////