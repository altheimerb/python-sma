#TO DO: run multiple sma analyses in parallel. Not as simple as I thought! Giving error:  
#PyEval_RestoreThread: NULL tstate. Hitting debug option gives more: unhandled win32 exception occured in python.exe
#I believe the problem is related to how resources are handled during threading?
#for now, use this as a simple batch file
#Remember: Don't be greedy!
#from multiprocessing import Process
from ffpdax import ffp_dax
from apdax import ap_dax

if __name__ == "__main__":
	#files = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma_data\movie_0003']
	#xmls = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma\ORBIT1']
#	p = Process(target=ffp_dax,args=(files[0],xmls[0]))
#	p.start()
#	p.join()

	#ffp_dax(PUTFILEHERE,PUTSETTINGSHERE)
	
	print 'all done!'