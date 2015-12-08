from multiprocessing import Pool
import itertools
from ffpdax import ffp_dax

def ffp_daxunpack(details):
	ffp_dax(*details)
	#http://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
	
#Remember: don't be greedy! but I think this'll only let you run four processes
if __name__ == "__main__":
	files = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma_data\movie_0003','C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma_data\movie_0003']
	xmls = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma\ORBIT1','C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma\ORBIT1']
	analysis = "C:\\Users\\B\\Documents\\Zhuang Lab\\Analysis Code\\storm-analysis-master\\sma\\ffpdax.py"
	p = Pool(processes = 4)#allow up to four processes
	p.map(ffp_daxunpack, itertools.izip(files, xmls))