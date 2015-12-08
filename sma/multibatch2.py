from apdax import ap_dax
from ffpdax import ffp_dax
import thread
import subprocess
#Remember: don't be greedy!
if __name__ == "__main__":
	files = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma_data\movie_0003','C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma_data\movie_0003']
	xmls = ['C:\Users\B\Documents\Zhuang Lab\Analysis Code\storm-analysis-master\sma\ORBIT1']
	analysis = "C:\\Users\\B\\Documents\\Zhuang Lab\\Analysis Code\\storm-analysis-master\\sma\\ffpdax.py"
	proc = subprocess.Popen(['python',analysis,files[0],xmls[0]])
	proc = subprocess.Popen(['python',analysis,files[1],xmls[0]])
	#proc = subprocess.Popen(['python',analysis,files[2],xmls[0]])
	#proc = subprocess.Popen(['python',analysis,files[3],xmls[0]])

