import numpy as np; import pylab as plt; 
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
import tktools.arraytools as arraytools
import tktools.tabletools as tabletools
warnings.simplefilter("once")
sys.path.append('/home/samuroff/GREAT-DES/tktools')
logging_level = logging.INFO; logger = logging.getLogger("proc_shear"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import tktools.plotstools as plotstools 
import tktools.fitting as fitting

#-----------------------------------------------------------------------------------------------------#

def plot_scatter(res):

	e1=np.array([]); ne1= np.array([]) 
	for i in range(len(res)):
		e1=np.hstack((e1,res[i]['e1'])) ; ne1= np.hstack((ne1,res[i]['ngmix_e1']))

	plt.errorbar(-1*e1[0:1000], ne1[0:1000], fmt='.')
	plt.show()
	plt.close()


