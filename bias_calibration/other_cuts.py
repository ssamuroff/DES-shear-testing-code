import numpy as np; import pylab as plt; 
import pyfits
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
import tktools.arraytools as arraytools
import tktools.tabletools as tabletools
import ngmix_recommended_cuts
warnings.simplefilter("once")
sys.path.append('/home/samuroff/GREAT-DES/tktools')
logging_level = logging.INFO; logger = logging.getLogger("proc_shear"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import tktools.plotstools as plotstools 
import tktools.fitting as fitting

#-----------------------------------------------------------------------------------------------------#

def get_bins(quantity):

	if quantity== 'Tsnr':	return np.array([0,1,1.5,2,2.5,3,3.5,4,5,6,7,8,10,12,14,16,18,20])
	if quantity== 'T_r':	return np.array([0,0.1,0.2,0.3,0.4,0.5,0.6, 0.75,0.9,1,1.25,1.5])
	if (quantity== 'snr_r') or (quantity== 'snr_w'):	return np.array([10,11,12,14,16,18,20,25,30,50,80,200])
	if quantity== 'g1_sens':	return np.array([-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
	if quantity== 'g2_sens':	return np.array([-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

#-----------------------------------------------------------------------------------------------------#

def get_label(quantity,mode):

	if quantity== 'Tsnr':	return 'T_s2n_r'
	if quantity== 'T_r':	return 'log_T_r'
	if quantity== 'snr_r' and mode=='ngmix':	return 's2n_r'
	if quantity== 'snr_r' and mode=='im3shape':	return 'round_snr'
	if quantity== 'snr_w' and mode=='ngmix':	return 's2n_w'
	if quantity== 'snr_w' and mode=='im3shape':        return 'snr'
	if quantity== 'g1_sens':	return 'g_sens'
	if quantity== 'g2_sens':	return 'g_sens'


#-----------------------------------------------------------------------------------------------------#

def get_xlim(quantity):

	if quantity== 'Tsnr':	return 0, 20
	if quantity== 'T_r':	return 0, 1.51
	if (quantity== 'snr_r') or (quantity== 'snr_w'):	return 0, 61
	if quantity== 'g1_sens' or quantity== 'g2_sens':	return -0.1, 1.1

#-----------------------------------------------------------------------------------------------------#

def get_fiducial_cut(quantity):

	if quantity== 'Tsnr':	return 2
	if quantity== 'T_r':	return 0
	if quantity== 'snr_r':	return 15
	if quantity== 'snr_w':	return 20
	if quantity== 'g1_sens' or quantity== 'g2_sens':	return 0

#-----------------------------------------------------------------------------------------------------#

def get_selection_string(i,lab_cut,quantity,mode, limit1, limit2):

	if limit2== -1.009804123:
		if lab_cut== 'log_T_r':	return "sel= (np.exp(res[%d]['%s'])>%f)" %(i, lab_cut, limit1)
		elif lab_cut== 'g_sens' and quantity== 'g1_sens' and mode=='ngmix':	return "sel= (np.transpose(res[%d]['%s'])[0]>%f)" %(i, lab_cut, limit1)
		elif lab_cut== 'g_sens' and quantity== 'g2_sens' and mode=='ngmix':	return "sel= (np.transpose(res[%d]['%s'])[1]>%f)" %(i, lab_cut, limit1)
		elif lab_cut== 'g_sens' and quantity== 'g1_sens' and mode=='im3shape':	return "sel= (res[%d]['%s']>%f)" %(i, lab_cut, limit1)
		elif lab_cut== 'g_sens' and quantity== 'g2_sens' and mode=='im3shape':	return "sel= (res[%d]['%s']>%f)" %(i, lab_cut, limit1)
		else:	return "sel= (res[%d]['%s']>%f)" %(i, lab_cut, limit1)	
	elif limit2!= -1.009804123:
		if lab_cut== 'log_T_r':	return "sel= (np.exp(res[%d]['%s'])>%f) & (np.exp(res[%d]['%s'])<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)
		elif lab_cut== 'g_sens' and quantity== 'g1_sens' and mode=='ngmix':	return "sel= (np.transpose(res[%d]['%s'])[0]>%f) & (np.transpose(res[%d]['%s'])[0]<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)
		elif lab_cut== 'g_sens' and quantity== 'g2_sens' and mode=='ngmix':	return "sel= (np.transpose(res[%d]['%s'])[1]>%f) & (np.transpose(res[%d]['%s'])[1]<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)
		elif lab_cut== 'g_sens' and quantity== 'g1_sens' and mode=='im3shape':	return "sel= (res[%d]['%s']>%f) & (res[%d]['%s']<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)
		elif lab_cut== 'g_sens' and quantity== 'g2_sens' and mode=='im3shape':	return "sel= (res[%d]['%s']>%f) & (res[%d]['%s']<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)
		else:	return "sel= (res[%d]['%s']>%f) & (res[%d]['%s']<%f)" %(i, lab_cut, limit1, i, lab_cut, limit2)	

