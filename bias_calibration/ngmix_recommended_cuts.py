import numpy as np

def apply_cuts(res, quantity, mode, truth):

#    0.4 < exp_arate < 0.6
## & (np.exp(res['log_T_r'])/res['psf_T']>0.15)
	#if mode=='ngmix':
	print 'Applying ngmix recommended cuts:'
	sel_string= "sel=(res['flags']==0) "
	if (quantity!= 'snr_r') and (quantity!= 'snr_w'):	sel_string+=         "& (res['s2n_r']>15) "  
	if mode=='im3shape':
		if quantity!= 'g1_sens':	sel_string+= "& (res['g1_sens']>0.0) " 
		if quantity!= 'g2_sens':	sel_string+= "& (res['g2_sens']>0.0) "
		if quantity!= 'snr_w':		sel_string+= "& (res['snr']>15) "
		if quantity!= 'rgpp_rp':	sel_string+= "& (res['mean_rgpp_rp']>15) "
		sel_string+= "& (np.inverse( np.isnan(res['psf_e1']) ) ) & (np.inverse( np.isnan(res['psf_e2']) ) )"
	elif mode=='ngmix':
		if quantity!= 'g1_sens':	sel_string+= "& (np.transpose(res['g_sens'])[0]>0.0) " 
		if quantity!= 'g2_sens':	sel_string+= "& (np.transpose(res['g_sens'])[1]>0.0) "
		if quantity!= 'T_r':	sel_string+= 	     "&  (np.exp(res['log_T_r'])/res['psf_T'] > 0.15) "
	
	print sel_string
	exec sel_string
	selected_res= res[sel]

	print 'Cut %f percent of %d' %( 100.0*float(len(res)-len(selected_res))/len(res), len(res) )

	if truth != None:
		truth= truth[sel]

	return selected_res, truth 

#	elif mode=='im3shape':
#		return res
