import numpy as np

def apply_cuts(res, quantity, mode, truth):

	#print res.dtype
	print 'Applying im3shape recommended cuts:'
	sel_string= "sel= "
	sel_string+= "(res['error_flag']==0) "
	if quantity!= 'snr_w':		sel_string+= "& (res['snr']>15) "
	if quantity!= 'rgpp_rp':	sel_string+= "& (res['mean_rgpp_rp']>1.2) "
	sel_string+= "& (np.invert( np.isnan(res['mean_psf_e1_sky']) ) ) & (np.invert( np.isnan(res['mean_psf_e2_sky']) ) ) "
	

	print sel_string
	exec sel_string
	selected_res= res[sel]

	print 'Cut %f percent of %d' %( 100.0*float(len(res)-len(selected_res))/len(res), len(res) )

	if truth != None:
		truth= truth[sel]

	return selected_res, truth 

