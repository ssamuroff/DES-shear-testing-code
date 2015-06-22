#import matplotlib
#matplotlib.use("agg")
import numpy as np; import pylab as plt; 
import cPickle as pickle
import pyfits
from scipy import interpolate 
import  sys, logging, yaml, argparse, time, copy, itertools, tktools, warnings, os, fitsio, pyfits;
import tktools.arraytools as arraytools
import tktools.tabletools as tabletools
import ngmix_recommended_cuts, im3shape_recommended_cuts
import other_cuts
import fitting_routines
warnings.simplefilter("once")
logging_level = logging.INFO; logger = logging.getLogger("proc_shear"); logger.setLevel(logging_level)  
log_formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s   %(message)s", "%Y-%m-%d %H:%M:%S")
stream_handler = logging.StreamHandler(sys.stdout); stream_handler.setFormatter(log_formatter)
if logger.handlers == [] : logger.addHandler(stream_handler); logger.propagate = False
import tktools.plotstools as plotstools 
import tktools.fitting as fitting

#-----------------------------------------------------------------------------------------------------#

def load_ngmix_measurements(load_truth, load_intrinsic_shapes):
	
	# Set up interpolator to calculate bias if needed
	#if args.calibration!='none':	interpolator_m, interpolator_alpha= initialise_calibration(mode='uncalibrated')
	#else:				interpolator_m= None ; interpolator_alpha= None		

	shape_file_string= config['ngmix']['results_string']
	truth_filelist= config['general']['truth_table_string']
	pickles_file_string= config['general']['truth_pickles_string']
	isim= 2; nfiles= 0; nobj= 0; dt=False ; dat= False
	ngmix_res= []
	ngmix_truth= []
	#truth_table=np.array([])
	res_cols= config['ngmix']['res_cols']
	for ig in config['general']['shear_groups_to_use']:
			filename= shape_file_string % (isim, ig)

			#Read the column headings from the first file 
			if not dt:	dtype_res=get_columns(filename) ; dt = True
		
			res= tabletools.loadTable(filename)[0:args.number*10000]

			# Reindex the results tables loaded
			res['number']*=0
			res['number']+= np.arange(0,len(res), 1)

			# Calculate the weights 
			# From Mike's prescription using the trace of the covariance matrix
			# See Jarvis et al (2015)
			if args.weighting==1:			
				c11= np.transpose(res['g_cov'])[0][0] ; c22= np.transpose(res['g_cov'])[1][1]	# Shear covariance matrix diagonal elements
				tr_cov= c11 + c22
				sn= 0.22										# Shape noise per component

				w= 1.0/(2 * sn**2 + tr_cov)
			
			# Split the results table to match the truth table
			l=len(res)/10000
			truth_table=None
			if load_truth and (args.use_truth_tables==1):
				for ifil in range(l):
					filename_truth= truth_filelist % (ifil, ig)
                                	truth= tabletools.loadTable(filename_truth)
                                	
					# Get the column headings from the first file
					if not dat:	dtype_truth= get_columns(filename_truth) ; dat=True

					# Copy the first truth table across
					# Then append subsequent results to it
					if (ifil== 0):
						truth_table= np.copy(truth)
					else:
						if len(truth_table.dtype.names)>len(truth.dtype.names):	
							dt1=set(truth_table.dtype.names) ; dt2=set(truth.dtype.names)
							missing_columns=list(dt1-dt2)
							print 'Warning: truth table columns vary between files. Missing columns: %s' %missing_columns
                                		truth_table= np.hstack((truth_table,truth.astype(truth_table.dtype)))
						
				truth_table= truth_table.astype(dtype=dtype_truth)

				# Check that the lengths of the truth and results arrays match
				if len(truth_table)!=len(res):
                                        	logger.error('Truth table (%d) does not match results catalogue (%d).' %(len(truth_table), len(res)))		

			# Apply selection if required
			if args.selection!= 0:	res, truth_table= select(res,truth=truth_table, num=0, mode='ngmix')
			res= res.astype(dtype_res)	

			# Look for a true shear column in the results file
			# If it's missing, add one
			try:
				print 'Found true shear column g1,2= %2.2f, %2.2f' %(res['g1_true'][0],res['g2_true'][1] )

			except:
				g1= config['general']['shear'][ig][0] ; g2= config['general']['shear'][ig][1] 
				res= tabletools.appendColumn(res,'g1_true', g1)
				res= tabletools.appendColumn(res,'g2_true', g2)
				print 'Added true shear column g1,2= %2.4f, %2.4f' %(-1*res['g1_true'][0],res['g2_true'][1] )

			# Append a zphot column to the results array if required #
			import cPickle as pickle
			ptruth=None
			global zbinning
			if zbinning:
				logger.info('Reading cosmos photo-zs.')
				if args.use_truth_tables==1:	zphot=truth_table['zphot']
				else:	
					ptruth=pickle.load(open(pickles_file_string %ig,'rb'))			
					zphot= ptruth[2][res['number']]

				res=tabletools.appendColumn(res,'zphot',zphot)
				print 'Added zphot column', res['zphot']				
				zupper= args.zphot[1]
				zlower= args.zphot[0]
				if zupper!=0:
					n1=len(res)
					sel= (res['zphot']>zlower) & (res['zphot']<zupper)	
					res=res[sel]
					logger.info('Selected galaxies in bin z=%2.2f-%2.2f',zlower,zupper)	
					n2=len(res)
			
					logger.info('Contains %f percent of objects', 100.0-100.0*(n1-n2)/n1)

			# Calculate the intrinsic shapes and add them to the results array
			if load_intrinsic_shapes:
				if args.use_truth_tables==1:	res, truth= get_intrinsic_shapes(res, truth_table)
				else:
					e1_sheared= ptruth[0][res['number']]
					e2_sheared= ptruth[1][res['number']]
					res=tabletools.appendColumn(res,'e1_sheared',e1_sheared)
					res=tabletools.appendColumn(res,'e2_sheared',e2_sheared)
					print 'Added true sheared e1 column', e1_sheared
				

			# Do the same for the weights
			if args.weighting==1:
				res= tabletools.appendColumn(res, 'w', w[ res['number'] ])
				print 'Added weights column', res['w']

			# Do the same for the bias corrections
			#if args.calibration!= 'none':
			#	m= get_calibration(interpolator_m, res) 
			#	res= tabletools.appendColumn(res,'nbc_m',m)
			#	print 'Added nbc_m column', res['nbc_m']

			#	a= get_calibration(interpolator_alpha, res)
                         #       res= tabletools.appendColumn(res,'nbc_a',a)
                          #      print 'Added nbc_a column', res['nbc_a']
			

			ngmix_res.append(res)
	
			#ngmix_truth.append(truth_table)

			nfiles+= 1
			nobj+= len(res['g'])
		
	ngmix_res= np.array(ngmix_res)
	logger.info('Obtained %2.2e objects from %d datafiles.', nobj, nfiles)

	return ngmix_res, ngmix_truth

#-----------------------------------------------------------------------------------------------------#

def load_im3shape_measurements(load_truth, load_ngmix_columns, load_intrinsic_shapes):
	
	import glob
	shape_filelist= config['im3shape']['results_string']
	truth_filelist= config['general']['truth_table_string']
	
	nfiles= 0; nobj= 0
	dt=False
	im3shape_res= []
	truth_table= [] 
	exec(config['im3shape']['res_cols'])
	exec(config['general']['truth_cols'])

	# Set up interpolator to calculate bias if needed
	if args.calibration!='none':	interpolator_m, interpolator_alpha= initialise_calibration(mode='uncalibrated')		
	else:				interpolator_m= None ; nterpolator_alpha= None	

	for ig in config['general']['shear_groups_to_use']:
		for isim in range(args.number):
			filename= shape_filelist % (isim, ig)
			#Read the column headings from the first file 
			if not dt:	
				try:	dtype= get_columns(filename)
				except:
					logger.error('Warning: file %s is missing or incompatible. Could not load.', filename)
					continue
		
			res= tabletools.loadTable(filename)

			# If required also load the associated truth table			
			if load_truth:	
				filename_truth= truth_filelist % (isim, ig)		
				truth= tabletools.loadTable(filename_truth)

				# Check that the lengths of the truth and results arrays match
				if len(truth)!=len(res):
					print 'Warning: Truth table does not match results catalogue.'
					truth= truth[res['coadd_objects_id']]

			elif not load_truth:	truth_table= None ; truth= None

			res= res.astype(dtype=dtype)

			# Load the ngmix parameter columns if specified
			if load_ngmix_columns:	res= match_ngmix_columns(res,ig,isim)

			# Now make the required selection
			if args.selection!= 0:	res, truth = select(res,truth=truth, num=0, mode='im3shape')

			if load_truth:	
				truth= truth.astype(dtype=dtype_truth)
				truth_table.append(truth)

			# Look for a true shear column in the results file
			# If it's missing, add one
			try:
				print 'Found true shear column g1,2= %2.2f, %2.2f' %(-1*res['g1_true'][0],res['g2_true'][1] )

			except:
				# Note the change of sign for im3shape
				g1= -1*config['general']['shear'][ig][0] ; g2= config['general']['shear'][ig][1] 
				res= tabletools.appendColumn(res,'g1_true', g1)
				res= tabletools.appendColumn(res,'g2_true', g2)
				print 'Added true shear column g1,2= %2.4f, %2.4f' %(-1*res['g1_true'][0],res['g2_true'][1] )	

			# Append a zphot column to the results array if required #
			global zbinning
			if zbinning:
				logger.info('Reading cosmos photo-zs.')
				zphot=truth['zphot']
				print truth['psf_e1']
				print res['mean_psf_e1_sky']
				#plt.errorbar(truth['psf_e1'],-1*res['mean_psf_e1_sky'] , fmt='.')
				#plt.show()


				res=tabletools.appendColumn(res,'zphot',zphot)
				res=tabletools.appendColumn(res,'shape_e1',truth['shape_e1'])
				res=tabletools.appendColumn(res,'shape_e2',truth['shape_e2'])
				print 'Added zphot column', res['zphot']
				zupper= args.zphot[1]
				zlower= args.zphot[0]
				n1=len(res)
				sel= (res['zphot']>zlower) & (res['zphot']<zupper)	
				#res=res[sel]
				logger.info('Selected galaxies in bin z=%2.2f-%2.2f',zlower,zupper)	
				n2=len(res)
			
				logger.info('Contains %f percent of objects', 100.0-100.0*(n1-n2)/n1)

			# Calculate the intrinsic shapes and add them to the results array
			if load_intrinsic_shapes:
				res, truth= get_intrinsic_shapes(res, truth)			

			# Do the same for the bias corrections
			if args.calibration!= 'none':
				m= get_calibration(interpolator_m, res) 
				res=tabletools.appendColumn(res,'nbc_m',m)
				print 'Added nbc_m column', res['nbc_m']
			
				a= get_calibration(interpolator_alpha, res)
                                res=tabletools.appendColumn(res,'nbc_a',a)
                                print 'Added nbc_a column', res['nbc_a']

			im3shape_res.append(res)
			#if isim==0 and ig==0:	tmp=np.copy(res)
			#else:	tmp=np.hstack((tmp,res))
			nfiles+= 1
			nobj+= len(res['coadd_objects_id'])
		
	im3shape_res= np.array(im3shape_res)
	logger.info('Obtained %2.2e objects from %d datafiles.', nobj, nfiles)
	#pyfits.writeto('im3shape_res_v006.fits',tmp,clobber=True)
	return im3shape_res, truth_table
	
#-----------------------------------------------------------------------------------------------------#

def load_matched_measurements(load_truth, load_ngmix_columns, load_intrinsic_shapes, repackage):
	"""
	Loads both im3shape and ngmix results files. Then takes the intersection.  
	Returns two 1d arrays of equal length. 
	"""
	# Load results files
	sel_opt= args.selection
	meth_opt= args.method
	cal_opt= args.calibration
	args.selection= 0
	
	args.calibration= 'none'
	args.method='ngmix'
	ngres,trng= load_ngmix_measurements(load_truth=load_truth, load_intrinsic_shapes=load_intrinsic_shapes)

	args.calibration= cal_opt
	args.method='im3shape'
	imres, trim= load_im3shape_measurements(load_truth=load_truth, load_ngmix_columns=load_ngmix_columns, load_intrinsic_shapes=load_intrinsic_shapes)

	logger.info('Matching im3shape and ngmix galaxies.')
	# Combine the files into a 1d array
	ng_all=np.array([])
	l=0
	for ng_file in ngres:
		if l==0:	ng_all=np.copy(ng_file) ; l=1
		else:		ng_all=np.hstack((ng_all,ng_file))
	
	# Match the im3shape/ngmix galaxies
	ngres=[]
	ntr=10000 
	for ir in range(len(imres)):
		ngres+= [ ng_all[ ir*ntr:(ir+1)*ntr ][imres[ir]['coadd_objects_id']] ]
	ngres=np.array(ngres)

	args.selection= sel_opt
	args.method= meth_opt

	# Apply the required cuts on both sets of results
	if args.selection== 1:	
				sel_string= "sel= (ngres['flags'==0]) & (ngres['flags_r'==0]) & (imres['error_flag']==0) " 
				logger.info('Applying flag cuts: %s', sel_string)
	if args.selection== 2:
		ngsel_string= ngmix_recommended_cuts.get_selection_string(ngres, quantity=args.quantity[0])
		imsel_string= im3shape_recommended_cuts.get_selection_string(imres, quantity=args.quantity[0])
	
		sel_string= imsel_string + ngsel_string.replace('sel=','&')

		logger.info('Applying recommended cuts: %s', sel_string)

	if args.selection!= 0:
		
		sel_ngres=[] ; sel_imres=[]
		n1= 0 ; n2= 0
		for ir in range(len(ngres)):
			exec sel_string
			n1+= len(ngres[ir])
			n2+= len(ngres[ir][sel])
			sel_ngres+= [ ngres[ir][sel] ] ; sel_imres+= [ imres[ir][sel] ]
		ngres=np.array(sel_ngres) ; imres=np.array(sel_imres)

		logger.info('Removed %d objects (%f percent)',n1-n2, 100.0*(n1-n2)/n1)
	else:	logger.info('No cuts applied.')
		 
	# Finally repackage the results into shear groups for consistency with the other loading functions if required
	repackage=False
	if repackage:
		ngres_all= np.copy(ngres)
		imres_all= np.copy(imres)
		ngres= [] ; imres= []
		logger.info('Repackaging files.')
		for i in config['shear_groups_to_use']:
			g1= config['shear'][i][0] ; g2= config['shear'][i][1]
			ngsel= (ngres_all['g1_true']== g1) & (ngres_all['g2_true']== g2)
			imsel= (imres_all['g1_true']== g1) & (imres_all['g2_true']== g2)
			ngres+= [ngres_all[ngsel]]
			imres+= [imres_all[imsel]]
		imres= np.array(imres) ; ngres= np.array(ngres)

	return imres, ngres
	

#-----------------------------------------------------------------------------------------------------#

def load_im3shape_measurements_des():
	
	import glob
	filelist= glob.glob(config['im3shape']['results_string_DES'])
	
	nfiles= 0; nobj= 0
	dt=False
	im3shape_res= np.array([])
	
	for filename in filelist and (nfiles<10):
		#Read the column headings from the first file 
		if not dt:	
			try:	dtype= get_columns(filename)
			except:
				logger.error('Warning: file %s is missing or incompatible. Could not load.', filename)
				continue
		
		res= tabletools.loadTable(filename)

		res= res.astype(dtype=dtype)

		# Now make the required selection
		if args.selection!= 0:	res, truth = select(res,truth=None, num=0, mode='im3shape')

		# Finally append it to the results array 
		if nfiles==0:	im3shape_res= np.copy(res)
		else:	im3shape_res=np.hstack((im3shape_res,res))
		nfiles+= 1
		nobj+= len(res['coadd_objects_id'])
		
	im3shape_res= np.array(im3shape_res)
	logger.info('Obtained %2.2e objects from %d datafiles.', nobj, nfiles)

	return im3shape_res

#-----------------------------------------------------------------------------------------------------#

def match_ngmix_columns(res,ig,isim):
	import cPickle as pickle

	# Load the required columns
	snr_r_filename= config['ngmix']['pickled_columns'] %('snr_r', ig)
	snr_T_filename= config['ngmix']['pickled_columns'] %('snr_T', ig)
	flags_filename= config['ngmix']['pickled_columns'] %('flags', ig)
	g1sens_filename= config['ngmix']['pickled_columns'] %('g1_sens', ig)
	g2sens_filename= config['ngmix']['pickled_columns'] %('g2_sens', ig)

	snr_r= pickle.load(open(snr_r_filename,'rb'))
	snr_T= pickle.load(open(snr_T_filename,'rb'))
	flags= pickle.load(open(flags_filename,'rb'))
	g1_sens= pickle.load(open(g1sens_filename,'rb'))
	g2_sens= pickle.load(open(g2sens_filename,'rb'))


	# Find the appropriate ngmix column corresponding to this file 
	# Then match the individual objects by id
	snr_r= snr_r[isim][res['coadd_objects_id']] ; 
	snr_T= snr_T[isim][res['coadd_objects_id']]
	flags= flags[isim][res['coadd_objects_id']] ; 
	g1_sens= g1_sens[isim][res['coadd_objects_id']] ; g2_sens= g2_sens[isim][res['coadd_objects_id']]	

	# Append the column to the im3shape results table
	res= tabletools.appendColumn(res,'s2n_r',snr_r)
	res= tabletools.appendColumn(res,'T_s2n_r',snr_T)
	res= tabletools.appendColumn(res,'flags',flags)
	res= tabletools.appendColumn(res,'g1_sens',g1_sens)
	res= tabletools.appendColumn(res,'g2_sens',g2_sens)

	print 'snr column', res['snr']
	print 'Matched to snr_r column', res['s2n_r']
	print 'For (%d %d)' %(ig, isim)

	# Finally add true shear columns if required
	try:
		test= res['g1_true']
		res['g1_true']*= -1
	except:
		# Initialise array and read true shear from the config file
				## IMPORTANT ##
		## The im3shape g1 is flipped relative to the ngmix definition ##
		g1_true= np.zeros( np.shape(res['e1']))+ -1.0*config['general']['shear'][ig][0]
		g2_true= np.zeros( np.shape(res['e2']))+ config['general']['shear'][ig][1]

		res= tabletools.appendColumn(res, 'g1_true', g1_true )
		res= tabletools.appendColumn(res, 'g2_true', g2_true )
	
	return res

#-----------------------------------------------------------------------------------------------------#

def load_truth_tables(res=None):
	logger.info('Loading truth tables')
	truth_filelist= config['general']['truth_table_string']
	
	truth= []
	dt=0

	# This loads files in order of shear grouping
	# The first 100 tables correspond to ig=0, the second to ig=1... up to ig=7
	for ig in config['general']['shear_groups_to_use']:
		#if ig > 0:	continue
		tmp_truth= []
		for isim in range(10):
			filename= truth_table_filelist % (isim, ig)
			#Read the column headings from the first file 
			if dt==0:	dtype=get_columns(filename)
		
			truth_table= tabletools.loadTable(filename)
			truth_table= truth_table.astype(dtype=dtype)			

			tmp_truth.append(truth_table)
			truth.append(truth_table)

			nfiles+= 1
			nobj+= len(res['id'])

		match_res_truth(res, tmp_truth, ig)

	return truth

#-----------------------------------------------------------------------------------------------------#

def select(res, truth, num, mode):
	# If using im3shape results the truth tables may also be required
	# These are in separate files so any cuts need to be applied to both the results and the truth arrays
	# Conveniently enough, the ngmix results files contain the true shears for each object
	# Remove flagged objects
	if mode== 'im3shape':	
		flags= 'error_flag'

	if mode== 'ngmix':	
		flags= 'flags'

		
	if truth!=None:	truth= truth[res[flags]==0]

	n0= len(res[flags])
	res= res[res[flags]==0]
	n1= len(res[flags])
	logger.info('Removed %d flagged objects (%f percent) \n',n0-n1, 100*float((n0-n1))/n1)
	
	# Take a random subsample if required
	if num != 0:
		res=np.random.choice(res,num)

	#Impose a selection cut as required	
	exec(config[mode]['cut'])
	if truth!= None:	truth= truth[sel]

	n_nocut=len(res)
	
	res= res[sel]
	n_cut=len(res)
	logger.info('Removed %d objects in selection cuts (%f percent)', n_nocut-n_cut, 100*float(n_nocut-n_cut)/n_nocut )
	
	# If specified, also impose the fiducial cuts for this catalogue
	if (args.selection==2): 
		if args.method=='ngmix':	res, truth= ngmix_recommended_cuts.apply_cuts(res,quantity= args.quantity[0], mode=args.method, truth=truth)
		if args.method=='im3shape':	res, truth= im3shape_recommended_cuts.apply_cuts(res,quantity= args.quantity[0], mode=args.method, truth=truth)

	return res, truth
		
#-----------------------------------------------------------------------------------------------------#

def get_columns(filename):
	fil=pyfits.open(filename)
	fil=fil[1].data
	cols=fil.dtype
	
	logger.debug('Read columns: %s',cols)

	return cols

#-----------------------------------------------------------------------------------------------------#

def get_intrinsic_shapes(res,truth):
	"""
	Computes the intrinsic galaxy shapes for the input selection. 
	Requires access to the COSMOS shape catalogue.
	@return galaxy ellipticity before and after cosmic shear   
	"""

	global match_res

	# Load the COSMOS shape catalogue
	filename_cosmos= config['general']['cosmos_shape_catalogue_string']
	cosmos_cat= tabletools.loadTable(filename_cosmos)

	# Match the truth tables for the specified GREAT-DES selection to the parent COSMOS galaxies 
	logger.info('Obtaining GREAT-DES intrinsic shapes.')
	shape_e1= cosmos_cat[truth['id_cosmos']]['e1']
	shape_e2= cosmos_cat[truth['id_cosmos']]['e2']
	cosmos_cat= cosmos_cat[truth['id_cosmos']]

	# Apply the distortion used by GREAT-DES
	eps= shape_e1+ shape_e2*1j
	q= np.sqrt( (1-np.abs(eps)) / (1+np.abs(eps)) )
	e= (1-q)/(1+q) * np.exp(1j*np.angle(eps))
	e1= e.real ; e2= e.imag

	# Rotate the galaxy by the same random angle used in GREAT-DES
	e_int= e*np.exp(2*1j*truth['rotation_angle'])
	e1_int= e.real ; e2_int= e.imag

	# Now shear the galaxy
	g1= truth['g1_true'] ; g2= truth['g2_true']
	g= g1+1j*g2

	e_sheared= (e_int+g)/(1+g.conjugate()*e_int) 
	e1_sheared= e_sheared.real ; e2_sheared= e_sheared.imag

	# Unpack the observed ellipticities as required   
	
	if (args.method=='ngmix'):	
		e1_obs= np.transpose(res['g'])[0]  ; e2_obs= np.transpose(res['g'])[1]  
	elif  (args.method=='im3shape'):
		e1_obs= res['e1'] ; e2_obs= res['e2']
	
	# Remove any objects flagged 'do_meas'== -1
	# The reason for this flagging is not apparent for the moment
	# But these objects have e1=e2= -10 in the COSMOS catalogue. 

	# This removes objects without usable ellipticities from the COSMOS catalogue
	# Disable this selection if using matching, as this causes problems taking the intersection later on
	#if not match_res:
	#	e1_int= e1_int[cosmos_cat['do_meas']==0] ; e2_int= e2_int[cosmos_cat['do_meas']==0] 
	#	e1_sheared= e1_sheared[cosmos_cat['do_meas']==0] ; e2_sheared= e2_sheared[cosmos_cat['do_meas']==0]
	#	truth= truth[cosmos_cat['do_meas']==0]
	#	res= res[cosmos_cat['do_meas']==0]
	#	e1_obs= e1_obs[cosmos_cat['do_meas']==0] 

	# Finally append the calculated shapes to the results table

	res= tabletools.appendColumn(res,'e1_sheared', e1_sheared)
	res= tabletools.appendColumn(res,'e2_sheared', e2_sheared) 	
	res= tabletools.appendColumn(res,'e1_intrinsic', e1_int)
	res= tabletools.appendColumn(res,'e2_intrinsic', e2_int)

#	print 'Measured e1=', e1_obs
#	print 'True sheared e1=', e1_sheared
	if len(e1_int[np.isnan(e1_int)])!= 0:
		logger.error('Invalid intrinsic ellipticity found.')
	print 'Added true sheared shape column:', res['e1_sheared']

#	plt.errorbar(e1_obs[0:5000],e1_sheared[0:5000],fmt='.')
#	plt.show()

	return res, truth

#-----------------------------------------------------------------------------------------------------#

def plot_obs_int_shape():
	if (args.method=='ngmix'):	res,truth= load_ngmix_measurements(load_truth=True, load_intrinsic_shapes=True)
	if (args.method=='im3shape'):	res,truth= load_im3shape_measurements(load_truth=True, load_ngmix_columns=False, load_intrinsic_shapes=True)

	res,truth = get_intrinsic_shapes(res,truth)			##NEEDS UPDATING FOR CHANGE IN FUNCTION INPUT##

	print
#-----------------------------------------------------------------------------------------------------#

def get_distributions():
	
	# The data structure and labelling in the im3shape and ngmix results files are slightly different
	# The plotted quantities here are broadly analogous but not identical
	if (args.method=='ngmix'):	res= load_ngmix_measurements(load_truth=False, load_intrinsic_shapes=False)
	if (args.method=='im3shape'):	res= load_im3shape_measurements(load_truth=False, load_ngmix_columns=False, load_intrinsic_shapes=False)

	lab_g1=0 ; lab_shear2=0 ; labflux=0 ; labsize=0 ; labsnr=0   
	if (args.method=='ngmix'): lab_g1= 'g' ; lab_flux= 'log_flux' ; lab_size= 'log_T' ; lab_snr= 's2n_r'
	if (args.method=='im3shape'): lab_g1= 'e1' ; lab_shear2= 'e2' ; lab_flux= 'mean_flux' ; lab_size= 'mean_rgpp_rp' ; lab_snr= 'snr'	

	
	g1=np.array([]) ; g2=np.array([]) ; size= np.array([]) ; flux=np.array([]); snr=np.array([]) ; 
	for i in range(len(res)):
		if (args.method=='ngmix'):	
			g1=np.hstack(( g1, np.transpose(res[i][lab_g1])[0] )) ; g2=np.hstack(( g2, np.transpose(res[i][lab_g1])[1] )) 
		elif  (args.method=='im3shape'):
			g1=np.hstack(( g1, res[i][lab_g1] )) ; g2=np.hstack(( g2, res[i][lab_shear2]))
		size=np.hstack(( size,res[i][lab_size] )) ; 
		flux=np.hstack((res[i][lab_flux] )) ; snr=np.hstack((snr,res[i][lab_snr]))

	# If using ngmix, we convert the size T into linear space
	if (args.method=='ngmix'):
		size= np.exp(size)
		flux= np.exp(flux)

	plt.subplot(2,2,1)
	plt.hist(size, bins= np.linspace(0,2,50), normed=True, histtype='step',color='m')
	plt.xlabel('size')
	plt.axvline(np.mean(size),color='k')
	plt.legend(framealpha=0.0,frameon=False, mode='expand')

	plt.subplot(2,2,2)
	plt.hist(g1, bins= np.linspace(-1.0,1.0,50), normed=True, histtype='step',color='m', label='$g_{1}$')
	plt.hist(g2, bins= np.linspace(-1.0, 1.0,50), normed=True, histtype='step',color='b', label='$g_{2}$')
	plt.xlabel('$g$')
	plt.axvline(np.mean(g1),color='m')
	plt.axvline(np.mean(g2),color='b')	
	plt.legend(framealpha=0.0,frameon=False, mode='expand', loc= 'upper right')

	plt.subplot(2,2,3)
	plt.hist(flux, bins= np.linspace(50,300,50), normed=True, histtype='step',color='m')
	#plt.hist(ngmix_res['flux_true'][0], bins= np.linspace(min(flux),max(flux),50), normed=True, histtype='step',label='true',color='b')
	plt.xlabel('flux')
	#plt.axvline(np.mean(flux[flux>0]),color='m')
	plt.legend(framealpha=0.0,frameon=False, mode='expand')

	plt.subplot(2,2,4)
	plt.hist(snr, bins= np.linspace(0,100,100), normed=True, histtype='step',color='m')
	
	print 100*float(len(snr[snr>1000]))/len(snr)
	#plt.hist(ngmix_res['s2n_true'][0], bins= np.linspace(0,1000,50), normed=True, histtype='step',label='true')
	plt.xlabel('SNR')
	plt.axvline(np.mean(snr[snr>0]),color='m')
	plt.legend(framealpha=0.0,frameon=False, mode='expand')


	print 'snr=',res[0][lab_snr]
	print 'flux=',res[0][lab_flux]

	plt.show()

#-----------------------------------------------------------------------------------------------------#

def get_meane_plots():
	# The parameters measured by the code are not identical
	# The closest analogues are selected for plotting 
	if (args.method=='ngmix'):
		res= load_ngmix_measurements(load_truth=False, load_intrinsic_shapes=False)
		lab_g='g' ; lab_psfg1= 'psf_g' ; lab_size= 'T_s2n_r' ; lab_snr= 's2n_r' ; lab_psf_size= 'psf_T'

	if (args.method=='im3shape'):
		res= load_im3shape_measurements(load_truth=False, load_ngmix_columns=False, load_intrinsic_shapes=False)
		lab_g1='e1' ; lab_g2='e2' ; lab_psfg1= 'mean_psf_e1_sky' ; lab_psfg2= 'mean_psf_e2_sky' ; lab_size= 'mean_rgpp_rp' ; 
		lab_snr= 'snr' ; lab_psf_size= 'mean_psf_fwhm'

	g1=np.array([]); g2=np.array([]); size= np.array([]) ; psf_g1= np.array([]); psf_g2= np.array([]); snr= np.array([]); psf_size= np.array([]) 
	# Cycle through arrays loaded from different files and combine the relevant fields for plotting		
	for i in range( len(res) ):
		if (args.method=='ngmix'):	
			g1=np.hstack(( g1, np.transpose(res[i][lab_g])[0] )) ; g2=np.hstack(( g2, np.transpose(res[i][lab_g])[1] )) 
			psf_g1=np.hstack(( psf_g1, np.transpose(res[i][lab_psfg1])[0] )) ; psf_g2=np.hstack(( psf_g2, np.transpose(res[i][lab_psfg1])[1] ))

		elif  (args.method=='im3shape'):
			g1=np.hstack(( g1, res[i][lab_g1] )) ; g2=np.hstack(( g2, res[i][lab_g2]))
			psf_g1= np.hstack((psf_g1,res[i][lab_psfg1])) ; psf_g2= np.hstack((psf_g2,res[i][lab_psfg2])) 
		size=np.hstack(( size,res[i][lab_size] )) ; 
		snr= np.hstack((snr, res[i][lab_snr])); psf_size= np.hstack((psf_size,res[i][lab_psf_size]))
		
	# If using ngmix convert the size parameter to linear space
	if args.method=='ngmix':
		size=np.exp(size)


	#Define PSF ellipticity bins
	ebins= np.linspace( min(psf_g1), max(psf_g1), 11 )
	x= []; 
	y= []; err= []
	mean_size= [] ; err_size= [] 
	mean_snr= [] ; err_snr= []
	#For each bin, select galaxies and find average g1 and size
	for i in range(len(ebins)-1):
		x+= [ (ebins[i]+ebins[i+1])/2.0 ]
		sel_g1= g1[psf_g1> ebins[i]] ; sel_size= size[psf_g1> ebins[i]] ; sel_snr= snr[psf_g1> ebins[i]]																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																															
		sel_psf= psf_g1[psf_g1> ebins[i]] 
		sel_g1= sel_g1[sel_psf< ebins[i+1]] ; sel_size= sel_size[sel_psf< ebins[i+1]] ; sel_snr= sel_snr[sel_psf< ebins[i+1]]
		sel_psf= sel_psf[sel_psf< ebins[i+1]]
		y+= [ np.mean(sel_g1) ] ; mean_size+= [ np.mean(sel_size) ] ; mean_snr+= [ np.mean(sel_snr) ]
		err+= [ np.std(sel_g1)/len(sel_g1)**0.5 ] ; err_size+= [ np.std(sel_size)/len(sel_size)**0.5 ] ; err_snr+= [ np.std(sel_snr)/len(sel_snr)**0.5 ]

	#np.delete(y,y[np.isnan(y)]) ; np.delete(mean_size,mean_size[np.isnan(mean_size)]) ; np.putmask(mean_snr,np.isnan(mean_snr),0)
	#np.putmask(err,np.isnan(err),0) ; np.putmask(err_size,np.isnan(err_size),0) ; np.putmask(err_snr,np.isnan(err_snr),0)  

	# Output diagnostic plots
	# These quantities should not be correlated 
	x= np.array(x); y= np.array(y); err= np.array(err)
	#font = {'family' : 'normal',
	#        'weight' : 'bold',
	#        'size'   : 18}
	#plt.rc('font', **font)

	plt.subplot(3,3,1)
	plt.axhline(0.0,color='k')
	plt.xlabel('$e_{1}^{psf}$') ; plt.ylabel('$<e^{gal}_{1}>$')
	print 'TEST:',np.isnan(y)
	print err
	#tmp=np.linspace(0,len(x),len(x)+1)
	#tmp= tmp[np.isnan(y)]
	#plt.errorbar(np.delete(x,np.isnan(y)), np.delete(y,np.isnan(y)), np.delete(err,np.isnan(y)), fmt='x')
	#plt.errorbar(np.delete(x,x[np.isnan(y)]), np.delete(y,y[np.isnan(y)]), np.delete(err,err[np.isnan(y)]))

	plt.errorbar(x, y, err, fmt='x')
	plt.errorbar(x, y, err)

	plt.subplot(3,3,2)
	plt.axhline(0.0,color='k')
	plt.xlabel('$e_{1}^{psf}$') ; plt.ylabel(lab_size)
	plt.errorbar(x, mean_size, err_size, fmt='x')
	plt.errorbar(x, mean_size, err_size)

	plt.subplot(3,3,3)
	plt.axhline(0.0,color='k')
	plt.xlabel('$e_{1}^{psf}$') ; plt.ylabel(lab_snr)
	plt.errorbar(x, mean_snr, err_snr, fmt='x')
	plt.errorbar(x, mean_snr, err_snr)

	#Now do the same for SNR
	#Define SNR bins
	snrbins= np.linspace( min(snr), max(snr), 11 )
	x= []; 
	mean_g1= []; err_g1= []; 
	mean_size= [] ; err_size= []
	mean_psfg= [] ; err_psfg= [] 
	#For each bin, select galaxies and find average g1 and size
	for i in range(len(snrbins)-1):
		x+= [ (snrbins[i]+snrbins[i+1])/2.0 ]
		sel_g1= g1[snr> snrbins[i]] ; sel_size= size[snr> snrbins[i]] ; sel_psfg= psf_g1[snr> snrbins[i]]
		sel_snr= snr[snr> snrbins[i]] 
		sel_g1= sel_g1[sel_snr< snrbins[i+1]] ; sel_size= sel_size[sel_snr< snrbins[i+1]] ; sel_psfg= sel_psfg[sel_snr< snrbins[i+1]]
		sel_snr= sel_snr[sel_snr< snrbins[i+1]]
		mean_g1+= [ np.mean(sel_g1) ] ; mean_size+= [ np.mean(sel_size) ] ; mean_psfg+= [ np.mean(sel_psfg) ] 
		err_g1+= [ np.std(sel_g1)/((len(sel_g1))**0.5) ] ; err_size+= [ np.std(sel_size)/((len(sel_size))**0.5) ] ; err_psfg+= [ np.std(sel_psfg)/((len(sel_psfg))**0.5) ] 

		#print 'x=',x
	plt.subplot(3,3,4)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_snr)
	plt.ylabel('$<e^{gal}_{1}>$')
	plt.errorbar(x, mean_g1, err_g1, fmt='x')
	plt.errorbar(x, mean_g1, err_g1)

	plt.subplot(3,3,5)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_snr) ; plt.ylabel(lab_size)
	plt.errorbar(x, mean_size, err_size, fmt='x')
	plt.errorbar(x, mean_size, err_size)

	plt.subplot(3,3,6)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_snr) ; plt.ylabel('$<e^{PSF}_{1}>$')
	plt.errorbar(x, mean_psfg, err_psfg, fmt='x')
	plt.errorbar(x, mean_psfg, err_psfg)

	#And for PSF size
	#Define PSF size bins
	psfTbins= np.linspace( min(psf_size), max(psf_size), 11 )
	x= []; 
	mean_g1= []; err_g1= []; 
	mean_size= [] ; err_size= []
	mean_psfg= [] ; err_psfg= [] 
	mean_snr= [] ; err_snr= [] 
	#For each bin, select galaxies and find average g1 and size
	for i in range(len(psfTbins)-1):
		x+= [ (psfTbins[i]+psfTbins[i+1])/2.0 ]
		sel_g1= g1[psf_size> psfTbins[i]] ; sel_size= size[psf_size> psfTbins[i]] ; sel_psfg= psf_g1[psf_size> psfTbins[i]] ; sel_snr= snr[psf_size> psfTbins[i]]
		sel_psf_size= psf_size[psf_size> psfTbins[i]] 
		sel_g1= sel_g1[sel_psf_size< psfTbins[i+1]] ; sel_size= sel_size[sel_psf_size< psfTbins[i+1]] ; sel_psfg= sel_psfg[sel_psf_size< psfTbins[i+1]] ; sel_snr= sel_snr[sel_psf_size< psfTbins[i+1]]
		sel_psf_size= sel_psf_size[sel_psf_size< psfTbins[i+1]]
		mean_g1+= [ np.mean(sel_g1) ] ; mean_size+= [ np.mean(sel_size) ] ; mean_psfg+= [ np.mean(sel_psfg) ] ; mean_snr+= [ np.mean(sel_snr) ] 
		#print 'snr=',sel_snr
		err_g1+= [ np.std(sel_g1)/(len(sel_g1)**0.5) ] ; err_size+= [ np.std(sel_size)/(len(sel_size)**0.5) ] ; err_psfg+= [ np.std(sel_psfg)/(len(sel_psfg)**0.5) ] ; err_snr+= [ np.std(sel_snr)/(len(sel_snr)**0.5) ] 
		#print 'errsnr=',err_snr

	print 'x=',x
	plt.subplot(3,3,7)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_psf_size) ; plt.ylabel('$<e^{gal}_{1}>$')
	plt.errorbar(x, mean_g1, err_g1, fmt='x')
	plt.errorbar(x, mean_g1, err_g1)

	plt.subplot(3,3,8)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_psf_size) ; plt.ylabel(lab_size)
	plt.errorbar(x, mean_size, err_size, fmt='x')
	plt.errorbar(x, mean_size, err_size)

	plt.subplot(3,3,9)
	plt.axhline(0.0,color='k')
	plt.xlabel(lab_psf_size) ; plt.ylabel(lab_snr)
	#plt.ylim(min(mean_snr),max(mean_snr))
	plt.errorbar(x, mean_snr, err_snr, fmt='x')
	plt.errorbar(x, mean_snr, err_snr)

	plt.show()
			
	print 'Done'

#-----------------------------------------------------------------------------------------------------#

def get_single_shear_bias():
	g1=[] ; g2=[] ; g1_true=[] ; g2_true=[] 
	if args.method== 'im3shape':	
		res, truth_table= load_im3shape_measurements(load_truth=True, load_ngmix_columns=False, load_intrinsic_shapes=False)
		lab_shear1='e1'; lab_shear2='e2'; lab_shear1_true='g1_true' ; lab_shear2_true='g2_true'  

		for i in range( len(res) ):
			g1=np.hstack((g1, res[i][lab_shear1])) ; g2=np.hstack((g2, res[i][lab_shear2])) ;
			g1_true=np.hstack((g1_true, truth_table[i][lab_shear1_true])) ; g2_true=np.hstack((g2_true, truth_table[i][lab_shear2_true])) ; 

		# The g1 shear component is reversed in the im3shape catalogues
		g1*=-1

	if args.method== 'ngmix':	
		res= load_ngmix_measurements(load_truth=False, load_intrinsic_shapes=False)
		lab_shear='g' ; lab_shear_true='shear_true' ;

		# Cycle through arrays loaded from different files and combine the relevant fields 	
		for i in range( len(res) ):
			g1=np.hstack((g1, np.transpose(res[i][lab_shear])[0])) ; g2=np.hstack((g2, np.transpose(res[i][lab_shear])[1])) ;
			g1_true=np.hstack((g1_true, np.transpose(res[i][lab_shear_true])[0])) ; g2_true=np.hstack((g2_true, np.transpose(res[i][lab_shear_true][1])[1])) ; 

	# Find the mean shear of this galaxy sample 
	mean_g1= np.mean(g1) ; mean_g2= np.mean(g2)
	err_g1= np.std(g1)/(len(g1)**0.5) ; err_g2= np.std(g2)/(len(g2)**0.5)

	print 'g1=',g1
	print 'g2=',g2
	print 'g1true=',g1_true[0]
	print 'g2true=',g2_true[0]

	import matplotlib.pyplot as plt
	bins=np.linspace(-1.0,1.0,50)
	plt.hist(g1,bins=bins,histtype='step', normed=True)
	plt.axvline(mean_g1,color='k')
	plt.axvline(g1_true[0],color='m',label='true')
	plt.legend(loc='upper right')
	plt.show()

	plt.hist(g2,bins=bins,histtype='step', normed=True)
	plt.axvline(mean_g2)
	plt.axvline(g2_true[0],color='m')
	plt.show()


	c1=0.0
	y1= [c1]
	y1+= [mean_g1-g1_true[0]]
	erry1= [0.0]
	erry1+= [(err_g1/mean_g1)*(mean_g1-g1_true[0])]
	x1= [0.0]
	x1+= [g1_true[0]]
	if (g1_true[0]<0):
		x1=x1[::-1] ; y1=y1[::-1] ; erry1=erry1[::-1]
	
	#plt.subplot(1,2,1)
	if g1_true[0]!=0:
		m11=(y1[1]-y1[0])/(x1[1]-x1[0])
		errm= abs((err_g1/mean_g1)*m11)
		label= '$g_{1}^{obs} g_{1}^{true}$, $m_{1}=%f \pm %f$' %(m11,errm)
		plt.errorbar(x1,y1,yerr=erry1,color='m',label=label)
		outstr='m1=%f +- %f' %(m11,errm) 
		print (outstr)
	plt.axhline(0,color='k')
	plt.axvline(0,color='k')
	plt.xlabel('$g^{true}$')
	plt.ylabel('$<g^{obs}>-g^{true}$')



	c2=0.0
	y2= [c2]
	y2+= [mean_g2-g2_true[0]]
	erry2= [0.0]
	erry2+= [(err_g2/mean_g2)*(mean_g2-g2_true[0])]
	x2= [0.0]
	x2+= [g2_true[0]]
	if (g2_true[0]<0):
		x2=x2[::-1] ; y2=y2[::-1] ; erry2=erry2[::-1]
	
	if g2_true[0]!=0:
		m22=(y2[1]-y2[0])/(x2[1]-x2[0])
		errm= abs((err_g2/mean_g2)*m22)
		label= '$g_{2}^{obs} g_{2}^{true}$, $m_{2}=%f \pm %f$' %(m22,errm)
		plt.errorbar(x2,y2,yerr=erry2,color='b',label=label)

		outstr='m2=%f +- %f' %(m22,errm) 
		print (outstr)

	plt.legend(loc='lower left')

	# Also plot the cross shear components
	c12=0.0
	y12= [c12]
	y12+= [mean_g1]
	erry12= [0.0]
	erry12+= [err_g1]
	x12= [0.0]
	x12+= [g2_true[0]]
	if (g2_true[0]<0):
		x12=x12[::-1] ; y12=y12[::-1] ; erry12=erry12[::-1]
	
	#plt.subplot(1,2,2)
	if g2_true[0]!=0:
		m12=(y12[1]-y12[0])/(x12[1]-x12[0])
		errm= abs((err_g1/mean_g1)*m12)
		label= '$g^{obs}_{1}$ $g^{true}_{2}$, $m_{12}=%f \pm %f$' %(m12,errm)
		#plt.errorbar(x12,y12,yerr=erry12,color='m',label=label)
		outstr='m12=%f +- %f' %(m12,errm) 
		print (outstr)
	#plt.axhline(0,color='k')
	#plt.axvline(0,color='k')
	#plt.xlabel('$g^{true}$')
	#plt.ylabel('$<g^{obs}>$')

	

	c21=0.0
	y21= [c21]
	y21+= [mean_g2]
	erry21= [0.0]
	erry21+= [err_g2]
	x21= [0.0]
	x21+= [g1_true[0]]
	if (g1_true[0]<0):
		x21=x21[::-1] ; y21=y21[::-1]; erry21=erry21[::-1]
	
	if g1_true[0]!=0:
		m21=(y21[1]-y21[0])/(x21[1]-x21[0])
		errm= abs((err_g2/mean_g2)*m21)
		label= '$g^{obs}_{2}$ $g^{true}_{1}$, $m_{21}=%f \pm %f$' %(m21,errm)
		#plt.errorbar(x21,y21,yerr=erry21,color='b',label=label)
		outstr='m21=%f +- %f' %(m21,errm) 
		print (outstr)
	#plt.legend(loc='lower left')

	
	
	plt.show()

#-----------------------------------------------------------------------------------------------------#

def get_bias():
	if args.method== 'im3shape':	
		res, truth_table= load_im3shape_measurements(load_truth=True, load_ngmix_columns=False, load_intrinsic_shapes=False)

	if args.method== 'ngmix':	
		res= load_ngmix_measurements(load_truth=False, load_intrinsic_shapes=False)
		truth_table= None

	b1,b2,b12,b21= evaluate_bias(res,truth_table,lab='',selection=False)

#-----------------------------------------------------------------------------------------------------#

def evaluate_bias(res, truth_table, lab, selection):
	g1_mean=[] ; g2_mean=[] ; g1_true=[] ; g2_true=[] ; g1_err=[] ; g2_err=[]
	if args.method== 'im3shape':	
		lab_shear1='e1'; lab_shear2='e2'; lab_shear1_true='g1_true' ; lab_shear2_true='g2_true' ; lab_sens1= 'g1_sens_A' ; lab_sens2= 'g2_sens_A'

	if args.method== 'ngmix':	
		lab_shear1='g' ; lab_shear2='g' ; 
		lab_shear_true='shear_true' ; lab_sens1='g_sens' ; lab_sens2='g_sens'

	# If the selection option is turned on use the true sheared shapes
	if selection:	lab_shear1= 'e1_sheared' ; lab_shear2= 'e2_sheared' ;

	# Check whether the files contain a sensitivity correction column
	try:
		test= res[0][lab_sens1]
		use_sens=True
		if not selection:	print 'Using sensitivity corrections.'
		else:			print 'Not using sensitivity corrections.'
	except:
		use_sens=False
		print 'No sensitivities found'

	# Also check for a bias calibration column
	try:
		test= res[0]['nbc_m']
		use_cal=True
		print 'Using bias corrections'
	except:
		use_cal=False
		print 'No bias corrections found'

	# And a weights column
	try:
		test= res[0]['w']
		use_w=True
		if not selection:	print 'Using galaxy weights.'
		else:			print 'Not using galaxy weights.'
	except:
		use_w=False
		print 'No galaxy weights found'

	# If using intrinsic shapes, ignore any sensitivity and weight columns
	# These are not the result of real measurements and so the associated corrections are not needed
	if selection:	use_sens= False ; use_w= False

	# Obtain the first shear group loaded from the config file
	i= config['general']['shear_groups_to_use'][0]
	g1_san_0= config['general']['shear'][i][0]
	# Reverse the e1 shear component if using im3shape results
	#if selection:	import pdb ; pdb.set_trace()
	if args.method== 'im3shape':	g1_san_0*= -1.0
	g2_san_0= config['general']['shear'][i][1]

	g1= np.array([]) ; g2= np.array([]);				# Measured shear estimate
	g1_san= np.array([]) ; g2_san= np.array([])			# True shear
	g1_sen= np.array([]) ; g2_sen= np.array([])			# Shear sensitivities
	g1_nbc= np.array([]) ; g2_nbc= np.array([])  			# Noise bias corrections
	g_w= np.array([])						# Galaxy weights
	n_fil=0								# Number of files loaded
	for single_shear_res in res: 
		# Unpack shear measurements as required for the two codes
		if args.method== 'ngmix':	
			if not selection:	g1_tmp= np.transpose(single_shear_res[lab_shear1])[0] ; g2_tmp=np.transpose(single_shear_res[lab_shear2])[1]
			else:			g1_tmp= single_shear_res[lab_shear1] ; g2_tmp= single_shear_res[lab_shear2]
			g1_san_tmp= np.transpose(single_shear_res[lab_shear_true])[0] ; g2_san_tmp= np.transpose(single_shear_res[lab_shear_true])[1] 
			
		if args.method== 'im3shape':
			## IMPORTANT 								##
			## Note that there is a sign difference between im3shape/ngmix g1 	##	
			g1_tmp= single_shear_res[lab_shear1] ; g2_tmp= single_shear_res[lab_shear2]
			# If using intrinsic shapes, flip e1 to match the im3shape measured e1
			if selection:	g1_tmp*= -1.0
			g1_san_tmp= single_shear_res[lab_shear1_true] ; g2_san_tmp= single_shear_res[lab_shear2_true]

		if len(g1_san_tmp)==0:	continue

		# If the results files contain a sensitivity column, also read this
		# Otherwise set all the sensitivites to 1
		# Do the same with the nbc
		if use_sens:
			if args.method== 'ngmix':	
				g1_sen_tmp= single_shear_res[lab_sens1][:,0] ; g2_sen_tmp= single_shear_res[lab_sens2][:,1]
			if args.method== 'im3shape':
				g1_sen_tmp= single_shear_res[lab_sens1] ; g2_sen_tmp=single_shear_res[lab_sens2]
		if not use_sens:
			g1_sen_tmp= np.ones(np.shape(g1_tmp)) ; g2_sen_tmp= np.ones(np.shape(g2_tmp))

		if use_cal:
			if args.method== 'ngmix':	
				g1_nbc_tmp= np.transpose(single_shear_res['nbc_m'])[0] ; g2_nbc_tmp= np.transpose(single_shear_res['nbc_m'])[1]
			if args.method== 'im3shape':
				g1_nbc_tmp= single_shear_res['nbc_m'] ; g2_nbc_tmp=single_shear_res['nbc_m']
		if not use_cal:
			g1_nbc_tmp= np.zeros(np.shape(g1_tmp)) ; g2_nbc_tmp= np.zeros(np.shape(g2_tmp))

		if use_w:	g_w_tmp= single_shear_res['w']
		if not use_w:	g_w_tmp= np.ones(np.shape(g1_tmp))
		
		# If the true shear group for this datafile matches the previous one, append the data to the relevant arrays but don't evaluate the mean
		if (g1_san_tmp[0]== g1_san_0) and (g2_san_tmp[0]== g2_san_0):
			g1=np.hstack((g1,g1_tmp)) ; g2=np.hstack((g2,g2_tmp))
			g1_san=np.hstack((g1_san,g1_san_tmp)) ; g2_san=np.hstack((g2_san,g2_san_tmp))
			g1_sen=np.hstack((g1_sen,g1_sen_tmp)) ; g2_sen=np.hstack((g2_sen,g2_sen_tmp))
			g1_nbc=np.hstack((g1_nbc,g1_nbc_tmp)) ; g2_nbc=np.hstack((g2_nbc,g2_nbc_tmp))
			g_w=np.hstack((g_w, g_w_tmp))
		# If a new shear group is reached, evaluate the shear estimate
		# Then reinitialise the temporary arrays before writing the new entries to them
		if ((g1_san_tmp[0]!= g1_san_0) or (g2_san_tmp[0]!= g2_san_0)):
			# Find the mean shear of this galaxy sample
			g1_sen=g1_sen[ (np.invert(np.isnan(g1))) & (np.invert(np.isnan(g1_sen))) ] ; g2_sen=g2_sen[(np.invert(np.isnan(g2))) & (np.invert(np.isnan(g2_sen)))]
			g1_nbc=g1_nbc[(np.invert(np.isnan(g1))) & (np.invert(np.isnan(g1_nbc)))] ; g2_nbc=g2_nbc[(np.invert(np.isnan(g2))) & (np.invert(np.isnan(g2)))]
			g_w= g_w[ (np.invert(np.isnan(g1))) & (np.invert(np.isnan(g2))) & (np.invert(np.isnan(g_w))) ]
			g1= g1[np.invert(np.isnan(g1))] ; g2= g2[np.invert(np.isnan(g2))]	

			w= 1.0 ; nbc1= 1.0 ; nbc2= 1.0 ; sen1= 1.0 ; sen2= 1.0
			if use_w:	w= np.copy(g_w)
			if use_cal:	nbc1= 1.0 + g1_nbc ; nbc2= 1.0 + g2_nbc
			if use_sens:	sen1= np.copy(g1_sen) ; sen2= np.copy(g2_sen)
			
			if (use_w) or (use_cal) or (use_sens):	f1= np.sum(nbc1*sen1*w) 		; f2= np.sum(nbc2*sen2*w) 
			else:					f1= np.sum(np.ones(np.shape(g1))) 	; f2= np.sum(np.ones(np.shape(g2)))
						
			#if selection:	import pdb ; pdb.set_trace()
			g1_mean+= [ np.sum(g_w*g1)/f1 ] ; g2_mean+= [ np.sum(g_w*g2)/f2 ]
			g1_err+= [ np.std(g_w*g1)/(len(g1)**0.5) ] ; g2_err+= [ np.std(g_w*g2)/(len(g2)**0.5) ]
			g1_true+= [ g1_san[0] ] ; g2_true+= [ g2_san[0] ] 
			#g1_mean.append(g1) ; g2_mean.append(g2)
			#g1_err.append(np.sqrt(1./g_w)) ; g2_err.append(np.sqrt(1./g_w) )
			#g1_true.append( g1_san ) ; g2_true.append( g2_san ) 

			# Reset arrays for new shear group
			g1= np.array([]) ; g2= np.array([]);
			g1_san= np.array([]) ; g2_san= np.array([])
			g1_sen= np.array([]) ; g2_sen= np.array([])
			g1_nbc= np.array([]) ; g2_nbc= np.array([])
			g_w= np.array([])

			g1= np.hstack((g1,g1_tmp)) ; g2= np.hstack((g2,g2_tmp))
			g1_san= np.hstack((g1_san,g1_san_tmp)) ; g2_san= np.hstack((g2_san,g2_san_tmp))
			g1_sen= np.hstack((g1_sen,g1_sen_tmp)) ; g2_sen= np.hstack((g2_sen,g2_sen_tmp))
			g1_nbc= np.hstack((g1_nbc,g1_nbc_tmp)) ; g2_nbc= np.hstack((g2_nbc,g2_nbc_tmp))
			g_w= np.hstack((g_w,g_w_tmp))
		# If the end of the final shear group is reached, evaluate the mean ellipticity for the selected galaxies
		if (n_fil== len(res)-1):
                        g1_sen=g1_sen[ (np.invert(np.isnan(g1))) & (np.invert(np.isnan(g1_sen))) ] ; g2_sen=g2_sen[(np.invert(np.isnan(g2))) & (np.invert(np.isnan(g2_sen)))]
                        g1_nbc=g1_nbc[(np.invert(np.isnan(g1))) & (np.invert(np.isnan(g1_nbc)))] ; g2_nbc=g2_nbc[(np.invert(np.isnan(g2))) & (np.invert(np.isnan(g2)))]
                        g_w= g_w[ (np.invert(np.isnan(g1))) & (np.invert(np.isnan(g2))) & (np.invert(np.isnan(g_w))) ]
			g1= g1[np.invert(np.isnan(g1))] ; g2= g2[np.invert(np.isnan(g2))]			

			w= 1.0 ; nbc1= 1 ; nbc2= 1.0 ; sen1= 1.0 ; sen2= 1.0
			if use_w:	w= g_w
			if use_cal:	nbc1= 1.0 + g1_nbc ; nbc2= 1.0 + g2_nbc
			if use_sens:	sen1= g1_sen ; sen2= g2_sen
			
			if use_w or use_cal or use_sens:	f1= np.sum(nbc1*sen1*w) ; f2= np.sum(nbc2*sen2*w) 
			else:					f1= np.sum(np.ones(np.shape(g1))) 	; f2= np.sum(np.ones(np.shape(g2)))
			
			g1_mean+= [ np.sum(g_w*g1)/f1 ] ; g2_mean+= [ np.sum(g_w*g2)/f2 ]
			g1_err+= [ np.std(g_w*g1)/(len(g1)**0.5) ] ; g2_err+= [ np.std(g_w*g2)/(len(g2)**0.5) ]
			g1_true+=[ g1_san[0] ] ; g2_true+= [ g2_san[0] ]
			#g1_mean.append(g1) ; g2_mean.append(g2)
			#g1_err.append(np.sqrt(1./g_w)) ; g2_err.append(np.sqrt(1./g_w) )
			#g1_true.append( g1_san ) ; g2_true.append( g2_san ) 
		
			#if selection:	import pdb ; pdb.set_trace()
		# Finally update the true_shear_0 parameter 
		g1_san_0= g1_san_tmp[0]
		g2_san_0= g2_san_tmp[0]
		n_fil+= 1
			

	g1_mean= np.array(g1_mean) ; g2_mean= np.array(g2_mean)
	g1_err= np.array(g1_err) ; g2_err= np.array(g2_err)
	g1_true= np.array(g1_true) ; g2_true= np.array(g2_true)	
	print 'g1=',g1_mean
	print 'g2=',g2_mean
	print 'g1true=',g1_true
	print 'g2true=',g2_true

	# Plot <g1obs>-g1true with g1true and fit to find bias	
	
	x1= g1_true
	y1= g1_mean-g1_true
	erry1= abs(g1_err)

	# Fit a line to this to find m1
	plt.subplot(1,2,1)
	p1,cov1= np.polyfit(x1,y1,1,w=1/erry1,full=False,cov=True)
	m1=p1[0] ; c1=p1[1]
	errm1= np.sqrt(np.abs(cov1[0][0])) ; errc1= np.sqrt(np.abs(cov1[1][1]))
	
	outstr='m1=%f +- %f' %(m1,errm1) 
	print (outstr),
	outstr='c1=%f +- %f' %(c1,errc1) 
	print (outstr)

	label= '$e_{1}^{obs} g_{1}^{true}$, $m_{1}= %f\pm %f$' %(m1,errm1)
	plt.errorbar(x1,y1,yerr=erry1,color='m',fmt='x',label=label)
	plt.plot(x1,x1*m1+c1,color='m')

	# Do the same for g2
	x2= g2_true
	y2= g2_mean-g2_true
	erry2= abs(g2_err)

	p2,cov2= np.polyfit(x2,y2,1,w=1/erry2,cov=True)
	m2=p2[0] ; c2=p2[1]
	errm2= np.sqrt(np.abs(cov2[0][0])) ; errc2= np.sqrt(np.abs(cov2[1][1])) 

	label= '$e_{2}^{obs} g_{2}^{true}$, $m_{2}= %f\pm %f$' %(m2,errm2)
	plt.errorbar(x2,y2,yerr=erry2,color='b',fmt='x',label=label)
	plt.plot(x2,x2*m2+c2,color='b')
	plt.axhline(0,color='k')
	plt.axvline(0,color='k')
	plt.xlabel('$g^{true}$')
	plt.ylabel('$<e^{obs}>-g^{true}$')
	plt.legend(loc='lower right')

	outstr='m2=%f +- %f' %(m2,errm2) 
	print (outstr),
	outstr='c2=%f +- %f' %(c2,errc2) 
	print (outstr)

	# And for the cross terms
	plt.subplot(1,2,2)
	x21= g1_true
	y21= g2_mean
	erry21= abs(g2_err)

	xp21=[] ; yp21=[] ; erryp21=[]
	for i in np.unique(g1_true):
		xp21+= [i]
		yp21+= [ np.sum( g2_mean[g1_true==i] ) ]		
		erryp21+= [ np.sqrt(sum( g2_err[g1_true==i]**2 )) ]
	xp21= np.array(xp21) ; yp21= np.array(yp21) ; erryp21=np.array(erryp21)

	p21,cov21= np.polyfit(xp21,yp21,1,w=1/erryp21,cov=True)
	m21=p21[0] ; c21=p21[1]
	errm21= np.sqrt(np.abs(cov21[0][0])) ; errc21= np.sqrt(np.abs(cov21[1][1])) 

	label= '$e_{2}^{obs} g_{1}^{true}$, $m_{21}= %f\pm %f$' %(m21,errm21)
	plt.errorbar(x21,y21,yerr=erry21,color='m',fmt='x',label=label)
	plt.plot(x21,x21*m21+c21,color='m')

	outstr='m21=%f +- %f' %(m21,errm21) 
	print (outstr),
	outstr='c21=%f +- %f' %(c21,errc21) 
	print (outstr)

	x12= g2_true
	y12= g1_mean
	erry12= abs(g1_err)

	xp12=[] ; yp12=[] ; erryp12=[]
	for i in np.unique(g2_true):
		xp12+= [i]
		yp12+= [ np.sum( g1_mean[g2_true==i] ) ]		
		erryp12+= [ np.sqrt(sum( g1_err[g2_true==i]**2 )) ]
	xp12= np.array(xp12) ; yp12= np.array(yp12) ; erryp12=np.array(erryp12)


	p12,cov12= np.polyfit(xp12,yp12,1,w=1/erryp12,full=False,cov=True)
	m12=p12[0] ; c12=p12[1]
	errm12= np.sqrt(np.abs(cov12[0][0])) ; errc12= np.sqrt(np.abs(cov12[1][1])) 

	label= '$e_{1}^{obs} g_{2}^{true}$, $m_{12}= %f\pm %f$' %(m12,errm12)
	plt.errorbar(x12,y12,yerr=erry12,color='b',fmt='x',label=label)
	plt.plot(x12,x12*m12+c12,color='b')

	outstr='m12=%f +- %f' %(m21,errm21) 
	print (outstr),
	outstr='c12=%f +- %f' %(c21,errc21) 
	print (outstr)

	plt.axhline(0,color='k')
	plt.axvline(0,color='k')
	plt.xlabel('$g^{true}$')
	plt.ylabel('$<e^{obs}>$')
	plt.legend(loc='upper right')

	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/bias_plot_%s.png' %str(lab))
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)

	if args.output=='show_plots':
		plt.show()

	plt.close()
	
	return ((m1,errm1),(c1,errc1)),((m2,errm2),(c2,errc2)),((m12,errm12),(c12,errc12)),((m21,errm21),(c21,errc21))

#-----------------------------------------------------------------------------------------------------#

def evaluate_psf_leakage(res, truth_table, cut_lab):
	if args.method== 'im3shape':	
		lab_shear1='e1'; lab_shear2='e2'; lab_shear1_psf='mean_psf_e1_sky' ; lab_shear2_true='mean_psf_e2_sky'  

	if args.method== 'ngmix':	
		lab_shear='g' ; lab_shear_psf='psf_g' ; lab_sensitivity='g_sens'

	# Look for a sensitivity correction column
	try:
		use_sens=True
	except:
		use_sens=False

	# Unpack the results for all true shear groups into arrays
	g1=np.array([]) ; g2=np.array([]) ; 
	g1_psf= np.array([]) ; g2_psf=np.array([]);
	g1_sens= np.array([]) ; g2_sens=np.array([])  
	for i in range(len(res)):
		if (args.method=='ngmix'):	
			g1=np.hstack(( g1, np.transpose(res[i][lab_shear])[0] ))  
			g2=np.hstack(( g2, np.transpose(res[i][lab_shear])[1] ))
			g1_psf=np.hstack(( g1_psf, np.transpose(res[i][lab_shear_psf])[0] )) 
			g2_psf=np.hstack(( g2_psf, np.transpose(res[i][lab_shear_psf])[1] )) 
			g1_sens=np.hstack(( g1_sens, np.transpose(res[i][lab_sensitivity])[0] )) 
			g2_sens=np.hstack(( g2_sens, np.transpose(res[i][lab_sensitivity])[1] )) 
		elif  (args.method=='im3shape'):
			g1=np.hstack(( g1, res[i][lab_shear1] ))  
			g2=np.hstack(( g2, res[i][lab_shear2]))
			g1_psf=np.hstack(( g1_psf, res[i][lab_shear1_psf] ))
			g2_psf=np.hstack(( g2_psf, res[i][lab_shear2_psf]))
	
	# Plot <e1obs> in bins of e1psf	
	bins= np.linspace(-0.025,0.025,8)
	x1=[] ; y1=[] ; erry1=[]
	y2=[] ; erry2=[]
	
	for i in range(len(bins)-1):	
		# For each bin find the central PSF e1 
		x1+= [ (bins[i]+bins[i+1])/2 ]
		# Make selection and calculate <e1obs> 
		sel_string= 'sel=(g1_psf>%f) & (g1_psf<%f)' %( bins[i], bins[i+1] )
		exec sel_string 
		y1+= [ np.sum(g1[sel])/np.sum(g1_sens[sel]) ]

		tmp= np.std(g1[sel])/np.sqrt(len(g1[sel]))
		tmp= (tmp)/ ( np.mean(g1[sel])/np.sum(g1_sens[sel]) )
		tmp= tmp**2
		tmp+= ( ( np.std(g1_sens[sel])/np.sqrt(len(g1_sens[sel])) ) / np.mean(g1_sens[sel]) )**2 
		tmp= np.sqrt(tmp)*np.sum(g1[sel])/np.sum(g1_sens[sel])		
		erry1+= [tmp] 

		y2+= [ np.sum(g2[sel])/np.sum(g2_sens[sel]) ]

		tmp= np.std(g2[sel])/np.sqrt(len(g2[sel]))
		tmp= (tmp)/ (np.mean(g2[sel]) /np.sum(g2_sens[sel]))
		tmp= tmp**2
		tmp+= ( ( np.std(g2_sens[sel])/np.sqrt(len(g2_sens[sel])) ) / np.mean(g2_sens[sel]) )**2 
		tmp= np.sqrt(tmp)*np.sum(g2[sel])/np.sum(g2_sens[sel])
		erry2+= [tmp] 

	x1= np.array(x1) ; 
	y1= np.array(y1) ; erry1= np.array(erry1)
	y2= np.array(y2) ; erry2= np.array(erry2)

	#print x1
	#print y1
	#print y2
	# Fit lines to these to find alpha+beta
	plt.subplot(1,2,1)
	p1,cov1= np.polyfit(x1,y1,1,w=1/erry1,full=False,cov=True)
	a11=p1[0] ; b11=p1[1]
	erra11= np.sqrt(cov1[0][0]) ; errb11= np.sqrt(cov1[1][1])
	
	p2,cov2= np.polyfit(x1,y2,1,w=1/erry2,full=False,cov=True)
	a21=p2[0] ; b21=p2[1]
	erra21= np.sqrt(cov2[0][0]) ; errb21= np.sqrt(cov2[1][1])	

	outstr='(alpha+beta)= %f +- %f' %(a11,erra11) 
	print (outstr)
	outstr='e1obs(e1psf=0)= %f +- %f' %(b11,errb11) 
	print (outstr)

	outstr='gradient <e2obs> e1psf= %f +- %f' %(a21,erra21) 
	print (outstr)
	outstr='e2obs(e1psf=0)= %f +- %f' %(b21,errb21) 
	print (outstr)

	label= r'$e_{1}^{obs}$, $\alpha+\beta$='
	label+= '%f$\pm$ %f' %(a11,erra11)
	plt.errorbar(x1,y1,yerr=erry1,color='m',fmt='.',label=label)
	plt.plot(x1,x1*a11+b11,color='m')

	label= '$e_{2}^{obs}$ grad=%f$\pm$ %f' %(a21,erra21)
	plt.errorbar(x1,y2,yerr=erry2,color='b',fmt='.',label=label)
	plt.plot(x1,x1*a21+b21,color='b')

	plt.axhline(0,color='k')
	plt.axvline(0,color='k')
	plt.xlabel('$e_{1}^{psf}$')
	plt.ylabel('$<e_{i}^{obs}>$')
	plt.legend(loc='upper left')

	# Do the same for PSF e2
	x2=[] ; 
	y2=[] ; erry2=[]
	y1=[] ; erry1=[]
	
	for i in range(len(bins)-1):	
		x2+= [ (bins[i]+bins[i+1])/2 ]
		sel_string= 'sel=(g2_psf>%f) & (g2_psf<%f)' %( bins[i], bins[i+1] )
		exec sel_string 
		y1+= [ np.sum(g1[sel])/np.sum(g1_sens[sel]) ]

		tmp= np.std(g1[sel])/np.sqrt(len(g1[sel]))
		tmp= (tmp)/ (np.sum(g1[sel])/np.sum(g1_sens[sel]) )
		tmp= tmp**2
		tmp+= ( ( np.std(g1_sens[sel])/np.sqrt(len(g1_sens[sel])) ) / np.mean(g1_sens[sel]) )**2 
		tmp= np.sqrt(tmp)*np.sum(g1[sel])/np.sum(g1_sens[sel])		
		erry1+= [tmp]

		print 'bin', i 
		print 'mg=', np.sum(g2[sel])/np.sum(g2_sens[sel])
		print 'ngal=', len(g2[sel])
		print 'err_mg',np.std(g2[sel])/np.sqrt(len(g2[sel])),  

		y2+= [ np.mean(g2[sel])/np.mean(g2_sens[sel]) ]

		tmp= np.std(g2[sel])/np.sqrt(len(g2[sel]))
		tmp= (tmp)/ (np.sum(g2[sel])/np.sum(g2_sens[sel]) ) 
		tmp= tmp**2
		tmp+= ( ( np.std(g2_sens[sel])/np.sqrt(len(g2_sens[sel])) ) / np.mean(g2_sens[sel]) )**2 
		tmp= np.sqrt(tmp)*np.sum(g2[sel])/np.sum(g2_sens[sel])
		erry2+= [tmp] 
		print tmp

	x2= np.array(x2) ;
	y1= np.array(y1) ; erry1= np.array(erry1) 
	y2= np.array(y2) ; erry2= np.array(erry2)

	# Fit lines to find alpha-beta
	plt.subplot(1,2,2)
	p1,cov1= np.polyfit(x2,y1,1,w=1/erry1,full=False,cov=True)
	a12=p1[0] ; b12=p1[1]
	erra12= np.sqrt(cov1[0][0]) ; errb12= np.sqrt(cov1[1][1])
	
	p2,cov2= np.polyfit(x2,y2,1,w=1/erry2,full=False,cov=True)
	a22=p2[0] ; b22=p2[1]
	erra22= np.sqrt(cov2[0][0]) ; errb22= np.sqrt(cov2[1][1])	

	outstr='gradient <e1obs> e2psf= %f +- %f' %(a12,erra12) 
	print (outstr)
	outstr='e1obs(e2psf=0)= %f +- %f' %(b12,errb12) 
	print (outstr)

	outstr='(alpha-beta)= %f +- %f' %(a22,erra22) 
	print (outstr)
	outstr='e2obs(e2psf=0)= %f +- %f' %(b22,errb22) 
	print (outstr)

	label= '$e_{1}^{obs}$, grad=%f$\pm$ %f' %(a12,erra12)
	plt.errorbar(x1,y1,yerr=erry1,color='m',fmt='.',label=label)
	plt.plot(x2,x2*a12+b12,color='m')

	label= r'$e_{2}^{obs}$, $\alpha-\beta$='
	label+= '%f$\pm$ %f' %(a22,erra22)
	plt.errorbar(x2,y2,yerr=erry2,color='b',fmt='.',label=label)
	plt.plot(x2,x2*a22+b22,color='b')

	plt.axvline(0,color='k')
	plt.axhline(0,color='k')
	plt.xlabel('$e_{2}^{psf}$')
	plt.ylabel('$<e_{i}^{obs}>$')
	plt.legend(loc='lower right')

	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/psf_leakage_plot_%s.png' %str(cut_lab) )
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)

	if args.output=='show_plots':
		plt.show()

	plt.close()

	alpha= 0.5*(a11+a22)
	erralpha= 0.5*np.sqrt(erra11**2 + erra22**2)
	beta= 0.5*(a11-a22)
	errbeta= 0.5*np.sqrt(erra11**2 + erra22**2)
	
	return ((alpha,erralpha),(beta,errbeta))

#-----------------------------------------------------------------------------------------------------#

def plot_bias_cuts():
	
	# If this is true plot the bias alongside the histograms of e1 for the data after the cuts and for the removed objects 
	if args.selection_bias=='1':	show_selection_bias= True
	else:	show_selection_bias= False

	if args.method== 'im3shape':
		if not show_selection_bias:	
			res, truth_table= load_im3shape_measurements(load_truth=False, load_ngmix_columns=False, load_intrinsic_shapes=show_selection_bias)
			selected_truth= None
		else:
			res, truth_table= load_im3shape_measurements(load_truth=False, load_ngmix_columns=False, load_intrinsic_shapes=show_selection_bias)

	if args.method== 'ngmix':
		if not show_selection_bias:	
			res, truth_table= load_ngmix_measurements(load_truth=False, load_intrinsic_shapes=show_selection_bias)
			selected_truth= None
		else:
			res, truth_table= load_ngmix_measurements(load_truth=True, load_intrinsic_shapes=show_selection_bias)

	lab_cut= other_cuts.get_label(quantity= args.quantity[0], mode= args.method)	

	# Calculate the bias repeatedly using different minima 
	cuts= other_cuts.get_bins(quantity= args.quantity[0])
	m1=[] ; errm1=[]; m2=[]; errm2=[]; m12=[] ; errm12=[] ; m21=[] ; errm21=[]
	alpha=[] ; beta=[] ; erralpha=[] ; errbeta=[]
	 
	for j in range(len(cuts)):
		# Extract selection
		selected_res= np.copy(res)
		n1= 0
		n2= 0
		limit1= cuts[j]
		limit2=-1.009804123
		if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':
			# This is just here to avoid an error at the end of the loop
			try:
				limit2= cuts[j+1]
				logger.info('Applying cut %f<%s<%f',cuts[j],lab_cut,cuts[j+1])
			except:	continue
		elif args.actions[0]=='plot_bias_cuts':
				limit2= -1.009804123
				logger.info('Applying cut %s>%f',lab_cut,cuts[j])

		for i in range(len(res)):
			
			sel_string= other_cuts.get_selection_string(i,lab_cut, quantity= args.quantity[0], mode=args.method,limit1=limit1,limit2=limit2)
			exec sel_string
			n1+= len(res[i])
			selected_res[i]= res[i][sel]
			if show_selection_bias:	selected_truth= truth_table[i][sel]
			else:	selected_truth= None
			n2+= len(selected_res[i])
			
		logger.info('Removed %d objects (%f percent)',n1-n2, 100.0*(n1-n2)/float(n1))
		# Calculate the bias for this selection
		lab='%s_%2.3f' %(lab_cut, limit1)
		b1,b2,b12,b21= evaluate_bias(selected_res,selected_truth,lab=lab,selection=False)
		#psf_leakage= evaluate_psf_leakage(selected_res,selected_truth,cut_lab=lab)

		m1+=  [ b1[0][0] ]  ; errm1+=  [ b1[0][1] ]
		m2+=  [ b2[0][0] ]  ; errm2+=  [ b2[0][1] ]
		m12+= [ b12[0][0] ] ; errm12+= [ b12[0][1] ]
		m21+= [ b21[0][0] ] ; errm21+= [ b21[0][1] ]	

		#alpha+= [ psf_leakage[0][0] ] ; erralpha += [ psf_leakage[0][1] ] ; 
		#beta += [ psf_leakage[1][0] ]  ; errbeta += [ psf_leakage[1][1] ] ;   

	m1=np.array(m1) ; errm1=np.array(errm1)
	m2=np.array(m2) ; errm2=np.array(errm2)
	m12=np.array(m12) ; errm12=np.array(errm12)
	m21=np.array(m21) ; errm21=np.array(errm21) 
	#alpha= np.array(alpha) ; erralpha= np.array(erralpha)
	#beta= np.array(beta) ; errbeta= np.array(errbeta) 

	plot_m(cuts,m1,errm1,m2,errm2,m12,errm12,m21,errm21,res,truth_table,lab_cut)

	#plot_psf_leakage(cuts,alpha,erralpha,beta,errbeta,res,lab_cut,truth_table)

#-----------------------------------------------------------------------------------------------------#

def plot_m(cuts,m1,errm1,m2,errm2,m12,errm12,m21,errm21,res,truth,lab):

	m_av= (m1+m2)/2.0
	errm_av= (errm1+errm2)/2

	xmin,xmax= other_cuts.get_xlim(quantity=args.quantity[0])
	fid_cut= other_cuts.get_fiducial_cut(quantity=args.quantity[0]) 
	
	x= cuts	
	if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':	x= (cuts[1:]+cuts[:-1])/2.0

	# Plot the bias as a function of imposed cut
	fig1= plt.figure()
	frame1= fig1.add_axes((0.1,0.3,0.8,0.6))
	plt.errorbar(x,m1,yerr=errm1,color='m',label='$m_{1}$',fmt='.')
	plt.plot(x,m1,color='m')	
	plt.errorbar(x,m2,yerr=errm2,color='b',label='$m_{2}$',fmt='.')
	plt.plot(x,m2,color='b')	
	plt.ylabel('multiplicative bias $m$')
	#plt.ylim(-0.05,0.02)
	plt.xlim(xmin, xmax )
	plt.tick_params(
		axis= 'x',
		which= 'both',
		bottom= 'on',
		top= 'on',
		labelbottom= 'off')

	plt.axhline(0,color='k')
	plt.axhspan(-0.02, 0.02, facecolor='0.5', alpha=0.5)
	if args.actions[0]=='plot_bias_cuts':	plt.xlabel('Lower %s cut' %(lab))
	if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':	plt.xlabel('%s bin centre' %(lab))
	plt.legend(loc='lower right')
	plt.axvline(fid_cut,linestyle='--',color='k')

	# Plot the distribution of the quantity cut on for the full dataset for comparison
	frame2=fig1.add_axes((0.1,0.1,0.8,0.2))

	q=[]
	#import pdb
	#pdb.set_trace()
	for i in range(len(res)):
		if i==0:
			if args.quantity[0]== 'g1_sens':	q= np.transpose(res[i][lab])[0]
			elif args.quantity[0]=='g2_sens':	q= np.transpose(res[i][lab])[1]
			elif args.quantity[0]=='T_r':	q= np.exp(res[i][lab])
			else:	q= res[i][lab]
		if i!=0:
			if args.quantity[0]== 'g1_sens':	q=np.hstack((q,np.transpose(res[i][lab])[0]))	
			elif args.quantity[0]=='g2_sens':	q=np.hstack((q,np.transpose(res[i][lab])[1]))
			elif args.quantity[0]=='T_r':	q=np.hstack((q,np.exp(res[i][lab])))
			else:	q=np.hstack((q,res[i][lab]))
	q=np.array(q)

	bins=np.linspace(xmin,xmax,61)
	plt.hist(q,bins=bins,histtype='step',normed=True,color='k')
	plt.axvline(fid_cut,linestyle='--',color='k')
	if args.actions[0]=='plot_bias_cuts':	plt.xlabel('Lower %s cut' %(lab))
	if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':	plt.xlabel('%s bin centre' %(lab))
	plt.xlim(xmin,xmax)
	
	fig2= plt.figure()
	frame1= fig2.add_axes((0.1,0.3,0.8,0.6))
	plt.errorbar(x,m12,yerr=errm12,color='m',label='$m_{12}$',fmt='.')
	plt.plot(x,m12,color='m')	
	plt.errorbar(x,m21,yerr=errm21,color='b',label='$m_{21}$',fmt='.')
	plt.plot(x,m21,color='b')	
	plt.ylabel('multiplicative bias $m$')
	plt.xlim(xmin,xmax)

	plt.tick_params(
		axis= 'x',
		which= 'both',
		bottom= 'on',
		top= 'on',
		labelbottom= 'off')

	plt.axhline(0,color='k')
	plt.axvline(fid_cut,linestyle='--',color='k')	
	plt.legend(loc='lower right')

	# Plot the distribution of the quantity cut on for the full dataset for comparison
	frame2=fig2.add_axes((0.1,0.1,0.8,0.2))
	bins=np.linspace(xmin,xmax,61)
	plt.hist(q,bins=bins,histtype='step',normed=True,color='k')
	plt.axvline(fid_cut,linestyle='--',color='k')
	if args.actions[0]=='plot_bias_cuts':	plt.xlabel('Lower %s cut' %(lab))
	if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':	plt.xlabel('%s bin centre' %(lab))
	plt.xlim(xmin,xmax)

	fig3= plt.figure()
	frame1= fig3.add_axes((0.1,0.3,0.8,0.6))
	plt.errorbar(x,m_av,yerr=errm_av,color='g',fmt='.')
	plt.plot(x,m_av,color='g')	
	plt.ylabel('multiplicative bias $|m|$')
	plt.xlim(xmin,xmax)

	plt.tick_params(
		axis= 'x',
		which= 'both',
		bottom= 'on',
		top= 'on',
		labelbottom= 'off')

	plt.axhline(0,color='k')
	plt.axvline(fid_cut,linestyle='--',color='k')	
	plt.axhspan(-0.02, 0.02, facecolor='0.5', alpha=0.5) 
		
	# If required plot the distribution of the quantity cut on for the full dataset for comparison
	frame2=fig3.add_axes((0.1,0.1,0.8,0.2))
	if args.selection_bias== 0:
		bins=np.linspace(xmin,xmax,61)
		plt.hist(q,bins=bins,histtype='step',normed=True,color='k')
		plt.axvline(fid_cut,linestyle='--',color='k')
		plt.tick_params(
			axis= 'y',
			which= 'both',
			bottom= 'on',
			top= 'on',
			labelleft= 'off')
	# Otherwise plot the mean intrinsic shape as a funtion of the relevant cut
	elif args.selection_bias==1:
		logger.info('Plotting intrinsic ellipticities.')
		cuts=other_cuts.get_bins(args.quantity[0])
		excess_e1_int=[] ; err_excess_e1_int=[]
		
                truth_all=np.array([]) ; res_all=np.array([])
                for i in range(len(truth)):
			if i==0:	truth_all= truth[i] ; res_all= res[i]
                	else:	truth_all=  np.hstack(( truth_all, truth[i] )) ; res_all=  np.hstack(( res_all, res[i] ))
		
		mean_int_e1_all= np.mean(res_all['e1_intrinsic'])
		logger.info('Found mean intrinsic e1=%f +- %f', mean_int_e1_all, np.std(res_all['e1_intrinsic'])/np.sqrt(len(res_all)))
 	
		for threshold in cuts:
			print 'Cutting %s>%f' %(lab, threshold)
			sel_int_e1= res_all[ res_all[lab]>threshold ]['e1_intrinsic']

			excess_e1_int+= [ np.mean(sel_int_e1) ]
			err1= np.std(sel_int_e1)/np.sqrt(len(sel_int_e1))
			err2= np.std(res_all['e1_intrinsic'])/np.sqrt(len(res_all))
			err2=0
			err_excess_e1_int+= [ np.sqrt(err1**2 + err2**2) ]
		excess_e1_int= np.array(excess_e1_int) ; err_excess_e1_int= np.array(err_excess_e1_int)			

		plt.errorbar(x,excess_e1_int, yerr=err_excess_e1_int, fmt='.')
		plt.ylabel('$e^{int}_{1}$')
                plt.axhline(mean_int_e1_all,color='k')

	if args.actions[0]=='plot_bias_cuts':	plt.xlabel('Lower %s cut' %(lab))
	if args.actions[0]=='plot_bias_binned' or args.actions[0]=='get_polyfit':	plt.xlabel('%s bin centre' %(lab))
	plt.xlim(xmin,xmax)

	plt.show()
	# Plot the intrinsic shape histogram
	if args.selection_bias==1:
		med= np.median(q)

		fig4=plt.figure()	
		frame1= fig4.add_axes((0.1,0.1,0.8,0.6))
		bins=np.linspace(-1,1,50)
		plt.hist( res_all['e1_intrinsic'], bins=bins, histtype='step',normed=True, color='k', label='Full dataset' )
		plt.axvline(np.mean(res_all['e1_intrinsic']),color='k')
		plt.hist(res_all['e1_intrinsic'][res_all[lab]>med], bins=bins, histtype='step',normed=True, color='g', label='%s>%2.2f' %(lab,med) )
		plt.axvline(np.mean(res_all['e1_intrinsic'][res_all[lab]>med]), color='g')
		plt.legend(loc='upper right')
		plt.xlabel('$e^{int}_{1}$')		
	
	plt.show()
	plt.close()

#-----------------------------------------------------------------------------------------------------#

def plot_psf_leakage(cuts,alpha,erralpha,beta,errbeta,res,lab_cut):

	xmin,xmax= other_cuts.get_xlim(quantity= args.quantity[0])
	fid_cut= get_fiducial_cut(quantity= args.quantity[0])

	# Plot the bias as a function of lower cut
	fig= plt.figure()
	frame1= fig.add_axes((0.1,0.3,0.8,0.6))
	plt.errorbar(cuts,alpha,yerr=erralpha,color='m',label=r'$\alpha$',fmt='.')
	plt.plot(cuts,alpha,color='m')	
	plt.errorbar(cuts,beta,yerr=errbeta,color='b',label=r'$\beta$',fmt='.')
	plt.plot(cuts,beta,color='b')
	plt.axvline(10,linestyle='--',color='k')	
	plt.ylabel(r'PSF Leakage')
	plt.xlim(xmin,xmax)

	plt.tick_params(
		axis= 'x',
		which= 'both',
		bottom= 'on',
		top= 'on',
		labelbottom= 'off')

	plt.axhline(0,color='k')	
	plt.legend(loc='lower right')

	# Plot the histogram of the full dataset for comparison
	frame2=fig.add_axes((0.1,0.1,0.8,0.2))

	snr=[]
	for i in range(len(res)):
		snr=np.hstack((snr,res[i][lab_cut]))
	snr=np.array(snr)

	bins=np.linspace(xmin,xmax,61)
	plt.hist(snr,bins=bins,histtype='step',normed=True,color='k')
	plt.axvline(fid_cut,linestyle='--',color='k')
	plt.xlabel('Lower %s cut' %lab_cut)
	plt.xlim(xmin,xmax)

	plt.show()
	plt.close()

#-----------------------------------------------------------------------------------------------------#

def plot_bias_binned():
	logger.info('Calculating bias in bins of %s', args.quantity[0])
	plot_bias_cuts()

#-----------------------------------------------------------------------------------------------------#

def plot_bias_zbinned():
	"""
	Load the specified bias table with the files. Interpolate a bias correction for each galaxy. 
	Plot this bias in z bins. 
	Then recalculate the residual bias after correction, average over z bins and add the resulting points to the plot. 
	"""

	logger.info('Calculating bias in bins of redshift')

	global zbinning
	zbinning= True

	# Load the results files and add calibration columns
	if args.method=='im3shape':	res,truth= load_im3shape_measurements(load_truth=True, load_ngmix_columns=False, load_intrinsic_shapes=False)
	if args.method=='ngmix':	res,truth= load_ngmix_measurements(load_truth=True, load_intrinsic_shapes=False)

	zbins= np.array(config['general']['z_bins'])
	zbin_centres= (zbins[1:]+zbins[:-1])/2.0

	res_all=np.array([])
	for single_shear_res in res: 
		# Unpack shear measurements as required
		if len(res_all)==0: 	res_all= np.copy(single_shear_res) 
		else:			res_all= np.hstack((res_all, single_shear_res))			

	# First plot the mean bias as calculated using the interpolation table in bins of redshift
	logger.info('Calculating bias in z bins %f', zbin_centres) 
	m_zbinned=[] ; m_z_err=[]
	for i in range(len(zbins)-1):
		zmin= zbins[i] ; zmax= zbins[i+1]   
		sel= (res_all['zphot']>zmin) & (res_all['zphot']<zmax)

		selected_res= res_all[sel] 	
		m_zbinned+= [np.mean(selected_res['nbc_m'])]
		m_z_err+= [np.std(selected_res['nbc_m'])/np.sqrt(len(selected_res))]
	m_zbinned= np.array(m_zbinned) ; m_z_err= np.array(m_z_err) 

	plt.errorbar(zbin_centres, m_zbinned, yerr=m_z_err, fmt='x', label='Uncorrected.')
	print m_zbinned
	print m_z_err

	plt.xlabel('photo-z bin')
	plt.ylabel('multiplicative bias m')
	plt.legend(loc='upper right')
	plt.xticks(zbins)
	plt.xlim(zbins[0],zbins[-1])
	plt.grid()
	#plt.axhline(0,color='k')
	plt.axhspan(-0.02, 0.02, facecolor='0.5', alpha=0.5)

	plt.show()
	exit()

	# Now recalculate the residual bias after calibration
	res_all= get_residual_bias(res) 

	# Finally plot the bias in z bins as previously
	logger.info('Recalculating bias in z bins %f', zbin_centres) 
	m_zbinned=[] ; m_z_err=[]
	for i in range(len(zbins)-1):
		zmin= zbins[i] ; zmax= zbins[i+1]   
		sel= (res_all['zphot']>zmin) & (res_all['zphot']<zmax)

		selected_res= res_all[sel] 	
		m_zbinned+= [np.mean(selected_res['residual_nbc_m'])]
		m_z_err+= [np.std(selected_res['residual_nbc_m'])/np.sqrt(len(selected_res))]
	m_zbinned= np.array(m_zbinned) ; m_z_err= np.array(m_z_err) 

	plt.errorbar(zbin_centres, m_zbinned, yerr=m_z_err, fmt='x', label='Corrected.')
	print m_zbinned
	print m_z_err

	plt.xlabel('photo-z bin')
	plt.ylabel('multiplicative bias m')
	plt.legend(loc='upper right')
	plt.xticks(zbins)
	plt.xlim(zbins[0],zbins[-1])
	plt.grid()
	plt.axhline(0,color='k')
	plt.axhspan(-0.02, 0.02, facecolor='0.5', alpha=0.5)
	plt.show()

#-----------------------------------------------------------------------------------------------------#

def get_residual_bias(res):

	# Set up interpolator as above, but using residual bias tables
	interpolator= initialise_calibration(mode='calibrated')

	# Apply interpolator to calculate new bias corections
	res_all= np.array([])
	for res_file in res:
		m= get_calibration(interpolator, res_file)
		res_file= tabletools.appendColumn(res_file,'residual_nbc_m',m)
		if len(res_all)==0:	res_all= np.copy(res_file)
		else:			res_all= np.hstack((res_all,res_file))
	
	return res_all

#-----------------------------------------------------------------------------------------------------#

def export_ngmix_columns():
	## It may be useful to comment out the call to select() in the load function before running this function ##
	## Otherwise, the exported columns will include flag and fiducial cuts               			  ##
	
	if args.method!= 'ngmix':
		logger.error('Cannot run this command using im3shape results.')
		exit()

	res,truth= load_ngmix_measurements(load_truth= False, load_intrinsic_shapes=False)

	lab_snrT= 'T_s2n_r' ; lab_snr= 's2n_r' ; lab_flags='flags'; lab_sens='g_sens'

        snr_T= np.array([]) ; snr=np.array([]) ; flags=np.array([]) ; g1_sens=np.array([]) ; g2_sens=np.array([])
        for i in range(len(res)):
                snr= np.hstack( (snr, res[i][lab_snr] ) ) ; snr_T= np.hstack((snr_T, res[i][lab_snrT]))
		flags= np.hstack( (flags, res[i][lab_flags] ) ) ;
		g1_sens= np.hstack(( g1_sens, np.transpose(res[i][lab_sens])[0] ))
		g2_sens= np.hstack(( g2_sens, np.transpose(res[i][lab_sens])[1] ))


	l= len(snr)/10000
	snr= np.split(snr,l) ; snr_T= np.split(snr_T,l) ; flags= np.split(flags,l) ; g1_sens= np.split(g1_sens,l)  ; g2_sens= np.split(g2_sens,l) 

	ig= config['general']['shear_groups_to_use'][0]

	outfile1= config['ngmix']['pickled_columns'] %('snr_r', ig)
	outfile2= config['ngmix']['pickled_columns'] %('snr_T', ig)
	outfile3= config['ngmix']['pickled_columns'] %('flags', ig)
	outfile4= config['ngmix']['pickled_columns'] %('g1_sens', ig)
	outfile5= config['ngmix']['pickled_columns'] %('g2_sens', ig)

	pickle.dump( snr, open(outfile1,'wb') )
	logger.info('Wrote columns to %s',outfile1)
	pickle.dump( snr_T, open(outfile2,'wb') )
	logger.info('Wrote columns to %s',outfile2)
	pickle.dump( flags, open(outfile3,'wb') )
	logger.info('Wrote columns to %s',outfile3)
	pickle.dump( g1_sens, open(outfile4,'wb') )
	logger.info('Wrote columns to %s',outfile4)
	pickle.dump( g2_sens, open(outfile5,'wb') )
	logger.info('Wrote columns to %s',outfile5)


#-----------------------------------------------------------------------------------------------------#

def plot_matched_measurements():
	res,truth= load_im3shape_measurements(load_truth=False,load_ngmix_columns=True)

	matched_plot_routines.plot_scatter(res)

#-----------------------------------------------------------------------------------------------------#

def get_polyfit():

	load_truth= False
	if zbinning:	load_truth==True

	# Load files	
	if args.method== 'im3shape':
		res, truth_table= load_im3shape_measurements(load_truth=load_truth, load_ngmix_columns=False, load_intrinsic_shapes=False)
		selected_truth= None

	if args.method== 'ngmix':
		res, truth_table= load_ngmix_measurements(load_truth=load_truth, load_intrinsic_shapes=False)
		selected_truth= None

	# Define the snr and size bins to use
	lab_cut= other_cuts.get_label(quantity= args.quantity[0], mode= args.method)	
	size_bins=np.array([1.15, 1.2,1.23,1.25,1.3,1.4,1.5,1.75])

	res_all=np.copy(res)
	m_grid=[] ; errm=[] ; x=[]
	
	for si in range(len(size_bins)-1):
		res= np.copy(res_all)
		logger.info('Using rgpp_rp=%2.3f-%2.3f', size_bins[si], size_bins[si+1])

		# Calculate the bias repeatedly in different bins 
		cuts= other_cuts.get_bins(quantity= args.quantity[0])
		m1=[] ; errm1=[]; m2=[]; errm2=[]; m12=[] ; errm12=[] ; m21=[] ; errm21=[]
		indices=[]
		alpha=[] ; beta=[] ; erralpha=[] ; errbeta=[]
		 
		for j in range(len(cuts)):
			# Extract selection
			selected_res= np.copy(res)
			n1= 0
			n2= 0
			limit1= cuts[j]
			limit2=-1.009804123
			# This is just here to avoid an error at the end of the loop
			try:
				limit2= cuts[j+1]
				logger.info('Applying cut %f<%s<%f',cuts[j],lab_cut,cuts[j+1])
			except:	continue
			
			for i in range(len(res)):
			
				sel_string= other_cuts.get_selection_string(i,lab_cut, quantity= args.quantity[0], mode=args.method,limit1=limit1,limit2=limit2)
				sel_string+= "& (res[i]['mean_rgpp_rp']>size_bins[si]) & (res[i]['mean_rgpp_rp']<size_bins[si+1])"
				exec sel_string
				n1+= len(res[i])
				selected_res[i]= res[i][sel]
				selected_truth= None

				n2+= len(selected_res[i])
			
			logger.info('Removed %d objects (%f percent)',n1-n2, 100.0*(n1-n2)/float(n1))
			
			# Calculate the bias for this selection
			lab='%s_%2.3f' %(lab_cut, limit1)
			b1,b2,b12,b21= evaluate_bias(selected_res,selected_truth,lab=lab,selection=False)
			#except:	import pdb; pdb.set_trace() ; continue
			m1+=  [ b1[0][0] ]  ; errm1+=  [ b1[0][1] ]
			m2+=  [ b2[0][0] ]  ; errm2+=  [ b2[0][1] ]
			m12+= [ b12[0][0] ] ; errm12+= [ b12[0][1] ]
			m21+= [ b21[0][0] ] ; errm21+= [ b21[0][1] ]	
			indices+=[j]

		m1=np.array(m1) ; errm1=np.array(errm1)
		m2=np.array(m2) ; errm2=np.array(errm2)
		m12=np.array(m12) ; errm12=np.array(errm12)
		m21=np.array(m21) ; errm21=np.array(errm21) 
		indices=np.array(indices)

		if args.output[0]=='save_plots' or args.output[0]=='show_plots':
			plot_m(cuts,m1,errm1,m2,errm2,m12,errm12,m21,errm21,res,truth_table,lab_cut)
		else:
			x+= [(cuts[1:]+cuts[:-1])[indices]/2.0]			# snr bin centres
			m_grid+= [(m1+m2)/2]					# grid of m in snr-size space
			errm+= [(errm1+errm2)/2]		
			
	m_grid= np.array(m_grid)
	errm= np.array(errm)
	x=np.array(x)

	# Fit a polynomial function m=f(snr) in each size bin
	from scipy.optimize import curve_fit
	colourscale= tktools.plotstools.get_colorscale(len(size_bins),cmap_name='nipy_spectral')
	f_grid= []
	fig1= plt.figure()
	for i in range(len(m_grid)):
		plt.errorbar(x[i],m_grid[i], yerr=errm[i], fmt='.', color=colourscale[i], label='$R_{gpp}/R_{p}$=%2.2f-%2.2f' %(size_bins[i], size_bins[i+1]) )
		f_snr_par,f_cov= curve_fit(fitting_routines.snr_basis, x[i], m_grid[i])
		x_high_res= np.linspace(8,200,500)
		f_snr= fitting_routines.snr_basis(x_high_res, f_snr_par[0], f_snr_par[1], f_snr_par[2])
		f_grid+= [f_snr]						# Save the polynomial fits to a grid in 2d snr-size space
		plt.plot(x_high_res,f_snr, color=colourscale[i])
	
	f_grid= np.array(f_grid)

	size_bin_centres= (size_bins[1:]+size_bins[:-1])/2.0
	snr_bin_centres= np.copy(x[0])
	data_grid= np.copy(m_grid)
	noise_grid= np.copy(errm)

	# Now interpolate between them to produce a high resolution bias table
	y_high_res= np.linspace(size_bins[0],size_bins[-1],500)
	interpolator= interpolate.interp2d(x_high_res, size_bin_centres , f_grid, kind='cubic')
	bias_table= interpolator(x_high_res, y_high_res)

	# Try fitting the analytic function to the data
	logger.info('Fitting analytic bias function.')
	alpha, beta, gamma, delta= fitting_routines.get_fit_parameters(snr_bin_centres, size_bin_centres, data_grid, noise_grid)

	# Also use the analytic form from CFHTLenS
	bias_table_analytic=[]
	for size in y_high_res:
		bias_table_analytic+= [ fitting_routines.cfhtlens_bias(x_high_res, size, alpha=alpha, beta=beta, gamma=gamma, delta=delta) ]
	bias_table_analytic= np.array(bias_table_analytic)

	# For comparison also try the fiducial values of alpha,beta given in Miller et al (2013) 
	bias_table_analytic_fid=[]
	for size in y_high_res:
		bias_table_analytic_fid+= [ fitting_routines.cfhtlens_bias(x_high_res, size, alpha=0.057, beta=-0.37, gamma=0.0, delta=1.0) ]
	bias_table_analytic_fid= np.array(bias_table_analytic_fid)

	bias_out_dir= 'bias_tables_dir'	

	# Save the polynomials and the smoothed bias table as pickles
	if args.calibration=='none':
		outfile1= 'bias_data_%d.p' %args.number
		outfile1=os.path.join(config[args.method][bias_out_dir],outfile1)
		pickle.dump( [snr_bin_centres.astype( np.dtype('>f8') ), size_bin_centres.astype( np.dtype('>f8') ), m_grid, noise_grid], open(outfile1,'wb'))
		logger.info('Wrote binned bias data to %s',outfile1)	

		outfile1= 'bias_polynomials_%d.p'%args.number
		outfile1= os.path.join(config[args.method][bias_out_dir],outfile1)
		pickle.dump( [x_high_res.astype( np.dtype('>f8') ), size_bin_centres.astype( np.dtype('>f8') ), f_grid], open(outfile1,'wb'))
		logger.info('Wrote bias polynomial grid to %s',outfile1)

		outfile2='Kacprzak_1d_bias_table_%d.p'%args.number
		outfile2= os.path.join(config[args.method][bias_out_dir],outfile2)
		pickle.dump( [x_high_res.astype( np.dtype('>f8') ), y_high_res.astype( np.dtype('>f8') ), bias_table], open(outfile2,'wb'))
		logger.info('Wrote interpolated bias table to %s',outfile2)

		outfile3='analytic_bias_table_%d.p'%args.number
		outfile3= os.path.join(config[args.method][bias_out_dir],outfile3)
		pickle.dump( [x_high_res.astype( np.dtype('>f8') ), y_high_res.astype( np.dtype('>f8') ), bias_table_analytic], open(outfile3,'wb'))
		logger.info('Wrote analytic bias table to %s',outfile3)

		outfile4='analytic_bias_table_fid_%d.p'%args.number
		outfile4= os.path.join(config[args.method][bias_out_dir],outfile4)
		pickle.dump( [x_high_res.astype( np.dtype('>f8') ), y_high_res.astype( np.dtype('>f8') ), bias_table_analytic_fid], open(outfile4,'wb'))
		logger.info('Wrote analytic bias table to %s',outfile4)
	else:
		bias_out_dir= 'residual_' + bias_out_dir
		outfile= '%s_bias_table_fid_%d.p'%(args.number, args.number)
		outfile= os.path.join(config[args.method][bias_out_dir],outfile)
		if args.calibration== 'Kacprzak_1d':	bias_table_out= np.copy(bias_table)
		if args.calibration== 'Miller':		bias_table_out= np.copy(bias_table_analytic) 
		if args.calibration== 'Miller_fid':	bias_table_out= np.copy(bias_table_analytic_fid)
		pickle.dump( [x_high_res.astype( np.dtype('>f8') ), y_high_res.astype( np.dtype('>f8') ), bias_table_analytic_fid], open(outfile,'wb'))
		logger.info('Wrote analytic bias table to %s',outfile)


								## Various comparison plots ##

	# Plot the fitted polynomials
	plt.ylim(-0.4,0.08)
	plt.xlim(5,90)
	plt.xticks(cuts[:-1])
	plt.grid()
	plt.axhline(0,color='k')
	plt.xlabel('%s bin edge' %args.quantity[0])
	plt.ylabel('Multiplicative bias |m|')
	plt.title('Kacprzak interpolated polynomials')
	plt.legend(loc='lower right')
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/m-snr_1d_fit_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()
	#plt.show()

	# Do the equivalent plot for the fitted Miller function
	fig2= plt.figure()
	for i in range(len(m_grid)):
		plt.errorbar(x[i],m_grid[i], yerr=errm[i], fmt='.', color=colourscale[i], label='$R_{gpp}/R_{p}$=%2.2f-%2.2f' %(size_bins[i], size_bins[i+1]) )
		x_high_res= np.linspace(8,200,500)
		binned_miller= fitting_routines.cfhtlens_bias(x_high_res, size_bin_centres[i], alpha=alpha, beta=beta, gamma=gamma, delta=delta)
		plt.plot(x_high_res,binned_miller, color=colourscale[i])
	
	plt.ylim(-0.4,0.08)
	plt.xlim(5,90)
	plt.xticks(cuts[:-1])
	plt.grid()
	plt.axhline(0,color='k')
	plt.xlabel('%s bin edge' %args.quantity[0])
	plt.ylabel('Multiplicative bias |m|')
	plt.title('Miller function')
	plt.legend(loc='lower right')
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/m-snr_miller_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()
	#plt.show()

	# Do the equivalent plot for the fiducial Miller function
	fig2= plt.figure()
	for i in range(len(m_grid)):
		plt.errorbar(x[i],m_grid[i], yerr=errm[i], fmt='.', color=colourscale[i], label='$R_{gpp}/R_{p}$=%2.2f-%2.2f' %(size_bins[i], size_bins[i+1]) )
		x_high_res= np.linspace(8,200,500)
		binned_miller= fitting_routines.cfhtlens_bias(x_high_res, size_bin_centres[i], alpha=0.057, beta=-0.37, gamma=0, delta=1.0)
		plt.plot(x_high_res,binned_miller, color=colourscale[i])
	
	plt.ylim(-0.4,0.08)
	plt.xlim(5,90)
	plt.xticks(cuts[:-1])
	plt.grid()
	plt.axhline(0,color='k')
	plt.xlabel('%s bin edge' %args.quantity[0])
	plt.ylabel('Multiplicative bias |m|')
	plt.title('Miller function')
	plt.legend(loc='lower right')
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/m-snr_miller_fid_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()

	# Plot the bias surfaces for comparison
	fig, ax1= plt.subplots()
	plt.imshow(bias_table, origin='lower', interpolation='none')
	plt.title('Kacprzak interpolated polynomial')
	plt.xlabel('SNR')
	plt.ylabel('$R_{gp}/R_{p}$')
	plt.xticks(np.arange(1,len(bias_table[0]),len(bias_table[0])/11))
	plt.yticks(np.arange(1,len(bias_table),len(bias_table)/7))
	ax1.set_xticklabels(cuts)
	ax1.set_yticklabels(size_bins)
	plt.colorbar()
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/bias_surface_Kacprzak_1d_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()
	#plt.show()

	fig, ax2= plt.subplots()
	plt.imshow(bias_table_analytic, origin='lower', interpolation='none')
	plt.title('Miller bias function (fitted)')
	plt.xlabel('SNR')
	plt.ylabel('$R_{gp}/R_{p}$')
	plt.xticks(np.arange(1,len(bias_table_analytic[0]),len(bias_table_analytic[0])/11))
	plt.yticks(np.arange(1,len(bias_table_analytic),len(bias_table_analytic)/7))
	ax2.set_xticklabels(cuts)
	ax2.set_yticklabels(size_bins)
	plt.colorbar()
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/bias_surface_miller_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()
	#plt.show()
	
	fig, ax3= plt.subplots()
	plt.imshow(bias_table_analytic_fid, origin='lower', interpolation='none')
	plt.title('Miller bias function (unfitted)')
	plt.xlabel('SNR')
	plt.ylabel('$R_{gp}/R_{p}$')
	plt.xticks(np.arange(1,len(bias_table_analytic_fid[0]),len(bias_table_analytic_fid[0])/11))
	plt.yticks(np.arange(1,len(bias_table_analytic_fid),len(bias_table_analytic_fid)/7))
	ax3.set_xticklabels(cuts)
	ax3.set_yticklabels(size_bins)
	plt.colorbar()
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/bias_surface_miller_fid_%d.png'%args.number)
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
		plt.close()
	if args.output=='show_plots':	plt.show()

#-----------------------------------------------------------------------------------------------------#

def initialise_calibration(mode):
	
	# Load bias table and set up interpolator
	import cPickle as pickle
	bias_table_string=''
	alpha_table_string=''
	if mode=='uncalibrated':
		if args.calibration== 'Kacprzak_1d':	bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'Kacprzak_1d_bias_table_200.p')
		if args.calibration== 'Kacprzak_2d':	
							bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'Kacprzak_2d_bias_table.p')
							alpha_table_string= os.path.join(config[args.method]['bias_tables_dir'],'Kacprzak_2d_alpha_table.p')
		if args.calibration== 'Miller':		bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'analytic_bias_table_200.p')
		if args.calibration== 'Miller_fid':	bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'analytic_bias_table_fid_200.p')

	elif mode=='calibrated':
		if args.calibration== 'Kacprzak_1d':	bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'residual/Kacprzak_1d_bias_table.p')
		if args.calibration== 'Kacprzak_2d':	
							bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'residual/Kacprzak_2d_bias_table.p')
							alpha_table_string= os.path.join(config[args.method]['bias_tables_dir'],'residual/Kacprzak_2d_alpha_table.p')
		if args.calibration== 'Miller':		bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'residual/analytic_bias_table.p')
		if args.calibration== 'Miller_fid':		bias_table_string= os.path.join(config[args.method]['bias_tables_dir'],'residual/analytic_bias_table_fid.p')

	bias_table= pickle.load(open(bias_table_string, 'rb'))
	alpha_table= pickle.load(open(alpha_table_string, 'rb'))

	interpolator_m= interpolate.interp2d(bias_table[0], bias_table[1] , bias_table[2], kind='linear')
	interpolator_alpha= interpolate.interp2d(alpha_table[0], alpha_table[1] , alpha_table[2], kind='linear')
	logger.info('Constructed interpolator using %s.', bias_table_string)

	return interpolator_m, interpolator_alpha

#-----------------------------------------------------------------------------------------------------#

def get_calibration(interpolator,res):

	# Select the required columns
	if args.method=='im3shape':	snr_lab='snr' ; size_lab='mean_rgpp_rp'
	elif args.method=='ngmix':	snr_lab='snr' ; size_lab='log_T'
	
	snr= res[snr_lab] ; size= res[size_lab]

	if args.method=='ngmix':	size=np.exp(size)

	# Cycle through the entries for individual galaxies
	# Interpolate the bias table loaded above for the snr-size coordinates of each one
	b= np.array([]) 
	for i in range(len(snr)):	
		if i==0:	b= np.copy(interpolator(snr[i],size[i]))	
		else:		b= np.hstack(( b, interpolator(snr[i],size[i]) ))
	
	return b

#-----------------------------------------------------------------------------------------------------#

def relative_bias_test():
	"""
	Performs test of relative bias in ellipticity measurements from im3shape and ngmix.
	First loads both sets of results files, takes the intersection and does the same cuts on both.
	Then fits a line to the calibrated shape measurements from each code in bins of photo-z.
 	Plots this relative bias as a function of redshift.
	"""
	# Load matched results files
	global zbinning
	zbinning=True
	
	imres, ngres= load_matched_measurements(load_truth=True, load_ngmix_columns=False, load_intrinsic_shapes=False, repackage=False)

	# Apply cuts to both sets of results
	sel_ngres=[] ; sel_imres=[]
	for ir in range(len(ngres)):
		sel_ng= "sel= (ngres[ir]['flags']==0) & (np.transpose(ngres[ir]['g_sens'])[0]>0.0) & (np.transpose(ngres[ir]['g_sens'])[1]>0.0) & (np.exp(ngres[ir]['log_T_r'])/ngres[ir]['psf_T'] > 0.15) & (ngres[ir]['s2n_r']>15) & (ngres[ir]['flags']==0) "  
		sel_im= "& (imres[ir]['error_flag']==0) & (imres[ir]['snr']>15) & (imres[ir]['mean_rgpp_rp']>1.2) & (np.invert( np.isnan(imres[ir]['mean_psf_e1_sky']) ) ) & (np.invert( np.isnan(imres[ir]['mean_psf_e2_sky']) ) ) "
		sel_string= sel_ng + sel_im
		exec sel_string
		if args.selection!= 0:		sel_ngres+= [ ngres[ir][sel] ] ; sel_imres+= [ imres[ir][sel] ]
	ngres=np.array(sel_ngres) ; imres=np.array(sel_imres) 
		
	# Set up z bins and associated arrays
	zbins=np.linspace(0.3,1.3,40)
	bin_centres= (zbins[1:]+zbins[:-1])/2.0
	labcut='zphot'

	m1= [] ; errm1= []
	m2= [] ; errm2= []

	import scipy.optimize as op

	logger.info('Binning by photo-z.')
	for iz in range(len(zbins)-1):
		n1=0 ; n2=0  
		selected_imres=np.copy(imres) ; selected_ngres=np.copy(ngres) 
		for i in range(len(ngres)):
			selng= (ngres[i]['zphot']>zbins[iz]) & (ngres[i]['zphot']<zbins[iz+1])
			selim= (imres[i]['zphot']>zbins[iz]) & (imres[i]['zphot']<zbins[iz+1])  
			selected_ngres[i]= ngres[i][selng]
			selected_imres[i]= imres[i][selim]
			n1+= len(imres[i])
			n2+= len(selected_imres[i])  
		
		logger.info('z range:%2.3f-%2.3f. Contains %d objects (%f percent)',zbins[iz], zbins[iz+1], n1-n2, 100.0*n2/float(n1))

		selected_imres_all= np.array([]) ; l= 0
		selected_ngres_all=np.array([]) 
		for r in range(len(selected_imres)):	
			if l==0:
				selected_imres_all= np.copy(selected_imres[r])
				selected_ngres_all= np.copy(selected_ngres[r]) ; l= 1
			else:	
				selected_imres_all= np.hstack((selected_imres_all,selected_imres[r]))
				selected_ngres_all= np.hstack((selected_ngres_all,selected_ngres[r]))

		# Extract the data needed
		ng_e1= np.transpose(selected_ngres_all['g'])[0] ; im_e1= -1*selected_imres_all['e1']
		ng_e2= np.transpose(selected_ngres_all['g'])[1] ; im_e2= selected_imres_all['e2']
		c1= selected_imres_all['nbc_a']*selected_imres_all['mean_psf_e1_sky']
		c2= selected_imres_all['nbc_a']*selected_imres_all['mean_psf_e2_sky']

		x= (im_e1-c1)/(1+np.mean(selected_imres_all['nbc_m']))
		y= ng_e1/np.mean(np.transpose(selected_ngres_all['g_sens'])[0])

		print 'Applying calibration to im3shape results:'
		print '1+m=', (1+np.mean(selected_imres_all['nbc_m']))
		print 'c1=',np.mean(c1)
		print 'c2=',np.mean(c2)

		# Fit a line to e1^corr(im3shape) e1^corr(ngmix) 
		lab='%s_%2.3f' %(labcut, zbins[iz])
		a1,b1= op.curve_fit(func, x, y, [0,0])
		m1+=  [ a1[1] ]  ; errm1+=  [ np.sqrt(b1[1,1]) ]

		y= (im_e2-c2)/(1+np.mean(selected_imres_all['nbc_m']))
		x= ng_e2/np.mean(np.transpose(selected_ngres_all['g_sens'])[1])

		print 'Applying sensitivity corrections:'
		print 's1,2=',np.mean(np.transpose(selected_ngres_all['g_sens'])[0]), np.mean(np.transpose(selected_ngres_all['g_sens'])[1])

		# Fit a line to e2^corr(im3shape) e2^corr(ngmix)
		lab='%s_%2.3f' %(labcut, zbins[iz])
		a2,b2= op.curve_fit(func, x, y, [0,0])

		m2+=  [ a2[1] ]  ; errm2+=  [ np.sqrt(b2[1,1]) ]

	m1=np.array(m1) ; errm1=np.array(errm1)
	m2=np.array(m2) ; errm2=np.array(errm2)

	# Finally plot the relative bias as a function of z
	plt.errorbar(bin_centres, m1, yerr= errm1, fmt='.',color='b',label='e1')
	plt.plot(bin_centres, m1, color='b')
	plt.errorbar(bin_centres, m2, yerr= errm2, fmt='.',color='g',label='e2')
	plt.errorbar(bin_centres, m2, color='g')
	plt.ylim(1.01,1.09)
	plt.xlabel('photo-z')
	plt.ylabel('multiplicative bias m')
	plt.legend(loc='upper left')
	plt.xticks(zbins[0::4])
	plt.xlim(zbins[0],zbins[-1])
	plt.grid()
	plt.axhline(0,color='k')
	out_dir='im3shape-ngmix_biastest.png'
	plt.savefig(out_dir,bbox_inches='tight')
	logger.info('Saved figure to %s',out_dir)
	plt.show()

	#import pdb ; pdb.set_trace()

# Straight line to fit
def func(x, a, b):	return a+b*x

#-----------------------------------------------------------------------------------------------------#

def plot_m_z():
	"""
	Plot bias corrections to GREAT-DES in bins of photo-z.
	Not calibration scheme specific.
	"""

	logger.info('Calculating bias in bins of redshift')

	load_truth=True
	global zbinning
	zbinning=True

	global match_res ;	 
	if args.match==1:	match_res=True
	else:			match_res= False

	# Load files
	if not match_res:	
		if args.method== 'im3shape':
			res, truth_table= load_im3shape_measurements(load_truth=load_truth, load_ngmix_columns=False, load_intrinsic_shapes=True)
			selected_truth= None

		if args.method== 'ngmix':
			res, truth_table= load_ngmix_measurements(load_truth=load_truth, load_intrinsic_shapes=True)
			selected_truth= None

	elif match_res:	
			if args.method== 'ngmix':	
				imres,res= load_matched_measurements(load_truth=load_truth, load_ngmix_columns=False, load_intrinsic_shapes=True, repackage=True)
			if args.method== 'im3shape':		
				res,ngres= load_matched_measurements(load_truth=load_truth, load_ngmix_columns=False, load_intrinsic_shapes=True, repackage=True)

	#import pdb; pdb.set_trace()
	# Define the bins to use
	z_bins= np.array(config['general']['z_bins'])
	labcut='zphot'

	res_all=np.copy(res)
	m=[] ; errm=[] ; x=[]
	# Initialise the bias arrays 
	m1=[] ; errm1=[]; m2=[]; errm2=[]; m12=[] ; errm12=[] ; m21=[] ; errm21=[]
	c1=[] ; errc1=[]; c2=[]; errc2=[]; c12=[] ; errc12=[] ; c21=[] ; errc21=[]
	m1_sel= [] ; m2_sel= [] ; err_m1_sel= [] ; err_m2_sel= []
	c1_sel= [] ; c2_sel= [] ; err_c1_sel= [] ; err_c2_sel= []
	e1_int=[] ; 
	
	# Define the selection for each bin
	for iz in range(len(z_bins)-1):
		res= np.copy(res_all) ; selected_res= np.copy(res)
		logger.info('Using photo-z=%2.3f-%2.3f', z_bins[iz], z_bins[iz+1])

		alpha=[] ; beta=[] ; erralpha=[] ; errbeta=[]

		n1=0 ; n2=0
		limit1= z_bins[iz] ; limit2=-1.009804123
		# This is just here to avoid an error at the end of the loop
		try:
			limit2= z_bins[iz+1]
		except:	continue
		for i in range(len(res)):
			sel= (res[i]['zphot']>z_bins[iz]) & (res[i]['zphot']<z_bins[iz+1])
			n1+= len(res[i])
			selected_res[i]= res[i][sel]
			selected_truth= None
			if i==0:	e1_int_tmp=np.copy(selected_res[i]['e1_sheared'])
			else:		e1_int_tmp=np.hstack((e1_int_tmp, selected_res[i]['e1_sheared']))

			n2+= len(res[i][sel])
		
		logger.info('Contains %d objects (%f percent)',n2, 100.0*n2/float(n1))
		# import pdb ; pdb.set_trace()	
		# Calculate the bias for this selection
		lab='%s_%2.3f' %(labcut, z_bins[iz])
		b1,b2,b12,b21= evaluate_bias(selected_res,selected_truth,lab=lab,selection=False)
		
		m1+=  [ b1[0][0] ]  ; errm1+=  [ b1[0][1] ]
		m2+=  [ b2[0][0] ]  ; errm2+=  [ b2[0][1] ]
		m12+= [ b12[0][0] ] ; errm12+= [ b12[0][1] ]
		m21+= [ b21[0][0] ] ; errm21+= [ b21[0][1] ]	
		c1+=  [ b1[1][0] ]  ; errc1+=  [ b1[1][1] ]
		c2+=  [ b2[1][0] ]  ; errc2+=  [ b2[1][1] ]
		c12+= [ b12[1][0] ] ; errc12+= [ b12[1][1] ]
		c21+= [ b21[1][0] ] ; errc21+= [ b21[1][1] ]	
		e1_int+= [e1_int_tmp]

		if args.selection_bias==1:
			b1_sel,b2_sel,b12_sel,b21_sel= evaluate_bias(selected_res,selected_truth,lab=lab,selection=True)
			m1_sel+= [ b1_sel[0][0] ] ; err_m1_sel+=  [ b1_sel[0][1] ]
			m2_sel+= [ b2_sel[0][0] ] ; err_m2_sel+=  [ b2_sel[0][1] ]
			c1_sel+= [ b1_sel[1][0] ] ; err_c1_sel+=  [ b1_sel[1][1] ]
			c2_sel+= [ b2_sel[1][0] ] ; err_c2_sel+=  [ b2_sel[1][1] ]
			
		#import pdb ; pdb.set_trace()	
	m1=np.array(m1) ; errm1=np.array(errm1)
	m2=np.array(m2) ; errm2=np.array(errm2)
	m12=np.array(m12) ; errm12=np.array(errm12)
	m21=np.array(m21) ; errm21=np.array(errm21) 
	c1=np.array(c1) ; errc1=np.array(errc1)
	c2=np.array(c2) ; errc2=np.array(errc2)
	c12=np.array(c12) ; errc12=np.array(errc12)
	c21=np.array(c21) ; errc21=np.array(errc21) 
	e1_int=np.array(e1_int)

	z_bin_centres=(z_bins[1:]+z_bins[:-1])/2.0
	plt.subplot(3,1,1)
	plt.errorbar(z_bin_centres, m1, yerr= errm1, fmt='.',color='m',label='$m_{1}$')
	plt.errorbar(z_bin_centres, m2, yerr= errm2, fmt='.',color='r',label='$m_{2}$')
	if args.selection_bias==1:
		m1_sel= np.array(m1_sel) ; err_m1_sel= np.array(err_m1_sel)
		m2_sel= np.array(m2_sel) ; err_m2_sel= np.array(err_m2_sel)  
		plt.errorbar(z_bin_centres, m1_sel, yerr= err_m1_sel, fmt='.',color='b',label='$m_{1}$ sel. bias')
		plt.errorbar(z_bin_centres, m2_sel, yerr= err_m2_sel, fmt='.',color='g',label='$m_{2}$ sel. bias')	

	plt.xlabel('photo-z bin')
	plt.ylabel('multiplicative bias m')
	plt.legend(loc='upper right', fontsize=8)
	plt.xticks(z_bins)
	plt.xlim(z_bins[0],z_bins[-1])
	plt.axhspan(-0.03, 0.03, facecolor='0.5', alpha=0.5)
	plt.grid()
	plt.axhline(0,color='k')

	plt.subplot(3,1,2)
	plt.axhspan(-7.1e-4, 7.1e-4, facecolor='0.5', alpha=0.5)
	plt.errorbar(z_bin_centres, c1, yerr= errc1, fmt='.',color='m',label='$c_{1}$')
	plt.errorbar(z_bin_centres, c2, yerr= errc2, fmt='.',color='r',label='$c_{2}$')
	if args.selection_bias==1:
		c1_sel= np.array(c1_sel) ; err_c1_sel= np.array(err_c1_sel)
		c2_sel= np.array(c2_sel) ; err_c2_sel= np.array(err_c2_sel)  
		plt.errorbar(z_bin_centres, c1_sel, yerr= err_c1_sel, fmt='.',color='b',label='$c_{1}$ sel. bias')
		plt.errorbar(z_bin_centres, c2_sel, yerr= err_c2_sel, fmt='.',color='g',label='$c_{2}$ sel. bias')	
	plt.xlabel('photo-z bin')
	plt.ylabel('additive bias c')
	plt.legend(loc='upper right', fontsize=8)
	plt.xticks(z_bins)
	plt.xlim(z_bins[0],z_bins[-1])
	plt.grid()
	plt.axhline(0,color='k')

	plt.subplot(3,1,3)
	plt.xlabel('$e_{1}^{int}$')
	plt.xlim(-1.0,1.0)
	bins=np.linspace(-1,1,100)
	plt.hist(e1_int[2],bins=bins,histtype='step',label='upper z bin')
	plt.hist(e1_int[1],bins=bins,histtype='step',label='middle z bin')
	plt.hist(e1_int[0],bins=bins,histtype='step',label='lower z bin')
	plt.axhline(0,color='k')
	plt.legend(loc='upper right', fontsize=8)
	
	if args.output=='save_plots':
		title= os.path.join(os.path.dirname(config[args.method]['results_string']), 'out_fig/m_z_%s_unmatched.png' %args.method)
		if match_res:	title= title.replace('_unmatched.','_matched.')
		if args.selection_bias==1:	title= title.replace('.png','_selection_bias.png')
		plt.savefig(title,bbox_inches='tight')
		logger.info('Saved plot as %s',title)
	else:	plt.show()

#-----------------------------------------------------------------------------------------------------#

def main():
	valid_actions = ['get_distributions', 'get_meane_plots', 'get_single_shear_bias', 'get_bias', 'plot_bias_cuts','plot_bias_binned', 'plot_obs_int_shape', 'export_ngmix_columns', 'plot_matched_measurements', 'get_polyfit', 'plot_bias_zbinned', 'plot_m_z', 'relative_bias_test']
	valid_quantities = ['snr_r','snr_w','Tsnr','g1_sens','g2_sens','T_r']
	global logger , config , args

	description = 'Get statistics and plot results of shape measurement codes'
	parser = argparse.ArgumentParser(description=description, add_help=True)
	# Basic input and output instructions
	parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3 ), help='integer verbosity level: min=0, max=3 [default=2]')
	parser.add_argument('-a', '--actions', nargs='+' ,default=None, type=str, action='store',  help='which actions to run, available: %s' % str(valid_actions) )
	parser.add_argument('-c', '--filename_config', default='conf.yaml',type=str, action='store', help='name of the yaml config file')
	parser.add_argument('-m', '--method', default='im3shape',type=str, action='store', help='im3shape or ngmix')
	parser.add_argument('-o', '--output', default='show_plots',type=str, action='store', help='show_plots, save_plots or none')
	parser.add_argument('-N', '--number', default=600, type=int, action='store', help='Number of input results files to use.')


	# Selection options
	parser.add_argument('-q', '--quantity', nargs='+', default='none', type=str, action='store',  help='which quantity to cut on, available: %s' % str(valid_quantities))
	parser.add_argument('-s', '--selection', default=1,type=int, action='store', help='apply initial selection: 0=none, 1= flags only, 2=recommended cuts')
	parser.add_argument('-sb', '--selection_bias', default=0, type=int, action='store', help='Test for selection bias?: 0= no, 1= yes')
	parser.add_argument('-z', '--zphot', default=[0.0,0.0] , type=float, nargs=2, action='store', help='Impose photo-z cuts?: zmin,zmax')
	parser.add_argument('-w', '--weighting', default= 0 , type=int, action='store', help='Apply galaxy weights? 0= no, 1= yes')

	parser.add_argument('--calibration', default='none' , type=str, action='store', help='Calibration scheme to apply: none, Miller, Miller_fid, Kacprzak_1d, Kacprzak_2d')

	parser.add_argument('--use_truth_tables', default='0' , type=str, action='store', help='Load photo-z and intrinsic shape data from truth tables rather than pre-processed pickles. 0= no, 1= yes')
	parser.add_argument('--match', default='0' , type=int, action='store', help='Match im3shape/ngmix results, 0=no 1= yes')


	# Set up the logger
	args = parser.parse_args()
	logging_levels = { 0: logging.CRITICAL,1: logging.WARNING,2: logging.INFO,3: logging.DEBUG }
	logging_level = logging_levels[args.verbosity]; logger.setLevel(logging_level)  
	logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

	# Check the photo-z options are sensible
	global zbinning ; zbinning=False
	if (args.zphot[0]>=0.0) and (args.zphot[1]>0.0):	zbinning=True								
	if (args.zphot[0]<0.0) or (args.zphot[1]<0.0):								# Both are positive
		zbinning=False ; logger.error('One or more of the selected redshift limits is unphysical (negative). Ignoring.') 	
	if (args.zphot[0]>args.zphot[1]):									# Lower limit is higher than upper limit
		zbinning=False ; logger.error('Redshift limits not meaningful (perhaps look at the ordering). Ignoring.')
	if (args.zphot[0]>1.8):											# Lower limit is beyond the DES z range		  
		logger.error('Warning: redshift lower limit is high (z>1.8).')
		

	config = yaml.load(open(args.filename_config))	
    	if args.actions==None:
        	logger.error('no action specified, choose from %s' % valid_actions)
        	return
    	for action in valid_actions:
        	if action in args.actions:
        	    logger.info('executing %s' % action)
        	    exec action+'()'
    	for ac in args.actions:
        	if ac not in valid_actions:
        	    print 'invalid action %s' % ac
	for q in args.quantity:
        	if q not in valid_quantities:
        	    print 'invalid quantity %s' % q

 	print 'Done'
    	logger.info(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
main()

#-----------------------------------------------------------------------------------------------------#
