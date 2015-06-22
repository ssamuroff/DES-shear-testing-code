import numpy as np; import matplotlib.pyplot as plt; 
from scipy.optimize import minimize as minimise 

#-----------------------------------------------------------------------------------------------------#

def get_fit_parameters(snr_bins, size_bins, data, noise):
	# First define some coordinates globally so we can access them in the fitting process below
	global size_bin_centres, snr_bin_centres, data_grid, noise_grid
	size_bin_centres= size_bins ; snr_bin_centres= snr_bins ; data_grid= data ; noise_grid= noise
	
	# Set up a starting point for the minimiser
	alpha0= 0.048 ; beta0= -2.014 ; gamma0=6.247
	fiducial_parameters= np.array([alpha0, beta0, gamma0])

	# Fit the bias model to the GREAT-DES data by minimising the rms difference
	fit_results= minimise(log_likelihood_surface, fiducial_parameters) 

	# Return an estimate for the parameter space coordinates that minimise the model-data residuals
	print 'Fit complete:'
	print fit_results
	return fit_results.x[0], fit_results.x[1], fit_results.x[2]

#-----------------------------------------------------------------------------------------------------#

def log_likelihood_surface(fit_parameters):
	"""
	Evaluates the logarithmic likelihood surface at a specified point in parameter space.
	Takes the parameter space coordinates as a numpy array. 
	Uses the same snr-size grid as the binned GREAT-DES bias calculation.  
	"""
	global snr_bin_centres, size_bin_centres, data_grid, noise_grid

	# Use the analytical bias function to get a model bias grid 
	# using the input parameters, (alpha,beta)
	model_grid=[]
	for size in size_bin_centres:	
		model_grid+= [ cfhtlens_bias(snr_bin_centres, size, fit_parameters[0], fit_parameters[1], fit_parameters[2]) ]
	model_grid= np.array(model_grid)

	# Evaluate a difference grid
	residual_grid= (model_grid- data_grid)

	# Sum over the residual grid to compute a likelihood point
	npts= (len(residual_grid)*len(residual_grid[0]))
	loglike= 0.5* sum(sum(residual_grid**2)) / sum(sum(noise_grid**2))

	return loglike

#-----------------------------------------------------------------------------------------------------#

def cfhtlens_bias(SNR, R, alpha, beta, gamma):
	"""
	Returns a multiplicative bias estimate using the analyical expression used to calibrate lensfit results. 
	See Miller et al (2013), eq 14
	"""
	#alpha= 0.057 
	#beta= -0.37
	return (beta/np.log(SNR)) * np.exp(-1.0*alpha*SNR*R**2.55) + gamma*(SNR**-2)

#-----------------------------------------------------------------------------------------------------#

def snr_basis(x, a0, a2,a4):
	#a0=0
	return a0+ a2*x**-2 + a4*x**-4


