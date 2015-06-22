import pyfits
import numpy as np

cat=pyfits.getdata('/home/samuroff/GREAT-DES/truth/cosmos_shapes/real_galaxy_23.5_shapes.fits')
l='home/samuroff/GREAT-DES/results/ngmix/v006/e02/sfit-e02-g%02d.fits'

ig=0

for ig in range(7):

	res=pyfits.getdata(l%ig)
	fil='nbc.truth.%03d.'+'g%02d.fits.gz'%ig
	tr=np.array([])
	for i in range(600):
	     f=fil%i
	     if i==0:        tr=np.copy(tabletools.loadTable(f)); d=tr.dtype
	     else:           tr=np.hstack((tr,tabletools.loadTable(f).astype(d)))

	shape_e1= cat[tr['id_cosmos']]['e1']
	shape_e2= cat[tr['id_cosmos']]['e2']

	eps= shape_e1+ shape_e2*1j
	q= np.sqrt( (1-np.abs(eps)) / (1+np.abs(eps)) )
	e= (1-q)/(1+q) * np.exp(1j*np.angle(eps))
	e1= e.real ; e2= e.imag

	e_int= e*np.exp(2*1j*tr['rotation_angle'])
	e1_int= e.real ; e2_int= e.imag

	g1= tr['g1_true'] ; g2= tr['g2_true']
	g= g1+1j*g2

	e_sheared= (e_int+g)/(1+g.conjugate()*e_int) 
	e1_sheared= e_sheared.real ; e2_sheared= e_sheared.imag

	z= tr['zphot']
	import cPickle as pickle
	pickle.dump([e1_sheared, e2_sheared,z],open('sfit-e02-g%02d-intrinsic-shapes.p'%ig,'wb'))

	print 'Wrote shear group %d columns to %s' %(ig,'sfit-e02-g%02d-intrinsic-shapes.p'%ig )
