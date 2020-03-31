"""
Quiero comparar los P(k)s de un zelrec, glass, y un ccvt variando el Nran

Tengo que montar los discos euclid y clemente
"""

from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np
import os

box=1.
dk=0.7


#Nran=87**3

for Nran in [87**3,2*87**3,4*87**3,8*87**3]:

	if Nran==87**3: N='1x87'
	if Nran==2*87**3: N='2x87'
	if Nran==4*87**3: N='4x87'
	if Nran==8*87**3: N='8x87'
	
	"""
	ZELREC
	"""
	#zr_dir = '/home/fdavilakurban/mnt/euclid/Multipoles/data/zelrec/{}/zelrec_001.npy'.format(N)
	zr_dir = '/home/fdavilakurban/mnt/clemente2/{}/zelrec_001_clem.npy'.format(N)

	zr=Table(np.load(zr_dir),names=['x','y','z']) 
	zcat=ArrayCatalog(zr)
	zcat['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])
	real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkzelrec = r.power

	#plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='ZelRec',color='red')

	"""
	GLASS
	"""
	g_dir = '/home/fdavilakurban/mnt/euclid/Multipoles/data/glass/{}/glass_001.bin'.format(N)

	box=1500.
	dtype=[('Position', ('f4', 3))]
	gcat=BinaryCatalog(g_dir,dtype,header_size=4)
	gcat['Position']/=box
	box=1.
	glass=Table(gcat.compute(gcat['Position']),names=['x','y','z']) 
	real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkglass = r.power

	"""
	CCVT
	"""
#	if Nran==87**3:
#		ccvt=ascii.read('ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
#		ccat=ArrayCatalog(ccvt) 
#		ccat['Position'] = transform.StackColumns(ccat['x'], ccat['y'], ccat['z'])
#		real_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
#		r = FFTPower(real_mesh, mode='1d',dk=dk)
#		Pkccvt = r.power

#		plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
		
	plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass',color='green')
	plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='ZelRec',color='red')
	plt.legend()
	
	plt.ylim([1E-17,1E-5])	
	
	plt.savefig('../plots/Pk_comparison/Pk_{}'.format(N))
	#plt.show()
	plt.close()

