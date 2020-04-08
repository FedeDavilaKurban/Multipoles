"""
Hay una version anterior que hace lo mismo, se llama compare_pks_differentN.py

Voy a seguir usando esta
"""

from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np
import os

Nran='1x87'

zrdir = '/home/fede/mnt/clemente2/'
glassdir = '/home/fede/mnt/euclid/Multipoles/data/glass/'
ccvtdir = '/home/fede/mnt/clemente2/'


box=1.
dk=0.7

for Nran in ['1x87','2x87','4x87','8x87']:
	zr=Table(np.load(zrdir+Nran+'/'+'zelrec_001_clem.npy'),names=['x','y','z']) 
	zcat=ArrayCatalog(zr)
	zcat['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])
	real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkzelrec = r.power

	box=1500.
	dtype=[('Position', ('f4', 3))]
	gcat=BinaryCatalog(glassdir+Nran+'/'+'glass_001.bin',dtype,header_size=4)
	gcat['Position']/=box
	box=1.
	glass=Table(gcat.compute(gcat['Position']),names=['x','y','z']) 
	real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkglass = r.power

	if Nran=='1x87':
		ccvt=ascii.read(ccvtdir+'ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
		ccat=ArrayCatalog(ccvt) 
		ccat['Position'] = transform.StackColumns(ccat['x'], ccat['y'], ccat['z'])
		real_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
		r = FFTPower(real_mesh, mode='1d',dk=dk)
		Pkccvt = r.power

		plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')

	plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass',color='green')
	plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='ZelRec',color='red')
	plt.xlim(10,2000)
	plt.legend()
	plt.ylim([1E-17,1E-5])	
	plt.savefig('/home/fede/Proyectos/Multipoles/plots/compare_pks/{}'.format(Nran))
	plt.close()
#plt.show()


