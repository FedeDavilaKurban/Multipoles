"""
Check Pks of created ZelRecs
"""

from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np
import os

box=1500.
dk=0.7

Nran='1x87'

zr_dir = '/home/fede/mnt/clemente2/{}'.format(Nran)
plot_dir= os.path.join('/home/fdavilakurban/Proyectos/HomDist/Multipoles/plots/check_pks/zelrec/'+Nran+'/')


files = os.listdir(zr_dir)
files.sort()

for f in files[:2]:
	
	print(files.index(f),'/',len(files)-1)
	
	zr=Table(np.load(zr_dir+'/'+f),names=['x','y','z'])
	
	for col in ['x','y','z']:
		zr[col]*=box 
	zcat=ArrayCatalog(zr)
	zcat['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])
	real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
	r = FFTPower(real_mesh, mode='1d',dk=dk/box)
	Pkzelrec = r.power

	plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real*zcat.size/box**3,label='ZelRec',color='red')
	#plt.ylim([1E-16,1E-4])
	#plt.savefig(os.path.join(plot_dir+f[:-4]+'.png'))
	plt.savefig('pktest.png')
	plt.close()

	ascii.write(Table(np.column_stack((Pkzelrec['power'].real,Pkzelrec['k']))),'pktest.txt',format='no_header',overwrite=True)

