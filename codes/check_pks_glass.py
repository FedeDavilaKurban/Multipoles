"""
Check Pks of created Glass
Debo montar el directorio clemente
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

Nran='2x87'

zr_dir = '/home/fede/mnt/clemente/HomDist/Multipoles/data/glass/{}'.format(Nran)
plot_dir= os.path.join('/home/fede/Proyectos/Multipoles/plots/check_pks/glass/'+Nran+'/')


files = os.listdir(zr_dir)
files.sort()

box=1500.
dtype=[('Position', ('f4', 3))]


for f in files:
	
	print(files.index(f),'/',len(files)-1)

	gcat=BinaryCatalog(zr_dir+'/'+f,dtype,header_size=4)

	gcat['Position']/=1500.
	box = 1.

	real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkglass = r.power

	plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass',color='red')
	plt.ylim([1E-16,1E-4])
	plt.savefig(os.path.join(plot_dir+f[:-4]+'.png'))
	
	plt.close()

