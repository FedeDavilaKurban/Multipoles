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
import glob

box=1.
dk=0.7

nran = 1*87**3
zr_dir = '/home/fede/mnt/clemente2/8x87_testsm/'
#zr_dir = '/home/fede/mnt/clemente2/8x87/'
plot_dir= '/home/fede/Proyectos/Multipoles/plots/testsm/'



#files = glob.glob(zr_dir+'/*sm0.3*')
files = glob.glob(zr_dir+'/*')
files.sort()
#files = files[:2]

fig, axes = plt.subplots(nrows=1, ncols=10, sharey=True, figsize=(20, 4))

for i in range(10):
#for i in range(1):

	f1 = files[i]
	f2 = files[i+10]
	#f2 = files[i+1]
	sm = f1[-7:-4]
	
	print(f1,f2)
	
	zr1=Table(np.load(f1),names=['x','y','z']) 
	zcat1=ArrayCatalog(zr1)
	zcat1['Position'] = transform.StackColumns(zcat1['x'], zcat1['y'], zcat1['z'])
	real_mesh = zcat1.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkzelrec1 = r.power

	zr2 = Table(np.load(f2),names=['x','y','z']) 
	zcat2 = ArrayCatalog(zr2)
	zcat2['Position'] = transform.StackColumns(zcat2['x'], zcat2['y'], zcat2['z'])
	real_mesh = zcat2.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkzelrec2 = r.power



	axes[i].loglog(Pkzelrec1['k'], Pkzelrec1['power'].real,label='ZelRec1',color='red')
	axes[i].loglog(Pkzelrec2['k'], Pkzelrec2['power'].real,label='ZelRec2',color='red')
	axes[i].set_title('sm='+sm)
	if nran==8*87**3: axes[i].set_ylim([1E-17,1E-4])
	if nran==1*87**3: axes[i].set_ylim([1E-15,1E-4])
	

plt.tight_layout()
plt.savefig(os.path.join(plot_dir+'pks.png'))
	
plt.close()

