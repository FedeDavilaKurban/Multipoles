"""
Quiero comparar los P(k)s de un zelrec, glass, y un ccvt

Lo estoy ejecutando en clemente:/mnt/is0/fdavilakurban/
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

zr=Table(np.load('zelrec_001.npy'),names=['x','y','z']) 
zcat=ArrayCatalog(zr)
zcat['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])
real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkzelrec = r.power

box=1500.
dtype=[('Position', ('f4', 3))]
gcat=BinaryCatalog('glass_001.bin',dtype,header_size=4)
gcat['Position']/=box
box=1.
glass=Table(gcat.compute(gcat['Position']),names=['x','y','z']) 
real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkglass = r.power

ccvt=ascii.read('ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
ccat=ArrayCatalog(ccvt) 
ccat['Position'] = transform.StackColumns(ccat['x'], ccat['y'], ccat['z'])
real_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkccvt = r.power

plt.loglog(Pkccvt['k'], Pkccvt['power'].real,label='CCVT')
plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass',color='green')
plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='ZelRec',color='red')
plt.legend()

plt.show()


