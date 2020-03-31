"""
Quiero ver todos los P(k)s de los ZelRec que creo para chequear que no haya alguno fulero
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


#Tengo que ubicarme en la carpeta donde esten los zelrecs

files = os.listdir('.')

for f in files:
	print(files.index(f),'/',len(files))
	zr=Table(np.load(f),names=['x','y','z']) 
	zcat=ArrayCatalog(zr)
	zcat['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])
	real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
	r = FFTPower(real_mesh, mode='1d',dk=dk)
	Pkzelrec = r.power
	plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real)

