from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np

###
#Visualize xy projection of one Zelrec
###
box=1500.

zcat=Table(np.load("zelrec_001_clem.npy",allow_pickle=True),names=['x','y','z'])
for c in ['x','y','z']:
	zcat[c]*=box

zcat1=zcat[np.where(zcat['z']<150.)]

plt.scatter(zcat1['x'],zcat1['y'],s=.1)
plt.show()


###
#Visualize PowSpec of one Zelrec
###
zelrec = ArrayCatalog(zcat, BoxSize=[1,1,1])
zelrec['Position'] = transform.StackColumns(zcat['x'], zcat['y'], zcat['z'])

dk=7./box
real_mesh = zelrec.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkzelrec = r.power
plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='Zelrec')
plt.show()

###
#Comparison of P(k)s with glass and ccvt, both of 87**3
###

dtype=[('Position', ('f4', 3))]
gcat=BinaryCatalog('../data/glass_001.bin',dtype,header_size=4)
#gcat['Position']/=1500.
real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkglass = r.power

rcat = CSVCatalog('../data/ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
rcat['Position']*=box
real_mesh = rcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkccvt = r.power

plt.loglog(Pkzelrec['k'], Pkzelrec['power'].real,label='Zelrec')
plt.loglog(Pkzelrec['k'], Pkccvt['power'].real,label='CCVT')
plt.loglog(Pkzelrec['k'], Pkglass['power'].real,label='Glass')
plt.legend()
plt.show()


###
#Plot the Xi calculated using Zelrecs
###
names=['xi_0','xi_2','xi_4','xi_6']
nr = 2.5+np.linspace(5.,150.,30)[:-1]
nran = 87**3

for i in range(0,200,2):
    z = ascii.read('../data/out/zelrec/xi_l_zelrec_{}_{}-{}.txt'.format(nran,i,i+1),names=names)
    plt.plot(nr, z['xi_0'] )
plt.yscale('log')
plt.show()


