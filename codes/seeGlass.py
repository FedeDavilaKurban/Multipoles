from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np

###
#Visualize xy projection of one Glass
###
box=1500.

dtype=[('Position', ('f4', 3))]

gcat=BinaryCatalog('glass_001.bin',dtype,header_size=4)

gcat['Position']/=box
glass=Table(gcat.compute(gcat['Position']),names=['x','y','z']) 

glass1=glass[np.where(glass['z']<.1)]

plt.scatter(glass1['x'],glass1['y'],s=.1)
plt.show()

box=1.

###
#Visualize PowSpec of one glass 
###
#gcat['Position'] = transform.StackColumns(gcat['x'], gcat['y'], gcat['z'])
dk=7./box
dk=0.7
real_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
r = FFTPower(real_mesh, mode='1d',dk=dk)
Pkglass = r.power
plt.loglog(Pkglass['k'], Pkglass['power'].real,label='Glass')

plt.show()

###
#Plot the Xi calculated using Glasses
###
names=['xi_0','xi_2','xi_4','xi_6']
nr = 2.5+np.linspace(5.,150.,30)[:-1]
nran = 8*87**3

xi_true = ascii.read('../data/xi_l_noran.txt',names=names)

for i in range(0,400,2):
    g = ascii.read('../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(nran,i,i+1),names=names)
    plt.plot(nr, g['xi_0'] )
    #plt.plot(nr, g['xi_0']-xi_true['xi_0'])
plt.yscale('log')
plt.show()

