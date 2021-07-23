#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy
from astropy.io import ascii
import seaborn
from astropy.table import Table
from nbodykit.lab import *
from scipy import integrate
#%%
"""
Test
"""
P = 1
k = np.arange(1,1000)
R = 1./40

y = k*R
W = 3 * (np.sin(y)-y*np.cos(y)) / y**3

f = (1./(2*np.pi**2)) * P * W**2 * k**2

plt.plot(k,f)
# %%
gPk = ascii.read('../data/glass_pk_16.dat',names=['k','P'])
rPk = ascii.read('../data/random_pk_16.dat',names=['k','P'])
cPk = ascii.read('../data/ccvt_pk_16.dat',names=['k','P'])
zPk = ascii.read('../data/glass_pk_16.dat',names=['k','P'])
zPk = ascii.read('../../Voro/codes/zelrecon/Pk_128_50.dat',names=['k','P'])

#%%
factor = 2*np.pi*16.
meansep = (4.*np.pi*4096/3)**(-1./3)
#meansep = 1./16
Rs = np.geomspace(0.0004 ,10*meansep,500)

#Random
P = rPk['P'][3:]
k = rPk['k'][3:]
s2_ran = []
for R in Rs:
    y = k*R
    W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
    f = (1./(2*np.pi**2)) * P * W**2 * k**2
    s2_ran.append( integrate.simps(f,k) )

#Grav Glass
P = gPk['P'][3:]
k = gPk['k'][3:]
s2_gg = []
for R in Rs:
    y = k*R
    W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
    f = (1./(2*np.pi**2)) * P * W**2 * k**2
    s2_gg.append( integrate.simps(f,k) )

#CCVT
P = cPk['P'][3:]
k = cPk['k'][3:]
s2_ccvt = []
for R in Rs:
    y = k*R
    W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
    f = (1./(2*np.pi**2)) * P * W**2 * k**2
    s2_ccvt.append( integrate.simps(f,k) )

#ZR
P = zPk['P'][3:]
k = zPk['k'][3:]
s2_z = []
for R in Rs:
    y = k*R
    W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
    f = (1./(2*np.pi**2)) * P * W**2 * k**2
    s2_z.append( integrate.simps(f,k) )


Rs*=1./meansep

plt.loglog(Rs,s2_ran, label='Random')
plt.loglog(Rs,s2_gg, label='Grav. Glass')
plt.loglog(Rs,s2_ccvt, label='CCVT')
plt.loglog(Rs,s2_z, label='Zeldovich')


#plt.loglog(Rs,Rs**(-3)*2.8E-1,ls=':',c='k',label=r'$R^{-3}$')
#plt.loglog(Rs,Rs**(-4)*.5E-1,ls='--',c='k',label=r'$R^{-4}$')

#meansep= .5
plt.vlines(.5,5E-4,10E3,ls='-.',label='Half MeanSep')

plt.xlabel('R')
plt.ylabel(r'$\sigma^2(R)$')

plt.legend(ncol=2)
plt.tight_layout()
#plt.savefig('../plots/integration.png')

# %%
data = np.column_stack((Rs,s2_ran))
filename = '../data/SvR_ran.dat'
names = ['r','s2']
ascii.write(data,filename,names=names,overwrite=True)

data = np.column_stack((Rs,s2_gg))
filename = '../data/SvR_glass.dat'
names = ['r','s2']
ascii.write(data,filename,names=names,overwrite=True)

data = np.column_stack((Rs,s2_ccvt))
filename = '../data/SvR_ccvt.dat'
names = ['r','s2']
ascii.write(data,filename,names=names,overwrite=True)

data = np.column_stack((Rs,s2_z))
filename = '../data/SvR_zel.dat'
names = ['r','s2']
ascii.write(data,filename,names=names,overwrite=True)

# %%

# %%
