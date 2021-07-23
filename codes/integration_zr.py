#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy
from astropy.io import ascii
import seaborn
from astropy.table import Table
from nbodykit.lab import *
from scipy import integrate

# %%
rPk = ascii.read('../data/random_pk_16.dat',names=['k','P'])

filesdir='/home/fede/Proyectos/Voro/codes/zelrecon/'
pk1 = ascii.read(filesdir+'Pk_128_1.dat',names=['k','p'])
pk2 = ascii.read(filesdir+'Pk_128_2.dat',names=['k','p'])
pk3 = ascii.read(filesdir+'Pk_128_3.dat',names=['k','p'])
pk4 = ascii.read(filesdir+'Pk_128_4.dat',names=['k','p'])
pk5 = ascii.read(filesdir+'Pk_128_5.dat',names=['k','p'])
pk10 = ascii.read(filesdir+'Pk_128_10.dat',names=['k','p'])
pk20 = ascii.read(filesdir+'Pk_128_20.dat',names=['k','p'])
pk30 = ascii.read(filesdir+'Pk_128_30.dat',names=['k','p'])
pk40 = ascii.read(filesdir+'Pk_128_40.dat',names=['k','p'])
pk50 = ascii.read(filesdir+'Pk_128_50.dat',names=['k','p'])


#%%
factor = 2*np.pi*16.
meansep = (4.*np.pi*4096/3)**(-1./3)
Rs = np.geomspace(meansep/5,5*meansep,100)

#Random
P = rPk['P'][3:]
k = rPk['k'][3:]
s2_ran = []
for R in Rs:
    y = k*R
    W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
    f = (1./(2*np.pi**2)) * P * W**2 * k**2
    s2_ran.append( integrate.simps(f,k) )

#ZRs
s2_zr = []
for Pk in [pk1, pk2, pk3, pk4, pk5, pk10, pk20, pk30, pk40, pk50]:
    P = Pk['p'][3:]
    k = Pk['k'][3:]
    for R in Rs:
        y = k*R
        W = 3 * (np.sin(y)-y*np.cos(y)) / y**3
        f = (1./(2*np.pi**2)) * P * W**2 * k**2
        s2_zr.append( integrate.simps(f,k) )

s2_zr1 = s2_zr[:100]
s2_zr2 = s2_zr[100:200]
s2_zr3 = s2_zr[200:300]
s2_zr4 = s2_zr[300:400]
s2_zr5 = s2_zr[400:500]
s2_zr10 = s2_zr[500:600]
s2_zr20 = s2_zr[600:700]
s2_zr30 = s2_zr[700:800]
s2_zr40 = s2_zr[800:900]
s2_zr50 = s2_zr[900:1000]

#%%
plt.loglog(Rs,s2_ran, label='Random')

plt.loglog(Rs,s2_zr1, label='zr1')
plt.loglog(Rs,s2_zr2, label='zr2')
plt.loglog(Rs,s2_zr3, label='zr3')
plt.loglog(Rs,s2_zr4, label='zr4')
plt.loglog(Rs,s2_zr5, label='zr5')
plt.loglog(Rs,s2_zr10, label='zr10')
plt.loglog(Rs,s2_zr20, label='zr20')
plt.loglog(Rs,s2_zr30, label='zr30')
plt.loglog(Rs,s2_zr40, label='zr40')
plt.loglog(Rs,s2_zr50, label='zr50')


plt.loglog(Rs,Rs**(-3)*7E-5,ls=':',c='k',label=r'$R^{-3}$')
plt.loglog(Rs,Rs**(-4)*7E-7,ls='--',c='k',label=r'$R^{-4}$')

plt.vlines(meansep/2,5E-4,10E3,ls='-.',label='Half MeanSep')

plt.xlabel('R')
plt.ylabel(r'$\sigma^2(R)$')

plt.legend(ncol=4)
plt.tight_layout()
plt.savefig('../plots/integration_zr.png')

# %%
data = np.column_stack((Rs,s2_ran))
filename = '../data/SvR_ran.dat'
names = ['r','s2']
ascii.write(data,filename,names=names,overwrite=True)

idx = [1,2,3,4,5,10,20,30,40,50]
i = 0
for s2 in [s2_zr1, s2_zr2, s2_zr3, s2_zr4, s2_zr5, s2_zr10, s2_zr20,\
     s2_zr30, s2_zr40, s2_zr50]:
    data = np.column_stack((Rs,s2))
    filename = '../data/SvR_zel{}.dat'.format(str(idx[i]))
    names = ['r','s2']
    ascii.write(data,filename,names=names,overwrite=True)
    i+=1


# %%
