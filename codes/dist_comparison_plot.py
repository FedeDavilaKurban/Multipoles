import numpy as np
import matplotlib.pyplot as plt
import scipy
from astropy.io import ascii
import seaborn
from astropy.table import Table
from nbodykit.lab import *

np.random.seed(1)

glass = ascii.read('../data/glass_16.dat',names=['x','y','z'])
ccvt = ascii.read('../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
ran = Table(np.random.random((len(glass),3)),names=['x','y','z'])
zr = ascii.read('/home/fede/Proyectos/Voro/codes/zelrecon/zelrec_final.dat',names=['x','y','z'])

glass_ = glass[np.where(glass['z']<.1)]
ccvt_ = ccvt[np.where(ccvt['z']<.1)]
ran_ = ran[np.where(ran['z']<.1)]
zr_ = zr[np.where(zr['z']<.1)]

fig, ax = plt.subplots(1, 4, figsize=(17,4), sharex='col', sharey='row')

fs=20
lfs=14

ax[0].scatter(ran_['x'],ran_['y'])
ax[1].scatter(glass_['x'],glass_['y'],color=seaborn.color_palette()[2])
ax[2].scatter(ccvt_['x'],ccvt_['y'],color=seaborn.color_palette()[1])
ax[3].scatter(zr_['x'],zr_['y'],color=seaborn.color_palette()[3])

ax[0].set_title('a) Random',fontsize=fs)
ax[1].set_title('b) Gravitational Glass',fontsize=fs)
ax[2].set_title('c) CCVT',fontsize=fs)
ax[3].set_title('d) Zeldovich Reconstruction',fontsize=fs)

#aumento el tamaño de los ticks
ax[0].tick_params(labelsize=fs-4,axis='both')
ax[1].tick_params(labelsize=fs-4,axis='both')
ax[2].tick_params(labelsize=fs-4,axis='both')
ax[3].tick_params(labelsize=fs-4,axis='both')


#seteo a 3 el numero de ticks
ax[0].set_xticks(np.linspace(0, 1, 3))
ax[1].set_xticks(np.linspace(0, 1, 3))
ax[2].set_xticks(np.linspace(0, 1, 3))
ax[3].set_xticks(np.linspace(0, 1, 3))

ax[0].set_yticks(np.linspace(0, 1, 3))
ax[1].set_yticks(np.linspace(0, 1, 3))
ax[2].set_yticks(np.linspace(0, 1, 3))
ax[3].set_xticks(np.linspace(0, 1, 3))

for i in range(4):
    ax[i].set_ylim([0,1])
    ax[i].set_xlim([0,1])

plt.tight_layout()
plt.savefig('../plots/dist_comparison_plot.png')
plt.savefig('../plots/dist_comparison_plot.pdf')
plt.show()
plt.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

msep = 1./len(glass)
dk = 2

gcat = ArrayCatalog(glass,names=['x','y','z'])      
gcat['Position'] = transform.StackColumns(gcat['x'], gcat['y'], gcat['z'])
g_mesh = gcat.to_mesh(compensated=True, resampler='tsc', position='Position', Nmesh=256, BoxSize=1)
gr = FFTPower(g_mesh, mode='1d',dk=dk)
gPk = gr.power

rcat = ArrayCatalog(ran,names=['x','y','z'])      
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
r_mesh = rcat.to_mesh(compensated=True, resampler='tsc', position='Position', Nmesh=256, BoxSize=1)
rr = FFTPower(r_mesh, mode='1d',dk=dk)
rPk = rr.power

ccat = ArrayCatalog(ccvt,names=['x','y','z'])      
ccat['Position'] = transform.StackColumns(ccat['x'], ccat['y'], ccat['z'])
c_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', Nmesh=256, BoxSize=1)
cr = FFTPower(c_mesh, mode='1d',dk=dk)
cPk = cr.power

plt.figure(figsize=(7,5))
factor = 16.*np.pi*2
plt.loglog(rPk['k']/factor, rPk['power'].real*len(ran),label='Random')
plt.loglog(gPk['k']/factor, gPk['power'].real*len(glass),label='Glass',color=seaborn.color_palette()[2],ls='-.')
plt.loglog(cPk['k']/factor, cPk['power'].real*len(ccvt),label='CCVT',color=seaborn.color_palette()[1],dashes=[5,2])

gk=gPk['k'][:-200]
plt.hlines(1., 6/factor,rPk['k'][-1]/factor,linestyles=':',label=r"$P(k) = 1/\bar{n}$")
plt.loglog(gk/factor,10**8*(gk*msep/(2*np.pi))**4,dashes=[3,3,2,2],c='k',alpha=.8,label=r"$P(k)\propto  k^{4}$") 


plt.xlabel(r"$k/(2\pi\bar{n}^{-1/3})$",fontsize=fs)
plt.ylabel(r"$P(k)*\bar{n}$",fontsize=fs)
plt.legend(fontsize=fs-3,markerscale=1)
plt.yticks(ticks=[1,10E-7,10E-4])
plt.xlim([5E-2,1E1])
plt.tick_params(axis='both',labelsize=fs-4)
plt.tight_layout()
plt.savefig('../plots/pk_comparison_plot_sq.png')
plt.savefig('../plots/pk_comparison_plot_sq.pdf')
plt.show()
plt.close()

wdata = Table(np.column_stack((gPk['k'],gPk['power'].real)))
ascii.write(wdata,'../data/glass_pk_16.dat',names=['k','P'])

wdata = Table(np.column_stack((rPk['k'],rPk['power'].real)))
ascii.write(wdata,'../data/random_pk_16.dat',names=['k','P'])

wdata = Table(np.column_stack((cPk['k'],cPk['power'].real)))
ascii.write(wdata,'../data/ccvt_pk_16.dat',names=['k','P'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filesdir='/home/fede/Proyectos/Voro/codes/zelrecon/'
colors = seaborn.color_palette("RdBu_r", 5)
colors = seaborn.color_palette("coolwarm", 20)
#colors = colors[-6:]

#Nbins 128
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
pk70 = ascii.read(filesdir+'Pk_128_70.dat',names=['k','p'])
pk80 = ascii.read(filesdir+'Pk_128_80.dat',names=['k','p'])
pk90 = ascii.read(filesdir+'Pk_128_90.dat',names=['k','p'])
pk100 = ascii.read(filesdir+'Pk_128_100.dat',names=['k','p'])
pk200 = ascii.read(filesdir+'Pk_128_200.dat',names=['k','p'])
pk300 = ascii.read(filesdir+'Pk_128_300.dat',names=['k','p'])
pk400 = ascii.read(filesdir+'Pk_128_400.dat',names=['k','p'])
pk500 = ascii.read(filesdir+'Pk_128_500.dat',names=['k','p'])
pk600 = ascii.read(filesdir+'Pk_128_600.dat',names=['k','p'])
pk700 = ascii.read(filesdir+'Pk_128_700.dat',names=['k','p'])
pk800 = ascii.read(filesdir+'Pk_128_800.dat',names=['k','p'])
pk900 = ascii.read(filesdir+'Pk_128_900.dat',names=['k','p'])
pk1000 = ascii.read(filesdir+'Pk_128_1000.dat',names=['k','p'])
pk1500 = ascii.read(filesdir+'Pk_128_1500.dat',names=['k','p'])
pk2000 = ascii.read(filesdir+'Pk_128_2000.dat',names=['k','p'])
pk2500 = ascii.read(filesdir+'Pk_128_2500.dat',names=['k','p'])
pk3000 = ascii.read(filesdir+'Pk_128_3000.dat',names=['k','p'])

zr = ascii.read(filesdir+'Pk_zelrec.dat',names=['k','p'])

#Nbins256
pk1_256 = ascii.read(filesdir+'Pk_256_1.dat',names=['k','p'])
pk2_256 = ascii.read(filesdir+'Pk_256_2.dat',names=['k','p'])
pk3_256 = ascii.read(filesdir+'Pk_256_3.dat',names=['k','p'])
pk4_256 = ascii.read(filesdir+'Pk_256_4.dat',names=['k','p'])
pk5_256 = ascii.read(filesdir+'Pk_256_5.dat',names=['k','p'])
pk10_256 = ascii.read(filesdir+'Pk_256_10.dat',names=['k','p'])
pk50_256 = ascii.read(filesdir+'Pk_256_50.dat',names=['k','p'])
pk100_256 = ascii.read(filesdir+'Pk_256_100.dat',names=['k','p'])
pk200_256 = ascii.read(filesdir+'Pk_256_200.dat',names=['k','p'])
pk300_256 = ascii.read(filesdir+'Pk_256_300.dat',names=['k','p'])
pk400_256 = ascii.read(filesdir+'Pk_256_400.dat',names=['k','p'])
pk500_256 = ascii.read(filesdir+'Pk_256_500.dat',names=['k','p'])
pk600_256 = ascii.read(filesdir+'Pk_256_600.dat',names=['k','p'])
pk700_256 = ascii.read(filesdir+'Pk_256_700.dat',names=['k','p'])
pk800_256 = ascii.read(filesdir+'Pk_256_800.dat',names=['k','p'])
pk900_256 = ascii.read(filesdir+'Pk_256_900.dat',names=['k','p'])
pk1000_256 = ascii.read(filesdir+'Pk_256_1000.dat',names=['k','p'])

plt.figure(figsize=(7,5))


#plt.loglog(gk,9.5**8*(gk*msep/(2*np.pi))**4,ls='-.',c='k',alpha=.8) 
#plt.loglog(gk,9**8*(gk*msep/(2*np.pi))**4,ls='-.',c='k',alpha=.8) 
#plt.loglog(gk,8.5**8*(gk*msep/(2*np.pi))**4,ls='-.',c='k',alpha=.8) 

#plt.loglog(cPk['k'], cPk['power'].real*len(ccvt),ls=':',color=seaborn.color_palette()[1])
#plt.loglog(gPk['k'], gPk['power'].real*len(glass),ls=':',color=seaborn.color_palette()[2])

#plt.loglog(pk50['k'], pk50['p']*16**3,label='Pk50')

plt.loglog(rPk['k']/factor, rPk['power'].real*len(ran),label='Random')
plt.loglog(pk1['k']/factor, pk1['p']*16**3,color=colors[0],ls='-')
plt.loglog(pk2['k']/factor, pk2['p']*16**3,color=colors[1],ls='-')
plt.loglog(pk3['k']/factor, pk3['p']*16**3,color=colors[2],ls='-')
plt.loglog(pk4['k']/factor, pk4['p']*16**3,color=colors[3],ls='-')
plt.loglog(pk5['k']/factor, pk5['p']*16**3,color=colors[15],ls='-')
plt.loglog(pk10['k']/factor, pk10['p']*16**3,color=colors[16],ls='-')
plt.loglog(pk20['k']/factor, pk20['p']*16**3,color=colors[17],ls='-')
plt.loglog(pk30['k']/factor, pk30['p']*16**3,color=colors[18],ls='-')
plt.loglog(pk40['k']/factor, pk40['p']*16**3,color=colors[19],ls='-')

# plt.loglog(pk70['k'], pk70['p']*16**3,label='Pk70')
# plt.loglog(pk80['k'], pk80['p']*16**3,label='Pk80')
# plt.loglog(pk100['k'], pk100['p']*16**3,label='Pk100')
# plt.loglog(pk200['k'], pk200['p']*16**3,label='Pk200')
# plt.loglog(pk300['k'], pk300['p']*16**3,label='Pk300')
# plt.loglog(pk500['k'], pk500['p']*16**3,label='Pk500')
# plt.loglog(pk1000['k'], pk1000['p']*16**3,label='Pk1000')
# plt.loglog(pk1500['k'], pk1500['p']*16**3,label='Pk1500')
# plt.loglog(pk2000['k'], pk2000['p']*16**3,label='Pk2000')
# plt.loglog(pk2500['k'], pk2500['p']*16**3,label='Pk2500')
# plt.loglog(pk3000['k'], pk3000['p']*16**3,label='Pk3000')

# plt.loglog(pk50_256['k'], pk50_256['p']*16**3,label='Pk50_256')
# plt.loglog(pk100_256['k'], pk100_256['p']*16**3,label='Pk100_256')
# plt.loglog(pk200_256['k'], pk200_256['p']*16**3,label='Pk200_256')
# plt.loglog(pk300_256['k'], pk300_256['p']*16**3,label='Pk300_256')
# plt.loglog(pk400_256['k'], pk400_256['p']*16**3,label='Pk400_256')
# plt.loglog(pk500_256['k'], pk500_256['p']*16**3,label='Pk500_256')
# plt.loglog(pk600_256['k'], pk600_256['p']*16**3,label='Pk600_256')
# plt.loglog(pk700_256['k'], pk700_256['p']*16**3,label='Pk700_256')
# plt.loglog(pk800_256['k'], pk800_256['p']*16**3,label='Pk800_256')
# plt.loglog(pk900_256['k'], pk900_256['p']*16**3,label='Pk900_256')
# plt.loglog(pk1000_256['k'], pk1000_256['p']*16**3,label='Pk1000_256')

plt.loglog(pk50['k']/factor, pk50['p']*16**3,label='Zel. Rec.',color=seaborn.color_palette()[3])

plt.loglog(gk/factor,10**8*(gk*msep/(2*np.pi))**4,dashes=[3,3,2,2],c='k',alpha=.8,label=r"$P(k)\propto  k^{4}$") 
plt.hlines(1., 6/factor,rPk['k'][-1]/factor,linestyles=':',label=r"$P(k)=1/\bar{n}$")

plt.xlabel(r"$k/(2\pi\bar{n}^{-1/3})$",fontsize=fs)
plt.ylabel(r"$P(k)*\bar{n}$",fontsize=fs)
plt.legend(fontsize=fs-3,markerscale=1)
plt.yticks(ticks=[1,10E-7,10E-4])
plt.xlim([5E-2,1E1])
plt.tick_params(axis='both',labelsize=fs-4)
plt.tight_layout()
plt.savefig('../plots/pk_zr_iter.png')
plt.savefig('../plots/pk_zr_iter.pdf')
plt.show()
plt.close()
