import numpy as np
from nbodykit.lab import *
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table


def get_xi0246(corr,nbins_m,nbins_ss):

    xi_sm = corr.corr.data['corr']
    
    ##Modificado para que lea los valores de entrada
    ##nbins_m=30 # number of bins in mu
    ##nbins_s=29 # number of bins in s
    dmu=1.0/nbins_m
    
    rs = corr.D1D2.coords['r']
    mu = corr.D1D2.coords['mu']
    
    xi_s0 = np.zeros(nbins_ss)
    xi_s2 = np.zeros(nbins_ss)
    xi_s4 = np.zeros(nbins_ss)
    xi_s6 = np.zeros(nbins_ss)
    
    sr = np.zeros(nbins_ss)
    rm = np.zeros(nbins_m)
    
    l0 = 0.0
    l1 = 1.0
    l2 = 2.0
    l3 = 3.0
    
    for i in range(nbins_ss):
    	sr[i] = rs[i]
    	for j in range(nbins_m):
    		rm[j]=mu[j]
    		xi_s0[i]  += (4.0*l0+1.0)*xi_sm[i,j]*1.0*dmu 
    		xi_s2[i]  += (4.0*l1+1.0)*xi_sm[i,j]*((3*rm[j]**2 - 1.0)/2.0)*dmu
    		xi_s4[i]  += (4.0*l2+1.0)*xi_sm[i,j]*((35*rm[j]**4 - 30*rm[j]**2 + 3.0)/8.0)*dmu
    		xi_s6[i]  += (4.0*l3+1.0)*xi_sm[i,j]*((231*rm[j]**6 - 315*rm[j]**4 + 105*rm[j]**2 - 5)/16.0)*dmu
    
    return xi_s0, xi_s2, xi_s4, xi_s6
##################################################################################

#niter = 100

bs = 1500.
#nran = .5 #*Ngal
nbins_m = 30
nbins_s = 29
nbins_ss = 14
names=['xi_0','xi_2','xi_4','xi_6']

dcat = CSVCatalog('../data/Galaxies_HOD_001_z0.57.dat',names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z','vz'])

#Creo un arreglo redshift 'rs' y le aplico condiciones preiodicas
rs = dcat['z'].compute()+dcat['vz'].compute() 
rs = [x+1500. if x<0. else x for x in rs ]
rs = [x-1500. if x>1500. else x for x in rs ]
dcat['rs']=rs

dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['rs'])

#s_array = np.linspace(5.,150.,nbins_s+1)
s_array = np.geomspace(.5,40.,nbins_ss+1)

corr_anal = SimulationBox2PCF('2d',dcat,s_array,Nmu=nbins_m,BoxSize=[bs,bs,bs], periodic=True)

xi_l = get_xi0246(corr_anal,nbins_m,nbins_ss)
    
ascii.write(xi_l,'../data/xi_l_noran_redshift_ss.txt',overwrite=True,names=['xi_0','xi_2','xi_4','xi_6'])
# #Leo los xi_l (si no los calcule ya)
xi_l =  ascii.read('../data/xi_l_noran_redshift_ss.txt',names=['xi_0','xi_2','xi_4','xi_6']) 
np.random.seed(9)
aran=bs*np.random.rand(1*dcat.size,3)                                                                                                                                                               
table_ran = Table(aran, names=['x','y','z'])    
rcat = ArrayCatalog(table_ran, BoxSize=[bs,bs,bs])
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])

corr_ran = SimulationBox2PCF('2d',dcat,s_array,randoms1=rcat,Nmu=nbins_m,BoxSize=[bs,bs,bs], periodic=True)

xi_l_ran = get_xi0246(corr_ran,nbins_m,nbins_ss)

ascii.write(xi_l_ran,'../data/xi_l_ran_redshift_ss.txt',overwrite=True,names=['xi_0','xi_2','xi_4','xi_6'])
#Leo los xi_l (si no los calcule ya)
xi_l_ran =  ascii.read('../data/xi_l_ran_redshift_ss.txt',names=['xi_0','xi_2','xi_4','xi_6']) 

nr = s_array[:-1]

plt.plot(nr,(nr**2)*xi_l['xi_0'],label=r'$\xi_0(k)$')
plt.plot(nr,(nr**2)*(-xi_l['xi_2']),label=r'$\xi_2(k)$')
plt.plot(nr,(nr**2)*xi_l['xi_4'],label=r'$\xi_4(k)$')

plt.plot(nr,(nr**2)*xi_l_ran['xi_0'],ls='--',label=r'$\xi_0 ran(k)$',color=seaborn.color_palette()[0])
plt.plot(nr,(nr**2)*(-xi_l_ran['xi_2']),ls='--',label=r'$\xi_2 ran(k)$',color=seaborn.color_palette()[1])
plt.plot(nr,(nr**2)*xi_l_ran['xi_4'],ls='--',label=r'$\xi_4 ran(k)$',color=seaborn.color_palette()[2])

plt.xscale('log')
plt.legend()
plt.savefig('../plots/xil_comparison_test.png')
plt.show()

