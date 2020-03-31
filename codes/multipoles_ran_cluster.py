"""
Script para correr en el cluster MPE

Calcula multipolos con randoms, estimador L-S
"""
import os
import sys
import numpy as np
from nbodykit.lab import *
from astropy.io import ascii
from astropy import table
from astropy.table import Table

def get_xi0246(corr,nbins_m,nbins_s):

    xi_sm = corr.corr.data['corr']
    
    ##Modificado para que lea los valores de entrada
    ##nbins_m=30 # number of bins in mu
    ##nbins_s=29 # number of bins in s
    dmu=1.0/nbins_m
    
    rs = corr.D1D2.coords['r']
    mu = corr.D1D2.coords['mu']
    
    xi_s0 = np.zeros(nbins_s)
    xi_s2 = np.zeros(nbins_s)
    xi_s4 = np.zeros(nbins_s)
    xi_s6 = np.zeros(nbins_s)
    
    sr = np.zeros(nbins_s)
    rm = np.zeros(nbins_m)
    
    l0 = 0.0
    l1 = 1.0
    l2 = 2.0
    l3 = 3.0
    
    for i in range(nbins_s):
    	sr[i] = rs[i]
    	for j in range(nbins_m):
    		rm[j]=mu[j]
    		xi_s0[i]  += (4.0*l0+1.0)*xi_sm[i,j]*1.0*dmu 
    		xi_s2[i]  += (4.0*l1+1.0)*xi_sm[i,j]*((3*rm[j]**2 - 1.0)/2.0)*dmu
    		xi_s4[i]  += (4.0*l2+1.0)*xi_sm[i,j]*((35*rm[j]**4 - 30*rm[j]**2 + 3.0)/8.0)*dmu
    		xi_s6[i]  += (4.0*l3+1.0)*xi_sm[i,j]*((231*rm[j]**6 - 315*rm[j]**4 + 105*rm[j]**2 - 5)/16.0)*dmu
    
    return xi_s0, xi_s2, xi_s4, xi_s6
##################################################################################

i = int(os.environ['SGE_TASK_ID'])-1

bs = 1500.
nran = int(sys.argv[1])*87**3
nbins_m = 30
nbins_s = 29
names=['xi_0','xi_2','xi_4','xi_6']

seed_path='../data/primos.txt'
data_path='../data/Galaxies_HOD_001_z0.57.dat'

seed=np.fromfile(seed_path,dtype='int',count=1000,sep=' ')

dcat = CSVCatalog(data_path,names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])


#Calculo y escribo los xi_l
    
print '~~~~~~~~~~~~~~~~~~~~~~~~',i,'~~~~~~~~~~~~~~~~~~~~~~~~~'
np.random.seed(seed[i])
aran = bs*np.random.rand(nran,3)
table_ran = Table(aran, names=['x','y','z'])
    
rcat = ArrayCatalog(table_ran, BoxSize=[bs,bs,bs])
rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])
    
corr = SimulationBox2PCF('2d',dcat,np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat,BoxSize=[bs,bs,bs], periodic=True)

xi_l = get_xi0246(corr,nbins_m,nbins_s)
    
ascii.write(xi_l,'../data/out/ran/xi_l_{}_{}.txt'.format(nran,i),overwrite=True,format='no_header')
    

