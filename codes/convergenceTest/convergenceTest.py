"""

"""
from __future__ import print_function, division
import os
import sys
import glob
import numpy as np
from nbodykit.lab import *
from scipy.fftpack import fftfreq
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import pyfftw
import fastmodules
from scipy.ndimage import gaussian_filter
import time
import matplotlib.pyplot as plt

t1=time.time()
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

#nran = int(sys.argv[1])*87**3 FOR EUCLID

#####ONLY TEMPORARY####
#import re
#l1=[] 
#for item in glob.glob('../data/out/zelrec/*{}*smallscale.txt'.format(nran)): 
#	num = re.search('../data/out/zelrec/xi_l_zelrec_{}_(.+?)_smallscale.txt'.format(nran),item) 
#	l1.append(int(num.group(1)))
#l1.sort() 
#l2=[]
#for item in range(0,500): 
#	l2.append(int('{:03d}'.format(item))) 
#l2.sort()

#set1=set(l1)
#set2=set(l2)
#set3=set.difference(set2,set1)
#missing=list(set3)
#missing.sort()
#########################

#FOR EUCLID
#i = int(os.environ['SGE_TASK_ID'])-1
##i = missing[int(os.environ['SGE_TASK_ID'])-1]

#i=0
#print(i)

nran = 16**3
bs = 1500.
nbins_m = 30
nbins_s = 29
nbins_ss = 14
names=['xi_0','xi_2','xi_4','xi_6']
nbins=128
Nmesh=128 #For Plot
niter=50
shiftmult=1.

#for plots
from pathlib import Path
pltfolder=os.getcwd()+"/pk_{}_{}_{}".format(nran,nbins,str(shiftmult))
Path(pltfolder).mkdir(parents=True, exist_ok=True)

datafile = '../../data/Galaxies_HOD_001_z0.57.dat'

dcat = CSVCatalog(datafile,names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z','vz'])

#Creo un arreglo redshift 'rs' y le aplico condiciones preiodicas
rs = dcat['z'].compute()+dcat['vz'].compute() 
rs = [x+1500. if x<0. else x for x in rs ]
rs = [x-1500. if x>1500. else x for x in rs ]
dcat['rs']=rs

dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['rs'])

#############################################################################
#Aca creo el catalogo random, lo divido en dos y le hago ZelRec a cada uno
############################################################################
#Creo el catalogo random
np.random.seed(0)
aran = bs*np.random.rand(nran,3)
table_ran = Table(aran, names=['x','y','z'])


#Lo divido en dos
#aran1 = aran[:int(nran/2)]
#aran2 = aran[int(nran/2):]
#print('Split files N=',len(aran1),len(aran2))
#print(aran1[:5])

print('Making first ZR')
data = aran
exec(open('zelrec_cluster.py').read())
#f1 = new
#close('zelrec_cluster.py')

# print('Making second ZR')
# data = aran2
# exec(open('zelrec_cluster.py').read())
# f2 = new
# #close('zelrec_cluster.py')

#############################################

#zr1 = Table(f1,names=['x','y','z']) 
#zr2 = Table(f2,names=['x','y','z'])
	
#zcat1 = ArrayCatalog(zr1)
#zcat1['Position'] = transform.StackColumns(zcat1['x'], zcat1['y'], zcat1['z'])
#zcat2 = ArrayCatalog(zr2)
#zcat2['Position'] = transform.StackColumns(zcat2['x'], zcat2['y'], zcat2['z'])

#Calculate and Write P(k)
#print('Calculating Pks')
# nameindex = 0
# for zcat in [zcat1,zcat2]:
# 	nameindex +=1
# 	real_mesh = zcat.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=512)
# 	r = FFTPower(real_mesh, mode='1d',dk=0.7/bs)
# 	Pkzelrec = r.power
# 	ascii.write(Table(np.column_stack((Pkzelrec['power'].real,Pkzelrec['k']))),'../data/out/zelrec/pk_{}_{}-{}.txt'.format(int(nran/2),i,nameindex),format='no_header',overwrite=True) #nran/2 porque es el espectro de uno de los zelrec, que tiene la mitad del nran que seteo arriba

# real_mesh = zcat1.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=nbins)
# r = FFTPower(real_mesh, mode='1d',dk=0.7/bs)
# Pkzelrec = r.power
# plt.loglog(Pkzelrec['k'], Pkzelrec['power']*nran/box**3)#,label='ZR1',color='red')
# plt.savefig('pk_{}_{}_{}'.format(nran,nbins,i))


#########################################
#De aca hasta abajo es para calcular y escribir la correlacion. 
#En este programa solo quiero ver como cambia el P(k) utilizando distintos parametros
#########################################
# print('Caltulating 2dCF...')
# corr = SimulationBox2PCF('2d',data1=dcat,data2=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=zcat1,randoms2=zcat2,BoxSize=[bs,bs,bs], periodic=True)

# print('Getting Multipoles')
# xi_l = get_xi0246(corr,nbins_m,nbins_s)

# wfile = '../data/out/zelrec/testLargeNbins/xi_l_zelrec_{}_{}.txt'.format(nran,i)

# print('Writing')
# ascii.write(xi_l,wfile,overwrite=True,format='no_header')


# ############################
# #Small Scale
# ############################
# #print('SMALLSCALE')
# #print('Calculating 2dCF...')
# #corr_ss = SimulationBox2PCF('2d',data1=dcat,data2=dcat,edges=np.geomspace(.5,90.,nbins_ss+1),Nmu=nbins_m,randoms1=zcat1,randoms2=zcat2,BoxSize=[bs,bs,bs], periodic=True)

# #print('Getting Multipoles')
# #xi_l_ss = get_xi0246(corr_ss,nbins_m,nbins_ss)

# #wfile_ss = '../data/out/zelrec/Niter50/xi_l_zelrec_{}_{}_smallscale.txt'.format(nran,i)

# #print('Writing')
# #ascii.write(xi_l_ss,wfile_ss,overwrite=True,format='no_header')

# #t2=time.time()
# #print('time: ',(t2-t1)/3600.,'hrs')
  
