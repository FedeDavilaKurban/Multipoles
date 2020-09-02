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
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from pathlib import Path

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

#FOR EUCLID
#i = int(os.environ['SGE_TASK_ID'])-1
##i = missing[int(os.environ['SGE_TASK_ID'])-1]

#i=0
#print(i)

#nran = 16**3
bs = 1500.
nbins_m = 30
nbins_s = 29
nbins_ss = 14
names=['xi_0','xi_2','xi_4','xi_6']
#nbins=128
Nmesh=128 #For Plot
niter=500
#shiftmult=1.



for nran in [87**3,2*87**3,4*87**3,8*87**3]:
	print('nran=',nran)
	
	for nbins in [512]:
		print('nbins=',nbins)
		
		for shiftmult in [1.0,1.5]:
		
			print('shiftmult=',shiftmult)
		
			#for plots
			pltfolder=os.getcwd()+"/pk_{}_{}_{}".format(nran,nbins,str(shiftmult))
			Path(pltfolder).mkdir(parents=True, exist_ok=True)
			print('Folder created: '+pltfolder)
						
			#Creo el catalogo random
			np.random.seed(0)
			aran = bs*np.random.rand(nran,3)
			table_ran = Table(aran, names=['x','y','z'])


			print('Making first ZR')
			data = aran
			exec(open('zelrec_cluster.py').read())
