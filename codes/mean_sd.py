
def get_statistic(method,stat,Nran,space):

    xi0 = [] # Monopole
    xi2 = [] # Quadrupole
    xi4 = [] # Hexadecapole
    if space=='real': N = len(glob.glob(fdir+'*{}*{}*.txt'.format(method,Nran))) # glob.glob works like "ls"
    if space=='redshift':  N = len(glob.glob(fdir+'*{}*{}*redshift.txt'.format(method,Nran))) # glob.glob works like "ls"
    if space=='smallscale':  N = len(glob.glob(fdir+'*{}*{}*smallscale.txt'.format(method,Nran))) # glob.glob works like "ls"

    print('Number of files is: {} (should be 200)'.format(N))
    if method=='ran':
        for iN in range(0,N,2):
            # "xil" meaning $\xi_l$, aka multipoles
            if space=='real': xil = ascii.read(fdir+'xi_l_{}_{}_{}.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])
            if space=='redshift': xil = ascii.read(fdir+'xi_l_{}_{}_{}_redshift.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])
            if space=='smallscale': xil = ascii.read(fdir+'xi_l_{}_{}_{}_smallscale.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])

            xi0.append( xil['xi_0'].data )
            xi2.append( xil['xi_2'].data )
            xi4.append( xil['xi_4'].data )

    else:
         for iN in range(0,2*N,2):
            # "xil" meaning $\xi_l$, aka multipoles
            if space=='real': xil = ascii.read(fdir+'xi_l_{}_{}_{}-{}.txt'.format(method,Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
            if space=='redshift': xil = ascii.read(fdir+'xi_l_{}_{}_{}-{}_redshift.txt'.format(method,Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
            if space=='smallscale': xil = ascii.read(fdir+'xi_l_{}_{}_{}-{}_smallscale.txt'.format(method,Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])

            xi0.append( xil['xi_0'].data )
            xi2.append( xil['xi_2'].data )
            xi4.append( xil['xi_4'].data )

    #Calculate Standard Deviation  
    stat_xi0 = []
    stat_xi2 = []	
    stat_xi4 = []		
    for ir in range(len(nr)):
        if stat=='mean':
            stat_xi0.append( np.mean([item[ir] for item in xi0]) )
            stat_xi2.append( np.mean([item[ir] for item in xi2]) )
            stat_xi4.append( np.mean([item[ir] for item in xi4]) )
        if stat=='sd':
            stat_xi0.append( np.std([item[ir] for item in xi0],ddof=1) )
            stat_xi2.append( np.std([item[ir] for item in xi2],ddof=1) )
            stat_xi4.append( np.std([item[ir] for item in xi4],ddof=1) )

    stat_table = Table([stat_xi0, stat_xi2, stat_xi4, nr], names=('xi0', 'xi2','xi4','r'))

    ascii.write(stat_table, fname, overwrite=True)
    print('Created {}'.format(fname))


import glob 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

#nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales

stat = 'sd' # 'mean', 'sd'
method = 'ransplit' # 'ran', 'rancross', 'ransplit', 'zelrec' 
Nran = 1*87**3
space = 'smallscale' # redhshift, smallscale

#for stat in ['mean','sd']:
#    for method in ['ran']:
#        for Nran in [16*87**3,2*87**3,4*87**3,8*87**3]:
fdir = '/home/fede/Proyectos/Multipoles/data/out/{}/'.format(method) #Directorio de donde leer archivos
if space=='real': fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}.txt'.format(stat,stat,method,Nran) #Nombre del archivo que quiero crear
if space=='redshift': fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}_redshift.txt'.format(stat,stat,method,Nran) #Nombre del archivo que quiero crear
if space=='smallscale': fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}_smallscale.txt'.format(stat,stat,method,Nran) #Nombre del archivo que quiero crear
if space=='smallscale':
    nr = np.geomspace(0.5,40.,15)[:-1]
else: 
    nr = 2.5+np.linspace(5.,150.,30)[:-1]

print('file directory: ' + fdir)
print('file name for creation: ' + fname)


get_statistic(method,stat,Nran,space)
