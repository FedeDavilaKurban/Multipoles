
def get_statistic(method,stat,Nran,space):

    xi0 = [] # Monopole
    xi2 = [] # Quadrupole
    xi4 = [] # Hexadecapole

    print('Number of files is: {}'.format(N))

    if space=='redshift':
        for iN in range(0,N):
            xil = ascii.read(fdir+'xi_l_{}_{}_{}.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])

            xi0.append( xil['xi_0'].data )
            xi2.append( xil['xi_2'].data )
            xi4.append( xil['xi_4'].data )

    if space=='smallscale': 
        for iN in range(0,N):
            xil = ascii.read(fdir+'xi_l_{}_{}_{}_smallscale.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])

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

#stat = 'sd' # 'mean', 'sd'
#method = 'ran' # 'ran', 'rancross', 'ransplit', 'zelrec' 
#Nran = 2*87**3 # (2,4,8,16) *87**3
#space = 'smallscale' # redhshift, smallscale

N = 500 #Num de archivos de cada metodo
Niter = 50 #Num de iteraciones de la reconstruccion de zeldovich. Cambia un poco el Pk de los ZelRec, nada mas.

for space in ['redshift','smallscale']:
    for stat in ['mean','sd']:
        for method in ['zelrec']:#,'ran','rancross','ransplit']:
            for Nran in [16*87**3]:#,2*87**3,4*87**3,8*87**3]:

                fdir = '/home/fede/Proyectos/Multipoles/data/out/{}/'.format(method) #Directorio de donde leer archivos
                if method=='zelrec': fdir = '/home/fede/Proyectos/Multipoles/data/out/Niter{}/'.format(Niter)

                if space=='redshift': 
                    fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}_Niter{}.txt'.format(stat,stat,method,Nran,Niter) #Nombre del archivo que quiero crear
                    nr = np.linspace(5.,150.,30)[:-1]

                if space=='smallscale': 
                    fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}_smallscale_Niter{}.txt'.format(stat,stat,method,Nran,Niter) #Nombre del archivo que quiero crear
                    #nr = np.geomspace(0.5,40.,15)[:-1]
                    nr = np.geomspace(0.5,90.,15)[:-1]

                print('file directory: ' + fdir)
                print('file name for creation: ' + fname)


                get_statistic(method,stat,Nran,space)
