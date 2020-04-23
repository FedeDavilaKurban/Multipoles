
def get_statistic(method,stat,Nran):

    xi0 = [] # Monopole
    xi2 = [] # Quadrupole
    xi4 = [] # Hexadecapole
    N = len(glob.glob(fdir+'*{}*{}*.txt'.format(method,Nran))) # glob.glob works like "ls"
    print('Number of xi_zelrec is: {} (should be 200)'.format(N))
    if method=='ran':
        for iN in range(0,N,2):
            # "xil" meaning $\xi_l$, aka multipoles
            xil = ascii.read(fdir+'xi_l_{}_{}_{}.txt'.format(method,Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])
            xi0.append( xil['xi_0'].data )
            xi2.append( xil['xi_2'].data )
            xi4.append( xil['xi_4'].data )

    else:
         for iN in range(0,2*N,2):
            # "xil" meaning $\xi_l$, aka multipoles
            xil = ascii.read(fdir+'xi_l_{}_{}_{}-{}.txt'.format(method,Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
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

nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales

stat = 'sd'
method = 'ccvt'
Nran = 1*87**3

fdir = '/home/fede/Proyectos/Multipoles/data/out/{}/'.format(method) #Directorio de donde leer archivos
fname = '/home/fede/Proyectos/Multipoles/data/out/{}_xi/{}_{}_{}.txt'.format(stat,stat,method,Nran) #Nombre del archivo que quiero crear
print('file directory: ' + fdir)
print('file name for creation: ' + fname)


get_statistic(method,stat,Nran)
