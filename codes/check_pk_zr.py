"""
Check Pks of created ZelRecs
"""
import os
import glob
from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np

box=1500.
dk=0.7

#Nran='1x87'
nran=2*87**3

zr_dir = '/home/fede/Proyectos/Multipoles/data/out/zelrec/pks/'
plot_dir= os.path.join('/home/fede/Proyectos/Multipoles/plots/check_pks/zelrec/')

files = glob.glob(zr_dir+'*{}*'.format(nran))
files.sort()

for f in files:
    print(f)
    #print(files.index(f),'/',len(files)-1)
    zr=ascii.read(os.path.join(f),names=['pk','k'])

    plt.loglog(zr['k'], zr['pk']*nran/box**3)#,label='ZR1',color='red')
    plt.ylim(1E-9,1E1)
    plt.savefig(os.path.join(plot_dir+'pk_{}_{}.png'.format(nran,files.index(f))))
    plt.close()

