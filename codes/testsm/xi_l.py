from astropy.io import ascii                                              
from astropy import table
from astropy.table import Table
from nbodykit.lab import *
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

names = ['xi0','xi2','xi4','xi6']
plot_dir= os.path.join('/home/fede/Proyectos/Multipoles/plots/testsm/')

nran = 8*87**3

xil_01 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.1.txt'.format(nran),names=names)
xil_02 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.2.txt'.format(nran),names=names)
#xil_03 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.3.txt'.format(nran),names=names)
xil_03 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_198-199_sm0.3.txt'.format(nran),names=names)
xil_04 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.4.txt'.format(nran),names=names)
xil_05 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.5.txt'.format(nran),names=names)
xil_06 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.6.txt'.format(nran),names=names)
xil_07 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.7.txt'.format(nran),names=names)
xil_08 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.8.txt'.format(nran),names=names)
xil_09 = ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm0.9.txt'.format(nran),names=names)
xil_1 =  ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_testsm1.0.txt'.format(nran),names=names)

xil =    ascii.read('/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/xi_l_zelrec_{}_0-1_clem.txt'.format(nran),names=names)
xil_true = ascii.read('/home/fede/Proyectos/Multipoles/data/xi_l_noran.txt',names=names)

nr = 2.5+np.linspace(5.,150.,30)[:-1]

"""
Puedo hacer la resta con un xi de referencia como el xil_1, o la true
"""

fig, axes = plt.subplots(nrows=1, ncols=10, sharey=True, figsize=(20, 4))

axes[0].plot(nr,xil_01['xi0']-xil_true['xi0'])
axes[1].plot(nr,xil_02['xi0']-xil_true['xi0'])
axes[2].plot(nr,xil_03['xi0']-xil_true['xi0'])
axes[3].plot(nr,xil_04['xi0']-xil_true['xi0'])
axes[4].plot(nr,xil_05['xi0']-xil_true['xi0'])
axes[5].plot(nr,xil_06['xi0']-xil_true['xi0'])
axes[6].plot(nr,xil_07['xi0']-xil_true['xi0'])
axes[7].plot(nr,xil_08['xi0']-xil_true['xi0'])
axes[8].plot(nr,xil_09['xi0']-xil_true['xi0'])
axes[9].plot(nr,xil_1['xi0']-xil_true['xi0'])
plt.tight_layout()
plt.savefig(os.path.join(plot_dir+'xi0test-xi0true1.png'))
plt.close()
