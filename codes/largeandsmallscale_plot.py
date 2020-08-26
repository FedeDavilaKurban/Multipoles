import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

red=seaborn.color_palette()[3]
blue=seaborn.color_palette()[0]
green=seaborn.color_palette()[2]
other=seaborn.color_palette()[4]

Niter=50

ran_ss_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(2*87**3))
ran_ss_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(4*87**3))
ran_ss_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(8*87**3))
ran_ss_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(16*87**3))

ran_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(2*87**3))
ran_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(4*87**3))
ran_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(8*87**3))
ran_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(16*87**3))

zr_ss_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(2*87**3))
zr_ss_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(4*87**3))
zr_ss_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(8*87**3))
zr_ss_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale_Niter{}.txt'.format(16*87**3,Niter))

zr_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(2*87**3))
zr_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(4*87**3))
zr_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(8*87**3))
zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_Niter{}.txt'.format(16*87**3,Niter))


fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))

for ax in axes[0,0],axes[0,1]:
    x, y = zr_8['r']/.695, zr_8['xi0']#/ran_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[zr]; large scale',color=red,ls='-')

    x, y = zr_ss_8['r']/.695, zr_ss_8['xi0']#/ran_ss_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[zr]; small scale',color=red,ls='--')

    x, y = ran_8['r']/.695, ran_8['xi0']#/ran_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[ran]; large scale',color=blue,ls='-')

    x, y = ran_ss_8['r']/.695, ran_ss_8['xi0']#/ran_ss_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[ran]; small scale',color=blue,ls='--')

    ax.legend()

for ax in axes[1,0],axes[1,1]:

    x, y = zr_8['r']/.695, zr_8['xi0']/ran_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[zr/ran]; large scale',color=red,ls='-')

    x, y = zr_ss_8['r']/.695, zr_ss_8['xi0']/ran_ss_8['xi0']
    ax.plot(x,y,label=r'$\sigma$'+'[zr/ran]; small scale',color=red,ls='--')

    ax.legend()

axes[0,1].set_yscale('log')
axes[1,1].set_yscale('log')
plt.tight_layout()
plt.savefig('largeandsmallscale_plot_Niter{}.png'.format(Niter))