import glob 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table


#Some Plot Parameters
fs=20
lfs=14
red=seaborn.color_palette()[3]
blue=seaborn.color_palette()[0]
green=seaborn.color_palette()[2]
other=seaborn.color_palette()[4]

Niter = 100

#######################################################################
#Reading
#######################################################################
#Ran
ran_1 = ascii.read('../data/out/mean_xi/mean_ran_{}.txt'.format(2*87**3))
ran_2 = ascii.read('../data/out/mean_xi/mean_ran_{}.txt'.format(4*87**3))
ran_4 = ascii.read('../data/out/mean_xi/mean_ran_{}.txt'.format(8*87**3))
ran_8 = ascii.read('../data/out/mean_xi/mean_ran_{}.txt'.format(16*87**3))

#RanCross
rc_1 = ascii.read('../data/out/mean_xi/mean_rancross_{}.txt'.format(2*87**3))
rc_2 = ascii.read('../data/out/mean_xi/mean_rancross_{}.txt'.format(4*87**3))
rc_4 = ascii.read('../data/out/mean_xi/mean_rancross_{}.txt'.format(8*87**3))
rc_8 = ascii.read('../data/out/mean_xi/mean_rancross_{}.txt'.format(16*87**3))

#RanSplit
rs_1 = ascii.read('../data/out/mean_xi/mean_ransplit_{}.txt'.format(2*87**3))
rs_2 = ascii.read('../data/out/mean_xi/mean_ransplit_{}.txt'.format(4*87**3))
rs_4 = ascii.read('../data/out/mean_xi/mean_ransplit_{}.txt'.format(8*87**3))
rs_8 = ascii.read('../data/out/mean_xi/mean_ransplit_{}.txt'.format(16*87**3))

#ZelRec
zr_1 = ascii.read('../data/out/mean_xi/mean_zelrec_{}.txt'.format(2*87**3))
zr_2 = ascii.read('../data/out/mean_xi/mean_zelrec_{}.txt'.format(4*87**3))
zr_4 = ascii.read('../data/out/mean_xi/mean_zelrec_{}.txt'.format(8*87**3))
zr_8 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_Niter{}.txt'.format(16*87**3,Niter))
#zr_8 = ascii.read('../data/out/mean_xi/mean_zelrec_{}.txt'.format(16*87**3))

#Analytic
xi_noran = ascii.read('../data/xi_l_noran_redshift.txt',names=['xi0','xi2','xi4','xi6'])

#######################################################################

#nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales
nr = np.linspace(5.,150.,30)[:-1] #Scales

fig, axes = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True, figsize=(18, 8))

for ax,Nran in zip(axes[0],[2*87**3,4*87**3,8*87**3,16*87**3]):

    #-----------------------------------------------------------------------
    #PLOT 1 - xi_0 Monopole
    #-----------------------------------------------------------------------
    x = nr
    if Nran==2*87**3: 
        yran = ran_1['xi0']-xi_noran['xi0']
        yrc  = rc_1['xi0']-xi_noran['xi0']
        yrs  = rs_1['xi0']-xi_noran['xi0']
        yzr  = zr_1['xi0']-xi_noran['xi0']
    if Nran==4*87**3: 
        yran = ran_2['xi0']-xi_noran['xi0']
        yrc  = rc_2['xi0']-xi_noran['xi0']
        yrs  = rs_2['xi0']-xi_noran['xi0']
        yzr  = zr_2['xi0']-xi_noran['xi0']
    if Nran==8*87**3: 
        yran = ran_4['xi0']-xi_noran['xi0']
        yrc  = rc_4['xi0']-xi_noran['xi0']
        yrs  = rs_4['xi0']-xi_noran['xi0']
        yzr  = zr_4['xi0']-xi_noran['xi0']
    if Nran==16*87**3: 
        yran = ran_8['xi0']-xi_noran['xi0']
        yrc  = rc_8['xi0']-xi_noran['xi0']
        yrs  = rs_8['xi0']-xi_noran['xi0']
        yzr  = zr_8['xi0']-xi_noran['xi0']


    ax.plot(x,yran,color='k', linestyle=':',linewidth=2, label='Random')
    ax.plot(x,yrc,color=green,linestyle='-.',label='Crossed Random')
    ax.plot(x,yrs,color=blue,linestyle='--',label='Split Random')
    ax.plot(x,yzr,color=red,linestyle='-',label='Glass')

    #if Nran==2*87**3: ax1.set_ylabel(r'$\Delta \bar{\xi_0(r)}$',fontsize=fs)


    #if Nran==2*87**3: ax1.legend(fontsize=lfs,loc='lower right')
    #ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    #if Nran != 2*87**3: ax1.set_yticklabels([])

    #ax.set_ylim(-.0002,.0002)
    ax.tick_params(axis='both',labelsize=fs-6)
    ax.set_yticks([-2E-4, 0., 2E-4])

for ax, Nran in zip(axes[1],[2*87**3,4*87**3,8*87**3,16*87**3]):

    #-----------------------------------------------------------------------
    #PLOT 2 - xi_2 Dipole
    #-----------------------------------------------------------------------
    if Nran==2*87**3: 
        yran = ran_1['xi2']-xi_noran['xi2']
        yrc  = rc_1['xi2']-xi_noran['xi2']
        yrs  = rs_1['xi2']-xi_noran['xi2']
        yzr  = zr_1['xi2']-xi_noran['xi2']
    if Nran==4*87**3: 
        yran = ran_2['xi2']-xi_noran['xi2']
        yrc  = rc_2['xi2']-xi_noran['xi2']
        yrs  = rs_2['xi2']-xi_noran['xi2']
        yzr  = zr_2['xi2']-xi_noran['xi2']
    if Nran==8*87**3: 
        yran = ran_4['xi2']-xi_noran['xi2']
        yrc  = rc_4['xi2']-xi_noran['xi2']
        yrs  = rs_4['xi2']-xi_noran['xi2']
        yzr  = zr_4['xi2']-xi_noran['xi2']
    if Nran==16*87**3: 
        yran = ran_8['xi2']-xi_noran['xi2']
        yrc  = rc_8['xi2']-xi_noran['xi2']
        yrs  = rs_8['xi2']-xi_noran['xi2']
        yzr  = zr_8['xi2']-xi_noran['xi2']


    ax.plot(x,yran,color='k', linestyle=':',linewidth=2)
    ax.plot(x,yrc,color=green,linestyle='-.')
    ax.plot(x,yrs,color=blue,linestyle='--')
    ax.plot(x,yzr,color=red,linestyle='-')
    #if Nran==2*87**3: ax2.set_ylabel(r'$\Delta \bar{\xi_2(r)}$',fontsize=fs)

    #ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    #if Nran != 2*87**3: ax2.set_yticklabels([])
    #ax.set_ylim(-.0004,.0007)
    ax.tick_params(axis='both',labelsize=fs-6)
    ax.set_yticks([-2E-4, 0., 2E-4])


for ax,Nran in zip(axes[2],[2*87**3,4*87**3,8*87**3,16*87**3]):

    #-----------------------------------------------------------------------
    #PLOT 3 - xi_4 Hexadecapole
    #-----------------------------------------------------------------------
    if Nran==2*87**3: 
        yran = ran_1['xi4']-xi_noran['xi4']
        yrc  = rc_1['xi4']-xi_noran['xi4']
        yrs  = rs_1['xi4']-xi_noran['xi4']
        yzr  = zr_1['xi4']-xi_noran['xi4']
    if Nran==4*87**3: 
        yran = ran_2['xi4']-xi_noran['xi4']
        yrc  = rc_2['xi4']-xi_noran['xi4']
        yrs  = rs_2['xi4']-xi_noran['xi4']
        yzr  = zr_2['xi4']-xi_noran['xi4']
    if Nran==8*87**3: 
        yran = ran_4['xi4']-xi_noran['xi4']
        yrc  = rc_4['xi4']-xi_noran['xi4']
        yrs  = rs_4['xi4']-xi_noran['xi4']
        yzr  = zr_4['xi4']-xi_noran['xi4']
    if Nran==16*87**3: 
        yran = ran_8['xi4']-xi_noran['xi4']
        yrc  = rc_8['xi4']-xi_noran['xi4']
        yrs  = rs_8['xi4']-xi_noran['xi4']
        yzr  = zr_8['xi4']-xi_noran['xi4']

    ax.plot(x,yran,color='k', linestyle=':',linewidth=2)
    ax.plot(x,yrc,color=green,linestyle='-.')
    ax.plot(x,yrs,color=blue,linestyle='--')
    ax.plot(x,yzr,color=red,linestyle='-')
    #if Nran==2*87**3: ax3.set_ylabel(r'$\Delta \bar{\xi_4(r)}$',fontsize=fs)

    #ax3.set_xlabel('$R\,[Mpc]$',fontsize=fs)
    ax.set_ylim(-.0003,.0003)
    ax.set_yticks([-2E-4, 0., 2E-4])

    ax.tick_params(axis='both',labelsize=fs-6)

    #if Nran != 2*87**3: ax3.set_yticklabels([])

axes[0,0].set_ylabel(r'$\Delta \xi_0(s)$',fontsize=fs)
axes[1,0].set_ylabel(r'$\Delta \xi_2(s)$',fontsize=fs)
axes[2,0].set_ylabel(r'$\Delta \xi_4(s)$',fontsize=fs)

axes[2,0].set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs)
axes[2,1].set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs)
axes[2,2].set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs)
axes[2,3].set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs)

axes[0,0].set_title(r'$\alpha=2$',fontsize=fs)
axes[0,1].set_title(r'$\alpha=4$',fontsize=fs)
axes[0,2].set_title(r'$\alpha=8$',fontsize=fs)
axes[0,3].set_title(r'$\alpha=16$',fontsize=fs)

#axes[0,0].legend(markerscale=1,ncol=2,loc='lower left')
#f.suptitle(r'$R = {}Mpc$'.format(nr[ir]),fontsize=15)
fig.subplots_adjust(hspace=0)
fig.tight_layout()
fig.savefig('../plots/mean/mean_Niter{}.png'.format(Niter),dpi=fig.dpi)
fig.savefig('../plots/mean/mean_Niter{}.pdf'.format(Niter),dpi=fig.dpi)
#fig.savefig('../plots/mean/mean.png',dpi=fig.dpi)
#fig.savefig('../plots/mean/mean.pdf',dpi=fig.dpi)

#plt.show()
plt.close()

