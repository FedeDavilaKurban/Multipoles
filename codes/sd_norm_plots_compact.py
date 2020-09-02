"""
Quiero compactar las 4 columnas en una
"""
import glob 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table
import matplotlib.ticker as mticker

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
ran_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(2*87**3))
ran_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(4*87**3))
ran_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(8*87**3))
ran_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(16*87**3))
#ran_16 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(16*87**3))

#RanCross
rc_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(2*87**3))
rc_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(4*87**3))
rc_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(8*87**3))
rc_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(16*87**3))

#RanSplit
rs_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(2*87**3))
rs_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(4*87**3))
rs_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(8*87**3))
rs_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(16*87**3))

#ZelRec
zr_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(2*87**3))
zr_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(4*87**3))
zr_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(8*87**3))
zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_Niter{}.txt'.format(16*87**3,Niter))
#zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(16*87**3))

#######################################################################

#nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales
nr = np.linspace(5.,150.,30)[:-1] #Scales

#nr = np.linspace(5.,150.,30)[:-1]  / 0.695 #Scales

#######################################################################
# STD_i vs NRAN
#######################################################################
x = [2*87**3,4*87**3,8*87**3,16*87**3]

fig = plt.figure(figsize=(10, 6))
#for xil in ['xi0','xi2','xi4']:
#rs  = [rs_1[xil]/ran_1[xil], rs_2[xil]/ran_2[xil],  rs_4[xil]/ran_4[xil],  rs_8[xil]/ran_8[xil]]

# rc0  = [rc_1['xi0']/ran_1['xi0'], rc_2['xi0']/ran_2['xi0'],  rc_4['xi0']/ran_4['xi0'],  rc_8['xi0']/ran_8['xi0']]
# rc2  = [rc_1['xi2']/ran_1['xi2'], rc_2['xi2']/ran_2['xi2'],  rc_4['xi2']/ran_4['xi2'],  rc_8['xi2']/ran_8['xi2']]
# rc4  = [rc_1['xi4']/ran_1['xi4'], rc_2['xi4']/ran_2['xi4'],  rc_4['xi4']/ran_4['xi4'],  rc_8['xi4']/ran_8['xi4']]

# zr0  = [zr_1['xi0']/ran_1['xi0'], zr_2['xi0']/ran_2['xi0'],  zr_4['xi0']/ran_4['xi0'],  zr_8['xi0']/ran_8['xi0']]
# zr2  = [zr_1['xi2']/ran_1['xi2'], zr_2['xi2']/ran_2['xi2'],  zr_4['xi2']/ran_4['xi2'],  zr_8['xi2']/ran_8['xi2']]
# zr4  = [zr_1['xi4']/ran_1['xi4'], zr_2['xi4']/ran_2['xi4'],  zr_4['xi4']/ran_4['xi4'],  zr_8['xi4']/ran_8['xi4']]

rc0  = [rc_1['xi0'], rc_2['xi0'],  rc_4['xi0'],  rc_8['xi0']]
rc2  = [rc_1['xi2'], rc_2['xi2'],  rc_4['xi2'],  rc_8['xi2']]
rc4  = [rc_1['xi4'], rc_2['xi4'],  rc_4['xi4'],  rc_8['xi4']]

zr0  = [zr_1['xi0'], zr_2['xi0'],  zr_4['xi0'],  zr_8['xi0']]
zr2  = [zr_1['xi2'], zr_2['xi2'],  zr_4['xi2'],  zr_8['xi2']]
zr4  = [zr_1['xi4'], zr_2['xi4'],  zr_4['xi4'],  zr_8['xi4']]

#yrc = []
#yrs = []
for ir in [21]:

    yrc0 = [rc0[0][ir],rc0[1][ir],rc0[2][ir],rc0[3][ir]] 
    yzr0 = [zr0[0][ir],zr0[1][ir],zr0[2][ir],zr0[3][ir]]

    yrc2 = [rc2[0][ir],rc2[1][ir],rc2[2][ir],rc2[3][ir]] 
    yzr2 = [zr2[0][ir],zr2[1][ir],zr2[2][ir],zr2[3][ir]]

    yrc4 = [rc4[0][ir],rc4[1][ir],rc4[2][ir],rc4[3][ir]] 
    yzr4 = [zr4[0][ir],zr4[1][ir],zr4[2][ir],zr4[3][ir]]

    ms=8

    plt.plot(x,yrc0,color=blue,ms=ms,marker='o',linestyle='-',label=r'$l=0$; C.R.')
    plt.plot(x,yrc2,color=blue,ms=ms,marker='*',linestyle='--',label=r'$l=2$; C.R.')
    plt.plot(x,yrc4,color=blue,ms=ms,marker='^',linestyle='-.',label=r'$l=4$; C.R.')

    plt.plot(x,yzr0,color=red,ms=ms,marker='o',linestyle='-',label=r'$l=0$; G.')
    plt.plot(x,yzr2,color=red,ms=ms,marker='*',linestyle='--',label=r'$l=2$; G.')
    plt.plot(x,yzr4,color=red,ms=ms,marker='^',linestyle='-.',label=r'$l=4$; G.')

    
    #Plot Parameters
    #ax.get_yaxis().get_major_formatter().labelOnlyBase = False
    plt.ylim(1E-5,1E-3)
    plt.xscale('log')
    plt.yscale('log')
    #plt.setp(plt.xticklabels('both'), visible=False)
    
    #if xil== 'xi0': ax.set_ylabel(r'$\sigma_{\xi_0}/\sigma_{\xi_{0_\mathrm{LS}}}$',fontsize=fs)
    #if xil== 'xi2': ax.set_ylabel(r'$\sigma_{\xi_2}/\sigma_{\xi_{2_\mathrm{LS}}}$',fontsize=fs)
    #if xil== 'xi4': 
    #    ax.set_ylabel(r'$\sigma_{\xi_4}/\sigma_{\xi_{4_\mathrm{LS}}}$',fontsize=fs)
    #    ax.set_xlabel(r'$N_r/N_d$',fontsize=fs-4)

    plt.title(r's$\cong 110Mpc\,h^{-1}$',fontsize=fs)
    plt.tick_params(axis='both',which='minor',bottom=False,top=False,left=True,right=True,labelleft=False,labelbottom=False)
    plt.xticks(x,labels=['1','2','4','8'])
    #plt.yticks([.2, .4, .6, 1.], labels=['0.2','0.4','0.6','1.0'])
    plt.ylabel(r'$\sigma_{\xi_\ell}$',fontsize=fs+10)
    plt.xlabel(r'$\alpha$',fontsize=fs+6)

    #ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    plt.tick_params(axis='both',labelsize=fs)



#first_legend = plt.legend(handles=[rcplot,rsplot], loc='upper right')
## Add the legend manually to the current Axes.
#ax = plt.gca().add_artist(first_legend)
## Create another legend for the second line.
#plt.legend(handles=[line2], loc='lower right')

plt.subplots_adjust(hspace=0)
plt.tight_layout()
#plt.legend(markerscale=1.1,ncol=2,fontsize=fs-4)#,loc='upper left')
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_nran_norm/stdnorm.png',dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_nran_norm/pdfs/stdnorm.pdf',dpi=fig.dpi)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_nran_norm/sd_vs_n.png',dpi=fig.dpi)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_nran_norm/pdfs/sd_vs_n.pdf',dpi=fig.dpi)

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

#######################################################################
# STD_i/STD_Ran vs R
#######################################################################
x = nr

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(6, 10))
for ax, xil in zip(axes.flatten(),['xi0','xi2','xi4']):
    rc  = [rc_1[xil]/ran_1[xil], rc_2[xil]/ran_2[xil],  rc_4[xil]/ran_4[xil],  rc_8[xil]/ran_8[xil]]
    rs  = [rs_1[xil]/ran_1[xil], rs_2[xil]/ran_2[xil],  rs_4[xil]/ran_4[xil],  rs_8[xil]/ran_8[xil]]
    zr  = [zr_1[xil]/ran_1[xil], zr_2[xil]/ran_2[xil],  zr_4[xil]/ran_4[xil],  zr_8[xil]/ran_8[xil]]
    
    ax.plot(x,rc[0],color=blue,linestyle='-',label='C.R.; '+r'$N_r$')
    ax.plot(x,rc[1],color=blue,linestyle='--',label='C.R.; '+r'$2N_r$')
    ax.plot(x,rc[2],color=blue,linestyle='-.',label='C.R.; '+r'$4N_r$')
    ax.plot(x,rc[3],color=blue,linestyle='-',label='C.R.; '+r'$8N_r$')

    # ax.plot(x,rs[0],color=blue,linestyle='-',label='Split Ran.')
    # ax.plot(x,rs[1],color=blue,linestyle='-')
    # ax.plot(x,rs[2],color=blue,linestyle='-')
    # ax.plot(x,rs[3],color=blue,linestyle='-')

    ax.plot(x,zr[0],color=red,linestyle=':',label='G.; '+r'$N_r$')
    ax.plot(x,zr[1],color=red,linestyle='--',label='G.; '+r'$2N_r$')
    ax.plot(x,zr[2],color=red,linestyle='-.',label='G.; '+r'$4N_r$')
    ax.plot(x,zr[3],color=red,linestyle='-',label='G.; '+r'8$N_r$')

    if xil== 'xi0': ax.set_ylabel(r'$\sigma_{\xi_0}/\sigma_{\xi_{0_\mathrm{LS}}}$',fontsize=fs+6)
    if xil== 'xi2': ax.set_ylabel(r'$\sigma_{\xi_2}/\sigma_{\xi_{2_\mathrm{LS}}}$',fontsize=fs+6)
    if xil== 'xi4': 
        ax.set_ylabel(r'$\sigma_{\xi_4}/\sigma_{\xi_{4_\mathrm{LS}}}$',fontsize=fs+6)
        ax.set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs-2)

    ax.set_yscale('log')
    ax.set_ylim(.15,1.35)
    ax.tick_params(axis='y',which='minor',bottom=False,top=False,left=True,right=True,labelleft=False,labelbottom=False)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([.2, .4, .6, 1.])
    ax.tick_params(axis='both',labelsize=fs-6)

plt.tight_layout()
#axes[0].legend(markerscale=.5,ncol=4)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_r_norm/snorm_vs_r_Niter{}.png'.format(Niter),dpi=fig.dpi)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_r_norm/pdfs/snorm_vs_r_Niter{}.pdf'.format(Niter),dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_r_norm/snorm_vs_r.png',dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_r_norm/pdfs/snorm_vs_r.pdf',dpi=fig.dpi)

