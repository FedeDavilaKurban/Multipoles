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

Niter=100

#######################################################################
#Reading
#######################################################################
#Ran
ran_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(2*87**3))
ran_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(4*87**3))
ran_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(8*87**3))
ran_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(16*87**3))
#ran_16 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}_smallscale.txt'.format(16*87**3))

#RanCross
rc_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}_smallscale.txt'.format(2*87**3))
rc_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}_smallscale.txt'.format(4*87**3))
rc_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}_smallscale.txt'.format(8*87**3))
rc_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}_smallscale.txt'.format(16*87**3))

#RanSplit
rs_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}_smallscale.txt'.format(2*87**3))
rs_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}_smallscale.txt'.format(4*87**3))
rs_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}_smallscale.txt'.format(8*87**3))
rs_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}_smallscale.txt'.format(16*87**3))

#ZelRec
#zr_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(2*87**3))
#zr_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(4*87**3))
#zr_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(8*87**3))
#zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale_Niter{}.txt'.format(16*87**3,Niter))
#zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}_smallscale.txt'.format(16*87**3))

zr_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/adaptativeNbins/sd_zelrec_{}_smallscale_Niter{}.txt'.format(2*87**3,Niter))
zr_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/adaptativeNbins/sd_zelrec_{}_smallscale_Niter{}.txt'.format(4*87**3,Niter))
zr_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/adaptativeNbins/sd_zelrec_{}_smallscale_Niter{}.txt'.format(8*87**3,Niter))
zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/adaptativeNbins/sd_zelrec_{}_smallscale_Niter{}.txt'.format(16*87**3,Niter))

#######################################################################

nr = np.geomspace(0.5,40.,15)[:-1] #Scales
#nr = np.geomspace(0.5,90.,15)[:-1]  / 0.695 #Scales

# #######################################################################
# # STD_i/STD_Ran vs NRAN
# #######################################################################
# x = [2*87**3,4*87**3,8*87**3,16*87**3]

# fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(6, 10))
# for ax, xil in zip(axes.flatten(),['xi0','xi2','xi4']):
#     rc  = [rc_1[xil]/ran_1[xil], rc_2[xil]/ran_2[xil],  rc_4[xil]/ran_4[xil],  rc_8[xil]/ran_8[xil]]
#     rs  = [rs_1[xil]/ran_1[xil], rs_2[xil]/ran_2[xil],  rs_4[xil]/ran_4[xil],  rs_8[xil]/ran_8[xil]]
#     zr  = [zr_1[xil]/ran_1[xil], zr_2[xil]/ran_2[xil],  zr_4[xil]/ran_4[xil],  zr_8[xil]/ran_8[xil]]
    
#     #yrc = []
#     #yrs = []
#     for ir in [1,10]:

#         yrc = [rc[0][ir],rc[1][ir],rc[2][ir],rc[3][ir]] 
#         yrs = [rs[0][ir],rs[1][ir],rs[2][ir],rs[3][ir]] 
#         yzr = [zr[0][ir],zr[1][ir],zr[2][ir],zr[3][ir]]

#         if ir==1: 
#             ax.plot(x,yrc,color=green,linestyle='-',label='Crossed Random')
#             ax.plot(x,yrs,color=blue,linestyle='-',label='Split Random')
#         if ir!=1: 
#             ax.plot(x,yrc,color=green,linestyle='-')
#             ax.plot(x,yrs,color=blue,linestyle='-')

#         if ir==1: ax.plot(x,yzr,color=red,marker='x',linestyle='-',label='Zel. Rec.; '+r'$r=14.4Mpc$')
#         if ir==10: ax.plot(x,yzr,color=red,marker='*',linestyle='-',label='Zel. Rec.; '+r'$r=79.1Mpc$')
#         if ir==19: ax.plot(x,yzr,color=red,marker='^',linestyle='-',label='Zel. Rec.; '+r'$r=143.9Mpc$')
#         if ir==28: ax.plot(x,yzr,color=red,marker='o',linestyle='-',label='Zel. Rec.; '+r'$r=208.6Mpc$')


#         #Plot Parameters
#         ax.get_yaxis().get_major_formatter().labelOnlyBase = False
#         ax.set_ylim(2E-1,15E-1)
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#         plt.setp(ax.get_xticklabels('both'), visible=False)
	    
#         if xil== 'xi0': ax.set_ylabel(r'$\sigma_{\xi_0(r)}$',fontsize=fs)
#         if xil== 'xi2': ax.set_ylabel(r'$\sigma_{\xi_2(r)}$',fontsize=fs)
#         if xil== 'xi4': 
#             ax.set_ylabel(r'$\sigma_{\xi_4(r)}$',fontsize=fs)
#             ax.set_xlabel(r'$N*\delta_{Minerva}$',fontsize=fs-4)

#         ax.tick_params(axis='both',which='minor',bottom=False,top=False,left=True,right=True,labelleft=False,labelbottom=False)
#         ax.set_xticks(x)
#         ax.set_xticklabels(['1','2','4','8'])
#         ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#         ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#         ax.tick_params(axis='both',labelsize=fs-6)

#         ax.set_yticks([.2, .4, .6, 1.])


# #first_legend = plt.legend(handles=[rcplot,rsplot], loc='upper right')
# ## Add the legend manually to the current Axes.
# #ax = plt.gca().add_artist(first_legend)
# ## Create another legend for the second line.
# #plt.legend(handles=[line2], loc='lower right')

# plt.subplots_adjust(hspace=0)
# plt.tight_layout()
# axes[0].legend(markerscale=1)
# plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_nran_norm/stdnorm_ss.png',dpi=fig.dpi)
# plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_nran_norm/pdfs/stdnorm_ss.pdf',dpi=fig.dpi)


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
    #rs  = [rs_1[xil]/ran_1[xil], rs_2[xil]/ran_2[xil],  rs_4[xil]/ran_4[xil],  rs_8[xil]/ran_8[xil]]
    zr  = [zr_1[xil]/ran_1[xil], zr_2[xil]/ran_2[xil],  zr_4[xil]/ran_4[xil],  zr_8[xil]/ran_8[xil]]
    
    ax.plot(x,rc[0],color=blue,linestyle=':',label='Cross. Ran.; '+r'$N_r$')
    ax.plot(x,rc[1],color=blue,linestyle='--',label='Cross. Ran.; '+r'$2N_r$')
    ax.plot(x,rc[2],color=blue,linestyle='-.',label='Cross. Ran.; '+r'$4N_r$')
    ax.plot(x,rc[3],color=blue,linestyle='-',label='Cross. Ran.; '+r'$8N_r$')

    #ax.plot(x,rs[0],color=blue,linestyle='-',label='Split Ran.')
    #ax.plot(x,rs[1],color=blue,linestyle='-')
    #ax.plot(x,rs[2],color=blue,linestyle='-')
    #ax.plot(x,rs[3],color=blue,linestyle='-')

    ax.plot(x,zr[0],color=red,linestyle=':',label='Glass; '+r'$N_r$')
    ax.plot(x,zr[1],color=red,linestyle='--',label='Glass; '+r'$2N_r$')
    ax.plot(x,zr[2],color=red,linestyle='-.',label='Glass; '+r'$4N_r$')
    ax.plot(x,zr[3],color=red,linestyle='-',label='Glass; '+r'8$N_r$')

    if xil== 'xi0': ax.set_ylabel(r'$\sigma_{\xi_0}/\sigma_{\xi_{0_\mathrm{LS}}}$',fontsize=fs+6)
    if xil== 'xi2': ax.set_ylabel(r'$\sigma_{\xi_2}/\sigma_{\xi_{2_\mathrm{LS}}}$',fontsize=fs+6)
    if xil== 'xi4': 
        ax.set_ylabel(r'$\sigma_{\xi_4}/\sigma_{\xi_{4_\mathrm{LS}}}$',fontsize=fs+6)
        ax.set_xlabel(r'$s\; [Mpc\,h^{-1}]$',fontsize=fs-2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(.1,2.5)
    ax.tick_params(axis='y',which='minor',bottom=False,top=False,left=True,right=True,labelleft=False,labelbottom=False)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([.2, .4, .6, 1.])
    ax.tick_params(axis='both',labelsize=fs-6)


plt.tight_layout()
#axes[0].legend(markerscale=1,ncol=2)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/snorm_vs_r_ss_Niter{}.png'.format(Niter),dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/pdfs/snorm_vs_r_ss_Niter{}.pdf'.format(Niter),dpi=fig.dpi)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/snorm_vs_r_ss_adaptativeNbins.png'.format(Niter),dpi=fig.dpi)
plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/pdfs/snorm_vs_r_ss_adaptativeNbins.pdf'.format(Niter),dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/snorm_vs_r_ss.png',dpi=fig.dpi)
#plt.savefig('/home/fede/Proyectos/Multipoles/plots/smallscale/sd_vs_r_norm/pdfs/snorm_vs_r_ss.pdf',dpi=fig.dpi)

