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

#######################################################################
#Reading
#######################################################################
#Ran
#ran_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(87**3))
ran_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(2*87**3))
ran_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(4*87**3))
ran_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(8*87**3))
ran_16 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ran_{}.txt'.format(16*87**3))

#RanCross
rc_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(87**3))
rc_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(2*87**3))
rc_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(4*87**3))
rc_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_rancross_{}.txt'.format(8*87**3))

#RanSplit
rs_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(87**3))
rs_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(2*87**3))
rs_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(4*87**3))
rs_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_ransplit_{}.txt'.format(8*87**3))

#Glass
g_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_glass_{}.txt'.format(87**3))
g_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_glass_{}.txt'.format(2*87**3))
g_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_glass_{}.txt'.format(4*87**3))
g_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_glass_{}.txt'.format(8*87**3))

#ZelRec
zr_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(87**3))
zr_2 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(2*87**3))
zr_4 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(4*87**3))
zr_8 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/sd_zelrec_{}.txt'.format(8*87**3))

#CCVT
#c_1 = ascii.read('/home/fede/Proyectos/Multipoles/data/out/sd_xi/std_c_{}.txt'.format(87**3))
#######################################################################

nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales



#######################################################################
# STD_i/STD_Ran vs NRAN
#######################################################################
for ir in range(len(nr)):
#for ir in range(1):

	#-----------------------------------------------------------------------
	#PLOT 1 - xi_0 Monopole
	#-----------------------------------------------------------------------
	x = [87**3,2*87**3,4*87**3,8*87**3]
	#yran = [ran_2['xi0'][ir], ran_4['xi0'][ir], ran_8['xi0'][ir], ran_16['xi0'][ir]]
	yrc  = [rc_1['xi0'][ir]/ran_2['xi0'][ir], rc_2['xi0'][ir]/ran_4['xi0'][ir],  rc_4['xi0'][ir]/ran_8['xi0'][ir],  rc_8['xi0'][ir]/ran_16['xi0'][ir]]
	yrs  = [rs_1['xi0'][ir]/ran_2['xi0'][ir], rs_2['xi0'][ir]/ran_4['xi0'][ir],  rs_4['xi0'][ir]/ran_8['xi0'][ir],  rs_8['xi0'][ir]/ran_16['xi0'][ir]]
	yg   = [g_1['xi0'][ir]/ran_2['xi0'][ir] , g_2['xi0'][ir]/ran_4['xi0'][ir],   g_4['xi0'][ir]/ran_8['xi0'][ir],   g_8['xi0'][ir]/ran_16['xi0'][ir]]
	yzr  = [zr_1['xi0'][ir]/ran_2['xi0'][ir], zr_2['xi0'][ir]/ran_4['xi0'][ir],  zr_4['xi0'][ir]/ran_8['xi0'][ir],  zr_8['xi0'][ir]/ran_16['xi0'][ir]]
	#yc   = [c_1['xi0'][ir]/ran_2['xi0'][ir]]

	f = plt.figure(figsize=(6,10))
	ax1 = f.add_subplot(311)#, sharex=ax2)
	#ax1.plot(x,yran,color=blue,marker='^',linestyle='--',label='Random')
	ax1.plot(x,yrc,color=green,linestyle='-.',label='Crossed Random')
	ax1.plot(x,yrs,color='k',linestyle=':',label='Split Random')
	#ax1.plot(87**3,yc,color=other,marker='s',linestyle='-.',label='CCVT')
	ax1.plot(x,yg,color=blue,linestyle='--',label='Glass')
	ax1.plot(x,yzr,color=red,marker='o',linestyle='-',label='Zel. Rec.')
	ax1.set_ylabel(r'$\sigma_{\xi_0(r)}$',fontsize=fs)
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height*.8])
	#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),ncol=3, fancybox=True, shadow=True,fontsize=lfs)
	ax1.legend()
	ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	ax1.set_ylim(2E-1,14E-1)

	#-----------------------------------------------------------------------
	#PLOT 2 - xi_2 Dipole
	#-----------------------------------------------------------------------
	#yran = [ran_1['xi2'][ir], ran_2['xi2'][ir], ran_4['xi2'][ir], ran_8['xi2'][ir]]
	yrc  = [rc_1['xi2'][ir]/ran_2['xi2'][ir], rc_2['xi2'][ir]/ran_4['xi2'][ir], rc_4['xi2'][ir]/ran_8['xi2'][ir], rc_8['xi2'][ir]/ran_16['xi2'][ir]]
	yrs  = [rs_1['xi2'][ir]/ran_2['xi2'][ir], rs_2['xi2'][ir]/ran_4['xi2'][ir], rs_4['xi2'][ir]/ran_8['xi2'][ir], rs_8['xi2'][ir]/ran_16['xi2'][ir]]
	yg   = [g_1['xi2'][ir]/ran_2['xi2'][ir],  g_2['xi2'][ir]/ran_4['xi2'][ir],  g_4['xi2'][ir]/ran_8['xi2'][ir],  g_8['xi2'][ir]/ran_16['xi2'][ir]]
	yzr  = [zr_1['xi2'][ir]/ran_2['xi2'][ir], zr_2['xi2'][ir]/ran_4['xi2'][ir], zr_4['xi2'][ir]/ran_8['xi2'][ir], zr_8['xi2'][ir]/ran_16['xi2'][ir]]
	#yc   = [c_1['xi2'][ir]/ran_2['xi2'][ir]]


	ax2 = f.add_subplot(312,sharex=ax1)
	#ax2.plot(x,yran,color=blue,marker='^',linestyle='--')
	ax2.plot(x,yrc,color=green,linestyle='-.')
	#ax2.plot(87**3,yc,color=other,marker='s',linestyle='-.')
	ax2.plot(x,yrs,color='k',linestyle=':')
	ax2.plot(x,yg,color=blue,linestyle='--')
	ax2.plot(x,yzr,color=red,marker='o',linestyle='-')
	ax2.set_ylabel(r'$\sigma_{\xi_2(r)}$',fontsize=fs)
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	ax2.set_ylim(3E-1,15E-1)

	#-----------------------------------------------------------------------
	#PLOT 3 - xi_4 Hexadecapole
	#-----------------------------------------------------------------------
	#yran = [ran_1['xi4'][ir],ran_2['xi4'][ir],ran_4['xi4'][ir],ran_8['xi4'][ir]]
	yrc  = [rc_1['xi4'][ir]/ran_2['xi4'][ir], rc_2['xi4'][ir]/ran_4['xi4'][ir], rc_4['xi4'][ir]/ran_8['xi4'][ir], rc_8['xi4'][ir]/ran_16['xi4'][ir]]
	yrs  = [rs_1['xi4'][ir]/ran_2['xi4'][ir], rs_2['xi4'][ir]/ran_4['xi4'][ir], rs_4['xi4'][ir]/ran_8['xi4'][ir], rs_8['xi4'][ir]/ran_16['xi4'][ir]]
	yg   = [g_1['xi4'][ir]/ran_2['xi4'][ir],  g_2['xi4'][ir]/ran_4['xi4'][ir],  g_4['xi4'][ir]/ran_8['xi4'][ir],  g_8['xi4'][ir]/ran_16['xi4'][ir]]
	yzr  = [zr_1['xi4'][ir]/ran_2['xi4'][ir], zr_2['xi4'][ir]/ran_4['xi4'][ir], zr_4['xi4'][ir]/ran_8['xi4'][ir], zr_8['xi4'][ir]/ran_16['xi4'][ir]]
	#yc   = [c_1['xi4'][ir]]

	ax3 = f.add_subplot(313, sharex=ax2)
	#ax3.plot(x,yran,color=blue,marker='^',linestyle='--')
	ax3.plot(x,yrc,color=green,linestyle='-.')
	ax3.plot(x,yrs,color='k',linestyle=':')
	#ax3.plot(87**3,yc,color=other,marker='s',linestyle='-.')
	ax3.plot(x,yg,color=blue,linestyle='--')
	ax3.plot(x,yzr,color=red,marker='o',linestyle='-')
	ax3.set_ylabel(r'$\sigma_{\xi_4(r)}$',fontsize=fs)
	ax3.set_xscale('log')
	ax3.set_yscale('log')
	ax3.set_xlabel('$N_{ran}$',fontsize=fs)
	ax3.set_ylim(1E-1,1.3E0)

	f.suptitle(r'$R = {}Mpc$'.format(nr[ir]),fontsize=15)
	f.subplots_adjust(hspace=0)
	f.tight_layout()
	f.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_nran_norm/stdnorm_{}.png'.format(ir),dpi=f.dpi)
	#plt.show()
	plt.close()




#######################################################################
# STD_i/STD_Ran vs R
#######################################################################
#Nran = 87**3

for Nran in [87**3,2*87**3,4*87**3,8*87**3]:

	#-----------------------------------------------------------------------
	#PLOT 1 - xi_0 Monopole
	#-----------------------------------------------------------------------
	x = nr
	if Nran==87**3: 
		yran = ran_2['xi0']
		yrc  = rc_1['xi0']
		yrs  = rs_1['xi0']
		yg   = g_1['xi0']
		yzr  = zr_1['xi0']
		#yc = c_1['xi0']
	if Nran==2*87**3: 
		yran = ran_4['xi0']
		yrc  = rc_2['xi0']
		yrs  = rs_2['xi0']
		yg   = g_2['xi0']
		yzr  = zr_2['xi0']
	if Nran==4*87**3: 
		yran = ran_8['xi0']
		yrc  = rc_4['xi0']
		yrs  = rs_4['xi0']
		yg   = g_4['xi0']
		yzr  = zr_4['xi0']
	if Nran==8*87**3: 
		yran = ran_16['xi0']
		yrc  = rc_8['xi0']
		yrs  = rs_8['xi0']
		yg   = g_8['xi0']
		yzr  = zr_8['xi0']


	f = plt.figure(figsize=(6,6))
	ax1 = f.add_subplot(311)#, sharex=ax2)
	#ax1.plot(x,yran,color=blue,marker='^',linestyle='--',label='Random')
	ax1.plot(x,yrc/yran,color=green,linestyle='-.',label='Crossed Random')
	ax1.plot(x,yrs/yran,color='k',linestyle=':',label='Split Random')
	#if Nran==87**3: ax1.plot(x,yc/yran,color=other,marker='s',linestyle='-.',label='CCVT')
	ax1.plot(x,yg/yran,color=blue,linestyle='--',label='Glass')
	ax1.plot(x,yzr/yran,color=red,marker='o',linestyle='-',label='Zel. Rec.')
	ax1.set_ylabel(r'$\sigma_{\xi_0(r)}$',fontsize=fs)
	#ax1.set_xscale('log')
	ax1.set_yscale('log')
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height*.8])
	#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),ncol=3, fancybox=True, shadow=True,fontsize=lfs)
	ax1.legend()
	ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	ax1.set_ylim(2E-1,13E-1)

	#-----------------------------------------------------------------------
	#PLOT 2 - xi_2 Dipole
	#-----------------------------------------------------------------------
	if Nran==87**3: 
		yran = ran_2['xi2']
		yrc  = rc_1['xi2']
		yrs  = rs_1['xi2']
		yg   = g_1['xi2']
		yzr  = zr_1['xi2']
		#yc = c_1['xi2']
	if Nran==2*87**3: 
		yran = ran_4['xi2']
		yrc  = rc_2['xi2']
		yrs  = rs_2['xi2']
		yg   = g_2['xi2']
		yzr  = zr_2['xi2']
	if Nran==4*87**3: 
		yran = ran_8['xi2']
		yrc  = rc_4['xi2']
		yrs  = rs_4['xi2']
		yg   = g_4['xi2']
		yzr  = zr_4['xi2']
	if Nran==8*87**3: 
		yran = ran_16['xi2']
		yrc  = rc_8['xi2']
		yrs  = rs_8['xi2']
		yg   = g_8['xi2']
		yzr  = zr_8['xi2']


	ax2 = f.add_subplot(312,sharex=ax1)
	#ax2.plot(x,yran,color=blue,marker='^',linestyle='--')
	ax2.plot(x,yrc/yran,color=green,linestyle='-.')
	#if Nran==87**3: ax2.plot(x,yc/yran,color=other,marker='s',linestyle='-.')
	ax2.plot(x,yrs/yran,color='k',linestyle=':')
	ax2.plot(x,yg/yran,color=blue,linestyle='--')
	ax2.plot(x,yzr/yran,color=red,marker='o',linestyle='-')
	ax2.set_ylabel(r'$\sigma_{\xi_2(r)}$',fontsize=fs)
	#ax2.set_xscale('log')
	ax2.set_yscale('log')
	ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	ax2.set_ylim(2E-1,13E-1)

	#-----------------------------------------------------------------------
	#PLOT 3 - xi_4 Hexadecapole
	#-----------------------------------------------------------------------
	if Nran==87**3: 
		yran = ran_2['xi4']
		yrc  = rc_1['xi4']
		yrs  = rs_1['xi4']
		yg   = g_1['xi4']
		yzr  = zr_1['xi4']
		#yc = c_1['xi4']
	if Nran==2*87**3: 
		yran = ran_4['xi4']
		yrc  = rc_2['xi4']
		yrs  = rs_2['xi4']
		yg   = g_2['xi4']
		yzr  = zr_2['xi4']
	if Nran==4*87**3: 
		yran = ran_8['xi4']
		yrc  = rc_4['xi4']
		yrs  = rs_4['xi4']
		yg   = g_4['xi4']
		yzr  = zr_4['xi4']
	if Nran==8*87**3: 
		yran = ran_16['xi4']
		yrc  = rc_8['xi4']
		yrs  = rs_8['xi4']
		yg   = g_8['xi4']
		yzr  = zr_8['xi4']

	ax3 = f.add_subplot(313, sharex=ax2)
	#ax3.plot(x,yran,color=blue,marker='^',linestyle='--')
	ax3.plot(x,yrc/yran,color=green,linestyle='-.')
	ax3.plot(x,yrs/yran,color='k',linestyle=':')
	#if Nran==87**3: ax3.plot(x,yc/yran,color=other,marker='s',linestyle='-.')
	ax3.plot(x,yg/yran,color=blue,linestyle='--')
	ax3.plot(x,yzr/yran,color=red,marker='o',linestyle='-')
	ax3.set_ylabel(r'$\sigma_{\xi_4(r)}$',fontsize=fs)
	#ax3.set_xscale('log')
	ax3.set_yscale('log')
	ax3.set_xlabel('$r\, [Mpc]$',fontsize=fs)
	ax3.set_ylim(2E-1,13E-1)

	#f.suptitle(r'$R = {}Mpc$'.format(nr[ir]),fontsize=15)
	f.subplots_adjust(hspace=0)
	f.tight_layout()
	f.savefig('/home/fede/Proyectos/Multipoles/plots/sd_vs_r_norm/snorm_vs_r_{}.png'.format(Nran),dpi=f.dpi)
	#plt.show()
	plt.close()
	