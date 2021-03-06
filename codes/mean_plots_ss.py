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
ran_1 = ascii.read('../data/out/mean_xi/mean_ran_{}_smallscale.txt'.format(2*87**3)) 
ran_2 = ascii.read('../data/out/mean_xi/mean_ran_{}_smallscale.txt'.format(4*87**3)) 
ran_4 = ascii.read('../data/out/mean_xi/mean_ran_{}_smallscale.txt'.format(8*87**3)) 
ran_8 = ascii.read('../data/out/mean_xi/mean_ran_{}_smallscale.txt'.format(16*87**3)) 

#RanCross
rc_1 = ascii.read('../data/out/mean_xi/mean_rancross_{}_smallscale.txt'.format(2*87**3)) 
rc_2 = ascii.read('../data/out/mean_xi/mean_rancross_{}_smallscale.txt'.format(4*87**3)) 
rc_4 = ascii.read('../data/out/mean_xi/mean_rancross_{}_smallscale.txt'.format(8*87**3)) 
rc_8 = ascii.read('../data/out/mean_xi/mean_rancross_{}_smallscale.txt'.format(16*87**3)) 

#RanSplit
rs_1 = ascii.read('../data/out/mean_xi/mean_ransplit_{}_smallscale.txt'.format(2*87**3)) 
rs_2 = ascii.read('../data/out/mean_xi/mean_ransplit_{}_smallscale.txt'.format(4*87**3)) 
rs_4 = ascii.read('../data/out/mean_xi/mean_ransplit_{}_smallscale.txt'.format(8*87**3)) 
rs_8 = ascii.read('../data/out/mean_xi/mean_ransplit_{}_smallscale.txt'.format(16*87**3)) 

#ZelRec
#zr_1 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_smallscale.txt'.format(2*87**3)) 
#zr_2 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_smallscale.txt'.format(4*87**3)) 
#zr_4 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_smallscale.txt'.format(8*87**3)) 
#zr_8 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_smallscale.txt'.format(16*87**3)) 
#zr_8 = ascii.read('../data/out/mean_xi/mean_zelrec_{}_smallscale_Niter{}.txt'.format(16*87**3,Niter)) 

zr_1 = ascii.read('../data/out/mean_xi/adaptativeNbins/mean_zelrec_{}_smallscale_Niter{}.txt'.format(2*87**3,Niter)) 
zr_2 = ascii.read('../data/out/mean_xi/adaptativeNbins/mean_zelrec_{}_smallscale_Niter{}.txt'.format(4*87**3,Niter)) 
zr_4 = ascii.read('../data/out/mean_xi/adaptativeNbins/mean_zelrec_{}_smallscale_Niter{}.txt'.format(8*87**3,Niter)) 
zr_8 = ascii.read('../data/out/mean_xi/adaptativeNbins/mean_zelrec_{}_smallscale_Niter{}.txt'.format(16*87**3,Niter)) 

#Analytic
xi_noran = ascii.read('../data/xi_l_noran_smallscale.txt',names=['xi0','xi2','xi4','xi6']) 

#######################################################################

nr = np.geomspace(0.5,40.,15)[:-1] #Scales

x=nr

for Nran in [2*87**3,4*87**3,8*87**3,16*87**3]:
#for Nran in [87**3]:
	#-----------------------------------------------------------------------
	#PLOT 1 - xi_0 Monopole
	#-----------------------------------------------------------------------
	x = nr
	if Nran==2*87**3: 
		yran = ran_1['xi0']-xi_noran['xi0']
		yrc  = rc_1['xi0']-xi_noran['xi0']
		yrs  = rs_1['xi0']-xi_noran['xi0']
		#yg   = g_1['xi0']-xi_noran['xi0']
		yzr  = zr_1['xi0']-xi_noran['xi0']
		#yc = c_1['xi0']-xi_noran['xi0']

		print(yran)
		print(yrs)
		print(xi_noran['xi0'])
	if Nran==4*87**3: 
		yran = ran_2['xi0']-xi_noran['xi0']
		yrc  = rc_2['xi0']-xi_noran['xi0']
		yrs  = rs_2['xi0']-xi_noran['xi0']
		#yg   = g_2['xi0']-xi_noran['xi0']
		yzr  = zr_2['xi0']-xi_noran['xi0']
	if Nran==8*87**3: 
		yran = ran_4['xi0']-xi_noran['xi0']
		yrc  = rc_4['xi0']-xi_noran['xi0']
		yrs  = rs_4['xi0']-xi_noran['xi0']
		#yg   = g_4['xi0']-xi_noran['xi0']
		yzr  = zr_4['xi0']-xi_noran['xi0']
	if Nran==16*87**3: 
		yran = ran_8['xi0']-xi_noran['xi0']
		yrc  = rc_8['xi0']-xi_noran['xi0']
		yrs  = rs_8['xi0']-xi_noran['xi0']
		#yg   = g_8['xi0']-xi_noran['xi0']
		yzr  = zr_8['xi0']-xi_noran['xi0']

	
	
	f = plt.figure(figsize=(6,10))
	ax1 = f.add_subplot(311)#, sharex=ax2)
	ax1.plot(x,yran,color='k',ls=':',linewidth=2.,label='Random')
	ax1.plot(x,yrc,color=green,linestyle='-.',label='Crossed Random')
	ax1.plot(x,yrs,color=blue,linestyle='--',label='Split Random')
	ax1.plot(x,yzr,color=red,linestyle='-',label='Zel. Rec.')
	#if Nran==87**3: ax1.plot(x,yc,color=other,marker='s',linestyle='-.',label='CCVT')
	#ax1.plot(x,yg,color=blue,linestyle='--',label='Glass')
	if Nran==2*87**3: ax1.set_ylabel(r'$\Delta \bar{\xi_0(r)}$',fontsize=fs)
	ax1.set_xscale('log')
	#ax1.set_yscale('log')
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width, box.height*.8])
	#if Nran==2*87**3: ax1.legend(fontsize=lfs,loc='lower right')
	ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	if Nran != 2*87**3: ax1.set_yticklabels([])

	#ax1.set_xlim(-0.01,.01)
	ax1.set_ylim(-.01,1.5)

	#-----------------------------------------------------------------------
	#PLOT 2 - xi_2 Dipole
	#-----------------------------------------------------------------------
	if Nran==2*87**3: 
		yran = ran_1['xi2']-xi_noran['xi2']
		yrc  = rc_1['xi2']-xi_noran['xi2']
		yrs  = rs_1['xi2']-xi_noran['xi2']
		#yg   = g_1['xi2']-xi_noran['xi2']
		yzr  = zr_1['xi2']-xi_noran['xi2']
		#yc = c_1['xi2']-xi_noran['xi2']
	if Nran==4*87**3: 
		yran = ran_2['xi2']-xi_noran['xi2']
		yrc  = rc_2['xi2']-xi_noran['xi2']
		yrs  = rs_2['xi2']-xi_noran['xi2']
		#yg   = g_2['xi2']-xi_noran['xi2']
		yzr  = zr_2['xi2']-xi_noran['xi2']
	if Nran==8*87**3: 
		yran = ran_4['xi2']-xi_noran['xi2']
		yrc  = rc_4['xi2']-xi_noran['xi2']
		yrs  = rs_4['xi2']-xi_noran['xi2']
		#yg   = g_4['xi2']-xi_noran['xi2']
		yzr  = zr_4['xi2']-xi_noran['xi2']
	if Nran==16*87**3: 
		yran = ran_8['xi2']-xi_noran['xi2']
		yrc  = rc_8['xi2']-xi_noran['xi2']
		yrs  = rs_8['xi2']-xi_noran['xi2']
		#yg   = g_8['xi2']-xi_noran['xi2']
		yzr  = zr_8['xi2']-xi_noran['xi2']


	ax2 = f.add_subplot(312,sharex=ax1)
	ax2.plot(x,yran,color='k',ls=':',linewidth=2)
	ax2.plot(x,yrc,color=green,linestyle='-.')
	#if Nran==87**3: ax2.plot(x,yc,color=other,marker='s',linestyle='-.')
	ax2.plot(x,yrs,color=blue,linestyle='--')
	ax2.plot(x,yzr,color=red,linestyle='-')
	#ax2.plot(x,yg,color=blue,linestyle='--')
	if Nran==2*87**3: ax2.set_ylabel(r'$\Delta \bar{\xi_2(r)}$',fontsize=fs)
	#ax2.set_xscale('log')
	#ax2.set_yscale('log')
	ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
	if Nran != 2*87**3: ax2.set_yticklabels([])
	#ax2.set_ylim(-.04,.07)
	ax2.set_ylim(-.04,1.3)

	#-----------------------------------------------------------------------
	#PLOT 3 - xi_4 Hexadecapole
	#-----------------------------------------------------------------------
	if Nran==2*87**3: 
		yran = ran_1['xi4']-xi_noran['xi4']
		yrc  = rc_1['xi4']-xi_noran['xi4']
		yrs  = rs_1['xi4']-xi_noran['xi4']
		#yg   = g_1['xi4']-xi_noran['xi4']
		yzr  = zr_1['xi4']-xi_noran['xi4']
		#yc = c_1['xi4']-xi_noran['xi4']
	if Nran==4*87**3: 
		yran = ran_2['xi4']-xi_noran['xi4']
		yrc  = rc_2['xi4']-xi_noran['xi4']
		yrs  = rs_2['xi4']-xi_noran['xi4']
		#yg   = g_2['xi4']-xi_noran['xi4']
		yzr  = zr_2['xi4']-xi_noran['xi4']
	if Nran==8*87**3: 
		yran = ran_4['xi4']-xi_noran['xi4']
		yrc  = rc_4['xi4']-xi_noran['xi4']
		yrs  = rs_4['xi4']-xi_noran['xi4']
		#yg   = g_4['xi4']-xi_noran['xi4']
		yzr  = zr_4['xi4']-xi_noran['xi4']
	if Nran==16*87**3: 
		yran = ran_8['xi4']-xi_noran['xi4']
		yrc  = rc_8['xi4']-xi_noran['xi4']
		yrs  = rs_8['xi4']-xi_noran['xi4']
		#yg   = g_8['xi4']-xi_noran['xi4']
		yzr  = zr_8['xi4']-xi_noran['xi4']

	ax3 = f.add_subplot(313, sharex=ax2,sharey=ax2)
	ax3.plot(x,yran,color='k',ls=':',linewidth=2)
	ax3.plot(x,yrc,color=green,linestyle='-.')
	ax3.plot(x,yrs,color=blue,linestyle='--')
	ax3.plot(x,yzr,color=red,linestyle='-')
	#if Nran==87**3: ax3.plot(x,yc,color=other,marker='s',linestyle='-.')
	#ax3.plot(x,yg,color=blue,linestyle='--')
	if Nran==2*87**3: ax3.set_ylabel(r'$\Delta \bar{\xi_4(r)}$',fontsize=fs)
	#ax3.set_xscale('log')
	#ax3.set_yscale('log')
	ax3.set_xlabel('$R\,[Mpc]$',fontsize=fs)
	#ax3.set_ylim(-.0004,.0006)
	if Nran != 2*87**3: ax3.set_yticklabels([])

	#f.suptitle(r'$R = {}Mpc$'.format(nr[ir]),fontsize=15)
	#f.subplots_adjust(hspace=0)
	f.tight_layout()
	f.savefig('../plots/smallscale/mean/mean_ss_{}.png'.format(Nran),dpi=f.dpi)
	if Nran==16*87**3: f.savefig('../plots/smallscale/mean/mean_ss_{}_Niter{}.png'.format(Nran,Niter),dpi=f.dpi)

	#plt.show()
	plt.close()

