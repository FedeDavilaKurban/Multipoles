"""
Multipoles, Errormeans, and Variances for all scales
---
LOG

-26Feb2020
Modifico este viejo codigo para incluir nuevos metodos (ransplit, zelrec, y ccvt)


"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

#nran1 = 87**3
nran = 87**3
Nran = nran

names=['xi_0','xi_2','xi_4','xi_6']


xi_l_ran = []
xi_l_rancross = []
xi_l_glass = []
xi_l_zelrec = []
xi_l_ransplit = []
xi_l_ccvt = []

for i in range(200):
    xi_l_ran.append( ascii.read('../data/out/ran/xi_l_{}_{}.txt'.format(nran,i),names=names))
for i in range(0,400,2):
    xi_l_rancross.append( ascii.read('../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
    xi_l_ransplit.append( ascii.read('../data/out/ransplit/xi_l_ransplit_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
    xi_l_glass.append( ascii.read('../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
    xi_l_zelrec.append( ascii.read('../data/out/zelrec/xi_l_zelrec_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
if nran==87**3:
	for i in range(1,120,2):
		xi_l_ccvt.append( ascii.read('../data/out/ccvt/xi_l_ccvt_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

xi_l = ascii.read('../data/xi_l_noran.txt',names=names)

#Agrupo los xi_l
xi_0_ran=[]
xi_2_ran=[]
xi_4_ran=[]
for i in range(len(xi_l_ran)):
    xi_0_ran.append( xi_l_ran[i]['xi_0'] )
    xi_2_ran.append( xi_l_ran[i]['xi_2'] )
    xi_4_ran.append( xi_l_ran[i]['xi_4'] )

xi_0_rancross=[]
xi_2_rancross=[]
xi_4_rancross=[]
for i in range(len(xi_l_rancross)):
    xi_0_rancross.append( xi_l_rancross[i]['xi_0'] )
    xi_2_rancross.append( xi_l_rancross[i]['xi_2'] )
    xi_4_rancross.append( xi_l_rancross[i]['xi_4'] )

xi_0_ransplit=[]
xi_2_ransplit=[]
xi_4_ransplit=[]
for i in range(len(xi_l_ransplit)):
    xi_0_ransplit.append( xi_l_ransplit[i]['xi_0'] )
    xi_2_ransplit.append( xi_l_ransplit[i]['xi_2'] )
    xi_4_ransplit.append( xi_l_ransplit[i]['xi_4'] )

xi_0_glass=[]
xi_2_glass=[]
xi_4_glass=[]
for i in range(len(xi_l_glass)):
    xi_0_glass.append( xi_l_glass[i]['xi_0'] )
    xi_2_glass.append( xi_l_glass[i]['xi_2'] )
    xi_4_glass.append( xi_l_glass[i]['xi_4'] )

xi_0_zelrec=[]
xi_2_zelrec=[]
xi_4_zelrec=[]
for i in range(len(xi_l_zelrec)):
    xi_0_zelrec.append( xi_l_zelrec[i]['xi_0'] )
    xi_2_zelrec.append( xi_l_zelrec[i]['xi_2'] )
    xi_4_zelrec.append( xi_l_zelrec[i]['xi_4'] )

if nran==87**3:
	xi_0_ccvt=[]
	xi_2_ccvt=[]
	xi_4_ccvt=[]
	for i in range(len(xi_l_ccvt)):
		xi_0_ccvt.append( xi_l_ccvt[i]['xi_0'] )
		xi_2_ccvt.append( xi_l_ccvt[i]['xi_2'] )
		xi_4_ccvt.append( xi_l_ccvt[i]['xi_4'] )
   

#Calculo las medias
xi_0mean_ran=[]
xi_2mean_ran=[]
xi_4mean_ran=[]

xi_0var_ran=[]
xi_2var_ran=[]
xi_4var_ran=[]

xi_0mean_rancross=[]
xi_2mean_rancross=[]
xi_4mean_rancross=[]

xi_0var_rancross=[]
xi_2var_rancross=[]
xi_4var_rancross=[]

xi_0mean_ransplit=[]
xi_2mean_ransplit=[]
xi_4mean_ransplit=[]

xi_0var_ransplit=[]
xi_2var_ransplit=[]
xi_4var_ransplit=[]

xi_0mean_glass=[]
xi_2mean_glass=[]
xi_4mean_glass=[]

xi_0var_glass=[]
xi_2var_glass=[]
xi_4var_glass=[]

xi_0mean_zelrec=[]
xi_2mean_zelrec=[]
xi_4mean_zelrec=[]

xi_0var_zelrec=[]
xi_2var_zelrec=[]
xi_4var_zelrec=[]

if nran==87**3:
	xi_0mean_ccvt=[]
	xi_2mean_ccvt=[]
	xi_4mean_ccvt=[]

	xi_0var_ccvt=[]
	xi_2var_ccvt=[]
	xi_4var_ccvt=[]

nbins_s = 29

for i in range(nbins_s):
	#Random
    xi_0mean_ran.append( np.mean([xi[i] for xi in xi_0_ran]) )
    xi_2mean_ran.append( np.mean([xi[i] for xi in xi_2_ran]) )
    xi_4mean_ran.append( np.mean([xi[i] for xi in xi_4_ran]) )

    xi_0var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_0_ran],ddof=1)) )
    xi_2var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_2_ran],ddof=1)) )
    xi_4var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_4_ran],ddof=1)) )
    
    #RanCross
    xi_0mean_rancross.append( np.mean([xi[i] for xi in xi_0_rancross]) )
    xi_2mean_rancross.append( np.mean([xi[i] for xi in xi_2_rancross]) )
    xi_4mean_rancross.append( np.mean([xi[i] for xi in xi_4_rancross]) )

    xi_0var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_0_rancross],ddof=1)) )
    xi_2var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_2_rancross],ddof=1)) )
    xi_4var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_4_rancross],ddof=1)) )

	#RanSplit
    xi_0mean_ransplit.append( np.mean([xi[i] for xi in xi_0_ransplit]) )
    xi_2mean_ransplit.append( np.mean([xi[i] for xi in xi_2_ransplit]) )
    xi_4mean_ransplit.append( np.mean([xi[i] for xi in xi_4_ransplit]) )

    xi_0var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_0_ransplit],ddof=1)) )
    xi_2var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_2_ransplit],ddof=1)) )
    xi_4var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_4_ransplit],ddof=1)) )

	#Glass
    xi_0mean_glass.append( np.mean([xi[i] for xi in xi_0_glass]) )
    xi_2mean_glass.append( np.mean([xi[i] for xi in xi_2_glass]) )
    xi_4mean_glass.append( np.mean([xi[i] for xi in xi_4_glass]) )

    xi_0var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_0_glass],ddof=1)) )
    xi_2var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_2_glass],ddof=1)) )
    xi_4var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_4_glass],ddof=1)) )
    
    #ZelRec
    xi_0mean_zelrec.append( np.mean([xi[i] for xi in xi_0_zelrec]) )
    xi_2mean_zelrec.append( np.mean([xi[i] for xi in xi_2_zelrec]) )
    xi_4mean_zelrec.append( np.mean([xi[i] for xi in xi_4_zelrec]) )

    xi_0var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_0_zelrec],ddof=1)) )
    xi_2var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_2_zelrec],ddof=1)) )
    xi_4var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_4_zelrec],ddof=1)) )
    
    if nran==87**3:
        #CCVT
        xi_0mean_ccvt.append( np.mean([xi[i] for xi in xi_0_ccvt]) )
        xi_2mean_ccvt.append( np.mean([xi[i] for xi in xi_2_ccvt]) )
        xi_4mean_ccvt.append( np.mean([xi[i] for xi in xi_4_ccvt]) )

        xi_0var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_0_ccvt],ddof=1)) )
        xi_2var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_2_ccvt],ddof=1)) )
        xi_4var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_4_ccvt],ddof=1)) )
    	
    
nr = 2.5+np.linspace(5.,150.,30)[:-1]

#########
#Writing
########

var_ran = Table([xi_0var_ran,xi_2var_ran,xi_4var_ran,nr], names=('xi0', 'xi2','xi4','r'))
ascii.write(var_ran, '../data/out/var_xi/var_ran_{}.txt'.format(nran), overwrite=True)
#Testing... faltaria escribir los otros

"""
MULTIPOLES

"""

mps=[0]

if 0 in mps:
    plt.errorbar(nr,(nr**2)*xi_l['xi_0'],ls='-',color='r',label=r'Analytic')
if 2 in mps:
    plt.errorbar(nr,(nr**2)*xi_l['xi_2'],ls='-.',color='r',label=r'Analytic')
if 4 in mps:
    plt.errorbar(nr,(nr**2)*xi_l['xi_4'],ls='--',color='r',label=r'Analytic')


if 0 in mps:
    plt.errorbar(nr,(nr**2)*xi_0mean_ran,yerr=(nr**2)*np.array(xi_0var_ran),ls='-',color=seaborn.color_palette()[0],label=r'$\xi_0(r)\, Ran$')
if 2 in mps:
    plt.errorbar(nr,(nr**2)*(-np.array(xi_2mean_ran)),yerr=(nr**2)*np.array(xi_2var_ran),ls='-.',color=seaborn.color_palette()[0],label=r'$\xi_2(r)\, Ran$')
if 4 in mps:
    plt.errorbar(nr,(nr**2)*xi_4mean_ran,yerr=(nr**2)*np.array(xi_4var_ran),ls='--',color=seaborn.color_palette()[0],label=r'$\xi_4(r)\, Ran$')

if 0 in mps:
    plt.errorbar(nr,(nr**2)*xi_0mean_rancross,yerr=(nr**2)*np.array(xi_0var_rancross),ls='-',color=seaborn.color_palette()[1],label=r'$\xi_0(r)\, RanCross$')
if 2 in mps:
    plt.errorbar(nr,(nr**2)*(-np.array(xi_2mean_rancross)),yerr=(nr**2)*np.array(xi_2var_rancross),ls='-.',color=seaborn.color_palette()[1],label=r'$\xi_2(r)\, RanCross$')
if 4 in mps:
    plt.errorbar(nr,(nr**2)*xi_4mean_rancross,yerr=(nr**2)*np.array(xi_4var_rancross),ls='--',color=seaborn.color_palette()[1],label=r'$\xi_4(r)\, RanCross$')

if 0 in mps:
    plt.errorbar(nr,(nr**2)*xi_0mean_glass,yerr=(nr**2)*np.array(xi_0var_glass),ls='-',color=seaborn.color_palette()[2],label=r'$\xi_0(r)\, Glass$')
if 2 in mps:
    plt.errorbar(nr,(nr**2)*(-np.array(xi_2mean_glass)),yerr=(nr**2)*np.array(xi_2var_glass),ls='-.',color=seaborn.color_palette()[2],label=r'$\xi_2(r)\, Glass$')
if 4 in mps:
    plt.errorbar(nr,(nr**2)*xi_4mean_glass,yerr=(nr**2)*np.array(xi_4var_glass),ls='--',color=seaborn.color_palette()[2],label=r'$\xi_4(r)\, Glass$')

#plt.xscale('log')
plt.legend()
plt.title(r'$Nran={}$'.format(nran))
plt.show()


"""
VARIANCES

"""

mps=[0]

#Random
if 0 in mps:
    plt.plot(nr,xi_0var_ran,label=r'$\sigma_{\xi_0(r)}\, Ran$',ls='-',color=seaborn.color_palette()[0])
if 2 in mps:
    plt.plot(nr,xi_2var_ran,label=r'$\sigma_{\xi_2(r)}\, Ran$',ls='-',color=seaborn.color_palette()[1])
if 4 in mps:
    plt.plot(nr,xi_4var_ran,label=r'$\sigma_{\xi_4(r)}\, Ran$',ls='-',color=seaborn.color_palette()[2])

#RanCross
if 0 in mps:
    plt.plot(nr,xi_0var_rancross,label=r'$\sigma_{\xi_0(r)}\, RanCross$',ls='--',color=seaborn.color_palette()[0])
if 2 in mps:
    plt.plot(nr,xi_2var_rancross,label=r'$\sigma_{\xi_2(r)}\, RanCross$',ls='--',color=seaborn.color_palette()[1])
if 4 in mps:
    plt.plot(nr,xi_4var_rancross,label=r'$\sigma_{\xi_4(r)}\, RanCross$',ls='--',color=seaborn.color_palette()[2])

#RanSplit
if 0 in mps:
    plt.plot(nr,xi_0var_ransplit,label=r'$\sigma_{\xi_0(r)}\, RanSplit$',ls=(0,(1,1)),color=seaborn.color_palette()[0])
if 2 in mps:
    plt.plot(nr,xi_2var_ransplit,label=r'$\sigma_{\xi_2(r)}\, RanSplit$',ls=(0,(1,1)),color=seaborn.color_palette()[1])
if 4 in mps:
    plt.plot(nr,xi_4var_ransplit,label=r'$\sigma_{\xi_4(r)}\, RanSplit$',ls=(0,(1,1)),color=seaborn.color_palette()[2])

#Glass
if 0 in mps:
    plt.plot(nr,xi_0var_glass,label=r'$\sigma_{\xi_0(r)}\, Glass$',ls='-.',color=seaborn.color_palette()[0])
if 2 in mps:
    plt.plot(nr,xi_2var_glass,label=r'$\sigma_{\xi_2(r)}\, Glass$',ls='-.',color=seaborn.color_palette()[1])
if 4 in mps:
    plt.plot(nr,xi_4var_glass,label=r'$\sigma_{\xi_4(r)}\, Glass$',ls='-.',color=seaborn.color_palette()[2])

#ZelRec
if 0 in mps:
    plt.plot(nr,xi_0var_zelrec,label=r'$\sigma_{\xi_0(r)}\, Zel. Rec.$',ls=(0,(5,5)),color=seaborn.color_palette()[0])
if 2 in mps:
    plt.plot(nr,xi_2var_zelrec,label=r'$\sigma_{\xi_2(r)}\, Zel. Rec.$',ls=(0,(5,5)),color=seaborn.color_palette()[1])
if 4 in mps:
    plt.plot(nr,xi_4var_zelrec,label=r'$\sigma_{\xi_4(r)}\, Zel. Rec.$',ls=(0,(5,5)),color=seaborn.color_palette()[2])

#CCVT
if nran==87**3:
	if 0 in mps:
		plt.plot(nr,xi_0var_ccvt,label=r'$\sigma_{\xi_0(r)}\, CCVT$',ls=(0,(1,10)),color=seaborn.color_palette()[0])
	if 2 in mps:
		plt.plot(nr,xi_2var_ccvt,label=r'$\sigma_{\xi_2(r)}\, CCVT$',ls=(0,(1,10)),color=seaborn.color_palette()[1])
	if 4 in mps:
		plt.plot(nr,xi_4var_ccvt,label=r'$\sigma_{\xi_4(r)}\, CCVT$',ls=(0,(1,10)),color=seaborn.color_palette()[2])


#plt.xscale('log')
plt.yscale('log')
plt.title(r'$Nran = {}$'.format(nran))
plt.legend()
plt.show()
