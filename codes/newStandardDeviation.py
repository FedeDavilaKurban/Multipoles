"""
INTRO

21Feb2020
Luego de un frustrante fracaso en el intento de modificar el plot de las varianzas (var_vs_nran_v5.5 en adelante)
decido rehacer el código que las calcula.

El objetivo es hacer un código que:
1) Calcule Varianzas
2) Escriba en un archivo

para despues levantar esos archivos con un código que sólo haga plots

---
LOG

-21Feb2020
Primero voy a hacer una prueba: calcular sólo para los Random(N=200), 87^3


"""
import glob 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales

Nran = 87**3 #Number of points in random/homogeneous catalogue

########################################################################################


########################################################################################	
xil_ran = [] # "xil" meaning $\xi_l$, aka multipoles
xi0_ran = [] # Monopole
xi2_ran = [] # Quadrupole
xi4_ran = [] # Hexadecapole
N = len(glob.glob('../data/out/ran/*{}*'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_Ran is: {} (should be 200)'.format(N))
for iN in range(N):
	xil_ran.append( ascii.read('../data/out/ran/xi_l_{}_{}.txt'.format(Nran,iN),names=['xi_0','xi_2','xi_4','xi_6']) )	
	xi0_ran.append( xil_ran[iN]['xi_0'].data )
	xi2_ran.append( xil_ran[iN]['xi_2'].data )
	xi4_ran.append( xil_ran[iN]['xi_4'].data )

#Calculate Standard Deviation  
std_xi0_ran = []
std_xi2_ran = []	
std_xi4_ran = []		
for ir in range(len(nr)):
	std_xi0_ran.append( np.sqrt( np.var([item[ir] for item in xi0_ran],ddof=1) ))
	std_xi2_ran.append( np.sqrt( np.var([item[ir] for item in xi4_ran],ddof=1) ))
	std_xi4_ran.append( np.sqrt( np.var([item[ir] for item in xi2_ran],ddof=1) ))

std_ran = Table([std_xi0_ran,std_xi2_ran,std_xi4_ran,nr], names=('xi0', 'xi2','xi4','r'))
ascii.write(std_ran, '../data/out/std_xi/std_ran.txt')
#########################################################################################	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
 
