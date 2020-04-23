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
Primero voy a hacer una prueba: calcular sólo para los Random(N=200), 87^3 -- Éxito
Le cambié el nombre de newVariances.py --> StandardDeviation.py
Probar extender el código para todos los métodos -- Éxito
Hacer el código para los plots: StandardDeviation_plots.py

-26Feb2020
Los resultados no son iguales a los que obtenía con var_vs_nran_v5... todavía no sé por qué

Comparando con variances.py. Xi_l_ran y Xi_l_rancross dan igual en ambos programas, pero xi_l_glass no... investigar -- Ok. En este no estaba leyendo los 400 Glass

"""
import glob 
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

nr = 2.5+np.linspace(5.,150.,30)[:-1] #Scales


# Nran == Number of points in random/homogeneous catalogue

for Nran in [87**3,2*87**3,4*87**3,8*87**3]:
	print('Nran=',Nran)
	
	#print('Randoms')
	#exec(open("./StandardDeviation_Ran.py").read())
	
	#print('RanCross')
	#exec(open("./StandardDeviation_RanCross.py").read())

	#print('RanSplit')
	#exec(open("./StandardDeviation_RanSplit.py").read())

	#print('Glass')
	#exec(open("./StandardDeviation_Glass.py").read())
	
	print('ZelRec')
	exec(open("./StandardDeviation_ZelRec.py").read())
	
	#if Nran==87**3:
	#	print('CCVT')
	#	exec(open("./StandardDeviation_Ccvt.py").read())

for Nran in [16*87**3]:
	print('Randoms16')
	exec(open("./StandardDeviation_Ran.py").read())

