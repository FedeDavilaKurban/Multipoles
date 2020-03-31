"""
INTRO
---
LOG

-28Feb2020
Quiero calcular las medias. Copio los programas de StandardDeviation.py para hacer los cálculos


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
	
	print('Randoms')
	exec(open("./Mean_Ran.py").read())
	
	print('RanCross')
	exec(open("./Mean_RanCross.py").read())

	print('RanSplit')
	exec(open("./Mean_RanSplit.py").read())

	print('Glass')
	exec(open("./Mean_Glass.py").read())
	
	print('ZelRec')
	exec(open("./Mean_ZelRec.py").read())
	
	#if Nran==87**3:
	#	print('CCVT')
	#	exec(open("./Mean_Ccvt.py").read())

