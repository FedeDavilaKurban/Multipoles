"""
Codigo para hacer plots de varianzas
Debe ser usado con var_vs_nran_v6.py
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy
from astropy.table import Table

#nr = 2.5+np.linspace(5.,150.,30)[:-1]

for nran in ['1x87']

	ran = ascii.read('variances_ran_{}.txt'.format(nran))
	zelrec = ascii.read('variances_zelrec_{}.txt'.format(nran))
	glass = ascii.read('variances_glass_{}.txt'.format(nran))

	ratio_g = glass['xi0']/ran['xi0']

	plt.plot(nr,ratio_g,marker='o',linestyle='--')


#plt.hlines(1.,7.,148,linestyles='dashed')
#plt.ylim([.8,1.4])
plt.xscale('log')
#plt.title('{}^3'.format(col))
plt.ylabel(r'$\sigma_{Glass}/\sigma_{ZR}$')
plt.xlabel('Nran')
#plt.savefig('ratio_{}.png'.format(col))
#plt.close()

plt.plot(nr,zelrec['xi0'],marker='o',linestyle='--')
plt.title('{}^3'.format(nran))
plt.savefig('zrvar{}test.png'.format(nran))
