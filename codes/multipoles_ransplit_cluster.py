"""
Probando hacer split de randoms

-Calcular la xi(r,mu)
-Construir el estimador xi_split
-Calcular los multipolos con get_xi0246 usando el estimador xi_split

Pre-cluster: haciendo cambios finales para mandar al cluster
"""
from __future__ import division
import os
import sys
import numpy as np
from nbodykit.lab import *
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import matplotlib.pyplot as pltDE


def get_xi0246_wsplit(nbins_m,nbins_s):

	xi_sm = corr_split
	dmu=1.0/nbins_m

	rs = DD.pairs.coords['r']
	mu = DD.pairs.coords['mu']

	xi_s0 = np.zeros(nbins_s)
	xi_s2 = np.zeros(nbins_s)
	xi_s4 = np.zeros(nbins_s)
	xi_s6 = np.zeros(nbins_s)

	sr = np.zeros(nbins_s)
	rm = np.zeros(nbins_m)

	l0 = 0.0
	l1 = 1.0
	l2 = 2.0
	l3 = 3.0

	for i in range(nbins_s):
		sr[i] = rs[i]
		for j in range(nbins_m):
			rm[j]=mu[j]
			xi_s0[i]  += (4.0*l0+1.0)*xi_sm[i,j]*1.0*dmu 
			xi_s2[i]  += (4.0*l1+1.0)*xi_sm[i,j]*((3*rm[j]**2 - 1.0)/2.0)*dmu
			xi_s4[i]  += (4.0*l2+1.0)*xi_sm[i,j]*((35*rm[j]**4 - 30*rm[j]**2 + 3.0)/8.0)*dmu
			xi_s6[i]  += (4.0*l3+1.0)*xi_sm[i,j]*((231*rm[j]**6 - 315*rm[j]**4 + 105*rm[j]**2 - 5)/16.0)*dmu
    
	return xi_s0, xi_s2, xi_s4, xi_s6
##################################################################################


i = int(os.environ['SGE_TASK_ID'])*2-2
#i=1 #For testing

bs = 1500.
nran = int(sys.argv[1])*87**3
#nran = 600000 #For testing
print nran

nbins_m = 30
nbins_s = 29
names=['xi_0','xi_2','xi_4','xi_6']

datafile = '../data/Galaxies_HOD_001_z0.57.dat'

dcat = CSVCatalog(datafile,names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])
dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])


seed_path='../data/primos.txt'
seed=np.fromfile(seed_path,dtype='int',count=1000,sep=' ')


###Calculo y escribo los xi_l

print '~~~~~~~~~~~~~~~~~~~~~~~~',i,i+1,'~~~~~~~~~~~~~~~~~~~~~~~~~'
print nran    
np.random.seed(seed[i])
aran1 = bs*np.random.rand(nran,3)
table_ran1 = Table(aran1, names=['x','y','z'])

np.random.seed(seed[i+1])
aran2 = bs*np.random.rand(nran,3)
table_ran2 = Table(aran2, names=['x','y','z'])

rcat1 = ArrayCatalog(table_ran1, BoxSize=[bs,bs,bs])
rcat1['Position'] = transform.StackColumns(rcat1['x'],rcat1['y'],rcat1['z'])

rcat2 = ArrayCatalog(table_ran2, BoxSize=[bs,bs,bs])
rcat2['Position'] = transform.StackColumns(rcat2['x'],rcat2['y'],rcat2['z'])


#corr1 = SimulationBox2PCF('2d',data1=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat1,BoxSize=[bs,bs,bs], periodic=True)

#DD = SimulationBoxPairCount('2d',first=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)


#corr2 = SimulationBox2PCF('2d',data1=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat2,BoxSize=[bs,bs,bs], periodic=True)

R1R1 = SimulationBoxPairCount('2d',first=rcat1,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)
R2R2 = SimulationBoxPairCount('2d',first=rcat2,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)


DD = SimulationBoxPairCount('2d',first=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)
DR1 = SimulationBoxPairCount('2d',first=dcat,second=rcat1,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)

DR2 = SimulationBoxPairCount('2d',first=dcat,second=rcat2,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,BoxSize=[bs,bs,bs],periodic=True)

DDpc = DD.pairs['npairs'] #pc: pair counts
DR1pc = DR1.pairs['npairs']
DR2pc = DR2.pairs['npairs']

DRpc = DR1pc+DR2pc


R1R1pc = R1R1.pairs['npairs']
R2R2pc = R2R2.pairs['npairs']
RRpc = R1R1pc+R2R2pc

Nr_ = rcat1.size
Nr = 2*Nr_ #No. of splits * Nr'
Nd = dcat.size        

DD_ = DDpc/Nd/(Nd-1)
R1R1_ = (R1R1pc+R2R2pc)/Nr_/(Nr_-1)/2
DR1_ = (DR1pc+DR2pc)/Nr_/Nd/2
corr_split = (DD_ - 2*DR1_ + R1R1_)/R1R1_ 

xi_l = get_xi0246_wsplit(nbins_m,nbins_s)

wfile = '../data/out/ransplit/xi_l_ransplit_{}_{}-{}.txt'.format(nran,i,i+1)

ascii.write(xi_l,wfile,overwrite=True,format='no_header')


