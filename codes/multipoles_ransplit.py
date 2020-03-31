"""
Probando hacer split de randoms

-Calcular la xi(r,mu)
-Construir el estimador xi_split
-Calcular los multipolos con get_xi0246 usando el estimador xi_split
"""
from __future__ import division
import os
import sys
import numpy as np
from nbodykit.lab import *
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import matplotlib.pyplot as plt

def get_xi0246(corr,nbins_m,nbins_s,est):

    if est=='ls': xi_sm = corr.corr.data['corr'] 
    if est=='split': xi_sm = corr_split
    ##Modificado para que lea los valores de entrada
    ##nbins_m=30 # number of bins in mu
    ##nbins_s=29 # number of bins in s
    dmu=1.0/nbins_m
    
    rs = corr.D1D2.coords['r']
    mu = corr.D1D2.coords['mu']
    
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


#i = int(os.environ['SGE_TASK_ID'])*2-2
i=1 #For testing

bs = 1500.
#nran = int(sys.argv[1])*87**3
nran = 87**3 #For testing
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

corr = SimulationBox2PCF('2d',data1=dcat,data2=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat1,randoms2=rcat2,BoxSize=[bs,bs,bs], periodic=True)

xi_l = get_xi0246(corr,nbins_m,nbins_s,'ls')


##############################################
#~~~~~~~~~~~~~~~~~~PRUEBA~~~~~~~~~~~~~~~~~~~~~
#Probando armar el estimador LS a mano y ver si da lo mismo
xi0_automatico=xi_l[0] #Guardo la correlacion que me da el nbodykit

DD = corr.D1D2.data['npairs']
R1R2 = corr.R1R2.data['npairs']
DR1 = corr.D2R1.data['npairs']
DR2 = corr.D1R2.data['npairs']

Nr1 = corr.randoms1.size
Nr2 = corr.randoms2.size
Nd = corr.data1.size

corr_ls_manual = DD*Nr1*Nr2/(R1R2*Nd*(Nd-1)) - DR1*Nr1*Nr2/(R1R2*Nd*Nr1) - DR2*Nr1*Nr2/(R1R2*Nd*Nr2) + 1. 


#Calculo "manual" de los xi_l
 
xi_sm = corr_ls_manual
dmu=1.0/nbins_m

rs = corr.D1D2.coords['r']
mu = corr.D1D2.coords['mu']

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

plt.plot(corr.D1D2.coords['r'],corr.D1D2.coords['r']**2*xi_l[0],label='automatico')
plt.plot(corr.D1D2.coords['r'],corr.D1D2.coords['r']**2*xi_s0,label='manual')
plt.legend()
plt.show()
 

#~~~~~~~~~~~~~~~FIN DE PRUEBA~~~~~~~~~~~~~~~~~
##############################################
 

##############################################
#~~~~~~~~~~~~~~~~~~PRUEBA2~~~~~~~~~~~~~~~~~~~~~
#Probando armar el estimador LS_split 

corr1 = SimulationBox2PCF('2d',data1=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat1,BoxSize=[bs,bs,bs], periodic=True)

corr2 = SimulationBox2PCF('2d',data1=dcat,edges=np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat2,BoxSize=[bs,bs,bs], periodic=True)


DD = corr1.D1D2.data['npairs'] #=corr2.D1D2
DR1 = corr1.D1R2.data['npairs'] #=corr1.D2R1
DR2 = corr2.D1R2.data['npairs'] #=corr2.D2R1
DR = DR1+DR2

R1R1 = corr1.R1R2.data['npairs']
R2R2 = corr2.R1R2.data['npairs']
RR = R1R1+R2R2

Nr_ = corr1.randoms1.size
Nr = 2*Nr_ #No. of splits * Nr'
Nd = corr1.data1.size        

#corr_split = DD*Nr*(Nr_-1)/(RR*Nd*(Nd-1)) - DR*(Nr_-1)/(RR*Nd) + 1

#corr_split = DD*Nr_*(Nr_-1)/(R1R1*Nd*(Nd-1)) - DR1*(Nr_-1)/(R1R1*Nd) + 1 
#corr_split = DD/RR - DR/RR + 1

DD_ = DD/Nd/(Nd-1)
R1R1_ = (R1R1+R2R2)/Nr_/(Nr_-1)/2
DR1_ = (DR1+DR2)/Nr_/Nd/2
corr_split = (DD_ - 2*DR1_ + R1R1_)/R1R1_ 


#DD_ = DD/(corr1.data1.size*(corr1.data1.size-1))

#DR_ = DR/(corr1.data1.size*corr1.randoms1.size+corr1.data1.size*corr2.randoms1.size)

#RR_ = RR/(corr1.randoms1.size*(corr1.randoms1.size-1)+corr2.randoms1.size*(corr2.randoms1.size-1))

#corr_split = DD_/RR_ - DR_/RR_ + 1


#~~~~Calculo "manual" de los xi_l~~~~~
 
xi_sm = corr_split
dmu=1.0/nbins_m

rs = corr1.D1D2.coords['r']
mu = corr1.D1D2.coords['mu']

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


plt.plot(corr.D1D2.coords['r'],corr.D1D2.coords['r']**2*xi_l[0],label='automatico')
plt.plot(corr1.D1D2.coords['r'],corr1.D1D2.coords['r']**2*xi_s0,label='split')
plt.legend()
plt.show()

#Plot the first mu just for comparison
plt.plot(corr1.corr.data['r'][:,0],corr1.corr.data['corr'][:,0],label='LS')
plt.plot(corr1.corr.data['r'][:,0],corr_split[:,0],label='Split')
plt.legend()
plt.show()


 

#~~~~~~~~~~~~~~~FIN DE PRUEBA2~~~~~~~~~~~~~~~~~
##############################################


wfile = '../data/out/rancross/xi_l_ransplit_{}_{}-{}.txt'.format(nran,i,i+1)

ascii.write(xi_l,wfile,overwrite=True,format='no_header')
