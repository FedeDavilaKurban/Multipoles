def get_xi0246(corr,nbins_m,nbins_s):

    xi_sm = corr.corr.data['corr']
    
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


bs = 1500.


f = CSVCatalog('/media/fdavila/TOSHIBA EXT/minerva/tracers/Minerva_HOD_galaxies_01/Galaxies_HOD_001_z0.57.dat',names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])

f['Position'] = transform.StackColumns(f['x'], f['y'], f['z'])
#rcat = UniformCatalog(nbar=f.size/bs**3, BoxSize=bs, seed=42)

nbins_m = 30
nbins_s = 29

corr = SimulationBox2PCF('2d',f,np.linspace(5.,150.,nbins_s+1),nbins_m,BoxSize=[bs,bs,bs], periodic=True)

xi_l = get_xi0246(corr,nbins_m,nbins_s)

nr = 2.5+np.linspace(5.,150.,30)[:-1]
for i in range(len(pl)-1):
    if i==1: plt.plot(nr,(nr**2)*(-xi_l[i]),label=r'$\xi_{}(k)$'.format(i*2))
    if i!=1: plt.plot(nr,(nr**2)*xi_l[i],label=r'$\xi_{}(k)$'.format(i*2))
#plt.xscale('log')
plt.legend()
plt.show()

