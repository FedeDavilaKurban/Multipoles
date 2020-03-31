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

#niter = 100

bs = 1500.
#nran = .5 #*Ngal
nbins_m = 30
nbins_s = 29
names=['xi_0','xi_2','xi_4','xi_6']
kt=np.logspace(-3,0)

#seed=np.fromfile('/home/fdavila/Proyectos/Multipoles/data/primos.txt',dtype='int',count=niter,sep=' ')




#Calculo y escribo los xi_l
for i in ['010','020','030','040','050','060','070','080','090','100']:
    
    print '~~~~~~~~~~~~~~~~~~~~~~~~',i,'~~~~~~~~~~~~~~~~~~~~~~~~~'
        
    dcat = CSVCatalog('/home/fdavila/Proyectos/Multipoles/data/iter/iter_{}.txt'.format(i),names=['x','y','z'])
    dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])
    
    real_mesh = dcat.to_mesh(compensated=True, window='tsc', position='Position', Nmesh=256, BoxSize=bs)
    
    r = FFTPower(real_mesh, mode='1d',dk=0.05)
    Pk = r.power

    plt.loglog(Pk['k'], Pk['power'].real,label=i)# - Pk.attrs['shotnoise'])
    #plt.loglog(Pk['k'], Pk['power'].real - Pk.attrs['shotnoise'],label=i)


plt.loglog(kt,(10**(17.5)/(float(dcat.size))**2.32)*kt**4,c=c,ls='-.') 

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
#plt.xlim(1E2, 1E3)
plt.legend()
plt.show()






