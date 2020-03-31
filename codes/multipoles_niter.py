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

niter = 100

bs = 1500.
nran = 87**3
nbins_m = 30
nbins_s = 29
names=['xi_0','xi_2','xi_4','xi_6']

seed_path='/home/fdavila/Proyectos/Multipoles/data/primos.txt'
data_path='/home/fdavila/Proyectos/Multipoles/data/Galaxies_HOD_001_z0.57.dat'


###Para correr en Clemente/Mirta
#seed_path='/home/fdavilakurban/Multipoles/data/primos.txt'
#data_path='/home/fdavilakurban/Multipoles/data/Galaxies_HOD_001_z0.57.dat'


seed=np.fromfile(seed_path,dtype='int',count=niter,sep=' ')


dcat = CSVCatalog(data_path,names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])

dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])


#Calculo y escribo los xi_l
for i in range(66,niter):
    
    print '~~~~~~~~~~~~~~~~~~~~~~~~',i,'~~~~~~~~~~~~~~~~~~~~~~~~~'
    
    rcat = UniformCatalog(nbar=nran/bs**3, BoxSize=[bs,bs,bs], seed=seed[i])

    corr = SimulationBox2PCF('2d',dcat,np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat,BoxSize=[bs,bs,bs], periodic=True)

    xi_l = get_xi0246(corr,nbins_m,nbins_s)
    
    ascii.write(xi_l,'/home/fdavila/Proyectos/Multipoles/data/xi_l_{}_{}_test.txt'.format(nran,seed[i]),overwrite=True,format='no_header')
    
    ###Para correr en Clemente/Mirta
    #ascii.write(xi_l,'/home/fdavilakurban/Multipoles/data/xi_l_{}_{}.txt'.format(nran,seed[i]),overwrite=True,format='no_header')
    

#Leo los xi_l (si no los calcule ya)

niter = 100

xi_l = []
for i in range(niter):
    xi_l.append( ascii.read('/home/fdavila/Proyectos/Multipoles/data/xi_l_{}_{}_test.txt'.format(nran,seed[i]),names=['xi_0','xi_2','xi_4','xi_6']) )
    

#Agrupo los xi_l
xi_0_ran=[]
xi_2_ran=[]
xi_4_ran=[]
xi_6_ran=[]
for i in range(niter):
    xi_0_ran.append( xi_l[i]['xi_0'] )
    xi_2_ran.append( xi_l[i]['xi_2'] )
    xi_4_ran.append( xi_l[i]['xi_4'] )
    xi_6_ran.append( xi_l[i]['xi_6'] )
    

#Calculo las medias
xi_0bar_ran=[]
xi_2bar_ran=[]
xi_4bar_ran=[]
xi_6bar_ran=[]

xi_0var_ran=[]
xi_2var_ran=[]
xi_4var_ran=[]
xi_6var_ran=[]

for i in range(nbins_s):
    xi_0bar_ran.append( np.mean([xi[i] for xi in xi_0_ran]) )
    xi_2bar_ran.append( np.mean([xi[i] for xi in xi_2_ran]) )
    xi_4bar_ran.append( np.mean([xi[i] for xi in xi_4_ran]) )
    xi_6bar_ran.append( np.mean([xi[i] for xi in xi_6_ran]) )

    xi_0var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_0_ran],ddof=1)) )
    xi_2var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_2_ran],ddof=1)) )
    xi_4var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_4_ran],ddof=1)) )
    xi_6var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_6_ran],ddof=1)) )
    
    
    

nr = 2.5+np.linspace(5.,150.,30)[:-1]


plt.errorbar(nr,(nr**2)*xi_0bar_ran,yerr=(nr**2)*np.array(xi_0var_ran),label=r'$\overline{\xi_0(r)}\,Random$')
plt.errorbar(nr,(nr**2)*(-np.array(xi_2bar_ran)),yerr=(nr**2)*np.array(xi_2var_ran),label=r'$\overline{\xi_2(r)}\,Random$')
plt.errorbar(nr,(nr**2)*xi_4bar_ran,yerr=(nr**2)*np.array(xi_4var_ran),label=r'$\overline{\xi_4(r)}\,Random$')
#plt.xscale('log')
plt.legend()
plt.title(r'$Nran=\sqrt{2}*87^{3}=931263.87$')
plt.show()

plt.plot(nr,xi_0var_ran,label=r'$\sigma_{\xi_0(r)}\, Random$',color=seaborn.color_palette()[0])
plt.plot(nr,xi_2var_ran,label=r'$\sigma_{\xi_2(r)}\, Random$',color=seaborn.color_palette()[1])
plt.plot(nr,xi_4var_ran,label=r'$\sigma_{\xi_4(r)}\, Random$',color=seaborn.color_palette()[2])
plt.legend()
plt.show()


