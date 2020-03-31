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

#seed=np.fromfile('/home/fdavila/Proyectos/Multipoles/data/primos.txt',dtype='int',count=niter,sep=' ')


dcat = CSVCatalog('/home/fdavila/Proyectos/Multipoles/data/Galaxies_HOD_001_z0.57.dat',names=['x','y','z','vx','vy','vz','cenflag','Mhalo','CeinID','Xcmhalo','Ycmhalo','Zcmhalo','VXcmhalo','VYcmhalo','VZcmhalo'],usecols=['x','y','z'])

dcat['Position'] = transform.StackColumns(dcat['x'], dcat['y'], dcat['z'])

idxs = ['001','002','003','004','005','006','007','008','009','010']
idxs=['001']
#Calculo y escribo los xi_l
for i in idxs:
    
    print('~~~~~~~~~~~~~~~~~~~~~~~~',i,'~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    rcat = CSVCatalog('/home/fdavila/Proyectos/Multipoles/data/glass/glass_{}.dat'.format(i),names=['x','y','z'])
    rcat['Position'] = transform.StackColumns(rcat['x'], rcat['y'], rcat['z'])

    corr = SimulationBox2PCF('2d',dcat,np.linspace(5.,150.,nbins_s+1),Nmu=nbins_m,randoms1=rcat,BoxSize=[bs,bs,bs], periodic=True)

    xi_l = get_xi0246(corr,nbins_m,nbins_s)
    
    ascii.write(xi_l,'/home/fdavila/Proyectos/Multipoles/data/glass/xi_l_{}_{}.txt'.format(rcat.size,i),overwrite=True,format='no_header')
    

#Leo los xi_l (si no los calcule ya)
xi_l = []
for i in idxs:
    xi_l.append( ascii.read('/home/fdavila/Proyectos/Multipoles/data/glass/xi_l_{}_{}.txt'.format(rcat.size,i),names=['xi_0','xi_2','xi_4','xi_6']))
    

#Agrupo los xi_l
xi_0=[]
xi_2=[]
xi_4=[]
xi_6=[]
for i in range(len(idxs)):
    xi_0.append( xi_l[i]['xi_0'] )
    xi_2.append( xi_l[i]['xi_2'] )
    xi_4.append( xi_l[i]['xi_4'] )
    xi_6.append( xi_l[i]['xi_6'] )
    

#Calculo las medias
xi_0bar=[]
xi_2bar=[]
xi_4bar=[]
xi_6bar=[]
for i in range(nbins_s):
    xi_0bar.append( np.mean([xi[i] for xi in xi_0]) )
    xi_2bar.append( np.mean([xi[i] for xi in xi_2]) )
    xi_4bar.append( np.mean([xi[i] for xi in xi_4]) )
    xi_6bar.append( np.mean([xi[i] for xi in xi_6]) )

xi_0var=[]
xi_2var=[]
xi_4var=[]
xi_6var=[]
for i in range(nbins_s):
    xi_0var.append( np.sqrt(np.var([xi[i] for xi in xi_0],ddof=1)) )
    xi_2var.append( np.sqrt(np.var([xi[i] for xi in xi_2],ddof=1)) )
    xi_4var.append( np.sqrt(np.var([xi[i] for xi in xi_4],ddof=1)) )
    xi_6var.append( np.sqrt(np.var([xi[i] for xi in xi_6],ddof=1)) )
    
    
    
    
    
    
    
nr = 2.5+np.linspace(5.,150.,30)[:-1]


plt.errorbar(nr,(nr**2)*xi_0bar,yerr=(nr**2)*np.array(xi_0var),label=r'$\overline{\xi_0(r)} Glass$')
plt.errorbar(nr,(nr**2)*(-np.array(xi_2bar)),yerr=(nr**2)*np.array(xi_2var),label=r'$\overline{\xi_2(r)} Glass$')
plt.errorbar(nr,(nr**2)*xi_4bar,yerr=(nr**2)*np.array(xi_4var),label=r'$\overline{\xi_4(r)} Glass$')
#plt.xscale('log')
plt.legend()
plt.show()

