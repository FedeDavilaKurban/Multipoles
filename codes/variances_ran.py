seed_path='../data/primos.txt'
seed=np.fromfile(seed_path,dtype='int',count=niter,sep=' ')
names=['xi_0','xi_2','xi_4','xi_6']


xi_l_1 = []
xi_l_2 = []
xi_l_3 = []

niter=50
for nran in [658503,931263.873463,658503.111]:
 #   print nran
    if nran==658503:
        for i in range(niter):
            xi_l_1.append( ascii.read('../data/ran/xi_l_{}_{}_test.txt'.format(nran,seed[i]),names=names))

    if nran==931263.873463:
        for i in range(niter):
            xi_l_2.append( ascii.read('../data/ran/xi_l_{}_{}.txt'.format(nran,seed[i]),names=names))


    if nran==658503.111:
        inran=658503
        for i in range(0,niter*2,2):
            xi_l_3.append( ascii.read('../data/ran/xi_l_rancross_{}_{}-{}.txt'.format(inran,i,i+1),names=names))
    

#Agrupo los xi_l
xi_0_1=[]
xi_2_1=[]
xi_4_1=[]
xi_6_1=[]
for i in range(len(xi_l_1)):
    xi_0_1.append( xi_l_1[i]['xi_0'] )
    xi_2_1.append( xi_l_1[i]['xi_2'] )
    xi_4_1.append( xi_l_1[i]['xi_4'] )
    xi_6_1.append( xi_l_1[i]['xi_6'] )

xi_0_2=[]
xi_2_2=[]
xi_4_2=[]
xi_6_2=[]
for i in range(len(xi_l_2)):
    xi_0_2.append( xi_l_2[i]['xi_0'] )
    xi_2_2.append( xi_l_2[i]['xi_2'] )
    xi_4_2.append( xi_l_2[i]['xi_4'] )
    xi_6_2.append( xi_l_2[i]['xi_6'] )

xi_0_3=[]
xi_2_3=[]
xi_4_3=[]
xi_6_3=[]
for i in range(len(xi_l_3)):
    xi_0_3.append( xi_l_3[i]['xi_0'] )
    xi_2_3.append( xi_l_3[i]['xi_2'] )
    xi_4_3.append( xi_l_3[i]['xi_4'] )
    xi_6_3.append( xi_l_3[i]['xi_6'] )
    

#Calculo las medias
xi_0bar_1=[]
xi_2bar_1=[]
xi_4bar_1=[]
xi_6bar_1=[]

xi_0var_1=[]
xi_2var_1=[]
xi_4var_1=[]
xi_6var_1=[]

xi_0bar_2=[]
xi_2bar_2=[]
xi_4bar_2=[]
xi_6bar_2=[]

xi_0var_2=[]
xi_2var_2=[]
xi_4var_2=[]
xi_6var_2=[]

xi_0bar_3=[]
xi_2bar_3=[]
xi_4bar_3=[]
xi_6bar_3=[]

xi_0var_3=[]
xi_2var_3=[]
xi_4var_3=[]
xi_6var_3=[]

nbins_s = 29

for i in range(nbins_s):
    xi_0bar_1.append( np.mean([xi[i] for xi in xi_0_1]) )
    xi_2bar_1.append( np.mean([xi[i] for xi in xi_2_1]) )
    xi_4bar_1.append( np.mean([xi[i] for xi in xi_4_1]) )
    xi_6bar_1.append( np.mean([xi[i] for xi in xi_6_1]) )

    xi_0var_1.append( np.sqrt(np.var([xi[i] for xi in xi_0_1],ddof=1)) )
    xi_2var_1.append( np.sqrt(np.var([xi[i] for xi in xi_2_1],ddof=1)) )
    xi_4var_1.append( np.sqrt(np.var([xi[i] for xi in xi_4_1],ddof=1)) )
    xi_6var_1.append( np.sqrt(np.var([xi[i] for xi in xi_6_1],ddof=1)) )
    
    xi_0bar_2.append( np.mean([xi[i] for xi in xi_0_2]) )
    xi_2bar_2.append( np.mean([xi[i] for xi in xi_2_2]) )
    xi_4bar_2.append( np.mean([xi[i] for xi in xi_4_2]) )
    xi_6bar_2.append( np.mean([xi[i] for xi in xi_6_2]) )

    xi_0var_2.append( np.sqrt(np.var([xi[i] for xi in xi_0_2],ddof=1)) )
    xi_2var_2.append( np.sqrt(np.var([xi[i] for xi in xi_2_2],ddof=1)) )
    xi_4var_2.append( np.sqrt(np.var([xi[i] for xi in xi_4_2],ddof=1)) )
    xi_6var_2.append( np.sqrt(np.var([xi[i] for xi in xi_6_2],ddof=1)) )

    xi_0bar_3.append( np.mean([xi[i] for xi in xi_0_3]) )
    xi_2bar_3.append( np.mean([xi[i] for xi in xi_2_3]) )
    xi_4bar_3.append( np.mean([xi[i] for xi in xi_4_3]) )
    xi_6bar_3.append( np.mean([xi[i] for xi in xi_6_3]) )

    xi_0var_3.append( np.sqrt(np.var([xi[i] for xi in xi_0_3],ddof=1)) )
    xi_2var_3.append( np.sqrt(np.var([xi[i] for xi in xi_2_3],ddof=1)) )
    xi_4var_3.append( np.sqrt(np.var([xi[i] for xi in xi_4_3],ddof=1)) )
    xi_6var_3.append( np.sqrt(np.var([xi[i] for xi in xi_6_3],ddof=1)) )
    
    
    
nr = 2.5+np.linspace(5.,150.,30)[:-1]


#plt.errorbar(nr,(nr**2)*xi_0bar_rancross,yerr=(nr**2)*np.array(xi_0var_rancross),label=label)
#plt.errorbar(nr,(nr**2)*(-np.array(xi_2bar_rancross)),yerr=(nr**2)*np.array(xi_2var_rancross),label=label)
#plt.errorbar(nr,(nr**2)*xi_4bar_rancross,yerr=(nr**2)*np.array(xi_4var_rancross),label=label)
##plt.xscale('log')
#plt.legend()
#plt.title(r'$Nran=87^{3} = 658503$')
#plt.show()



#plt.plot(nr,xi_0var_1,label=r'$\sigma_{\xi_0(r)}\, Nran=87^{3}$',ls='-',color=seaborn.color_palette()[0])
#plt.plot(nr,xi_2var_1,label=r'$\sigma_{\xi_2(r)}\, Nran=87^{3}$',ls='-',color=seaborn.color_palette()[1])
plt.plot(nr,xi_4var_1,label=r'$\sigma_{\xi_4(r)}\, Nran=87^{3}$',ls='-',color=seaborn.color_palette()[2])

#plt.plot(nr,xi_0var_2,label=r'$\sigma_{\xi_0(r)}\, Nran=\sqrt{2}87^{3}$',ls='--',color=seaborn.color_palette()[0])
#plt.plot(nr,xi_2var_2,label=r'$\sigma_{\xi_2(r)}\, Nran=\sqrt{2}87^{3}$',ls='--',color=seaborn.color_palette()[1])
plt.plot(nr,xi_4var_2,label=r'$\sigma_{\xi_4(r)}\, Nran=\sqrt{2}87^{3}$',ls='--',color=seaborn.color_palette()[2])

#plt.plot(nr,xi_0var_3,label=r'$\sigma_{\xi_0(r)}\, Nran=87^{3} cross$',ls='-.',color=seaborn.color_palette()[0])
#plt.plot(nr,xi_2var_3,label=r'$\sigma_{\xi_2(r)}\, Nran=87^{3} cross$',ls='-.',color=seaborn.color_palette()[1])
plt.plot(nr,xi_4var_3,label=r'$\sigma_{\xi_4(r)}\, Nran=87^{3} cross$',ls='-.',color=seaborn.color_palette()[2])

plt.legend()
plt.show()
