"""
Variances vs nran for a specific distance 

"""
names=['xi_0','xi_2','xi_4','xi_6']

nr = 2.5+np.linspace(5.,150.,30)[:-1]

id_r = np.where(nr==102.5)[0][0]

mps=[0,2,4]


xi_l_ran = []
xi_l_rancross = []
xi_l_glass = []

for nran in [87**3,2*87**3,4*87**3,8*87**3]:
    for i in range(100):
        xi_l_ran.append( ascii.read('../data/out/ran/xi_l_{}_{}.txt'.format(nran,i),names=names))
    for i in range(0,200,2):
        xi_l_rancross.append( ascii.read('../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

        

    #Agrupo los xi_l
    xi_0_ran=[]
    xi_2_ran=[]
    xi_4_ran=[]
    xi_6_ran=[]
    for i in range(len(xi_l_ran)):
        xi_0_ran.append( xi_l_ran[i]['xi_0'] )
        xi_2_ran.append( xi_l_ran[i]['xi_2'] )
        xi_4_ran.append( xi_l_ran[i]['xi_4'] )
        xi_6_ran.append( xi_l_ran[i]['xi_6'] )

    xi_0_rancross=[]
    xi_2_rancross=[]
    xi_4_rancross=[]
    xi_6_rancross=[]
    for i in range(len(xi_l_rancross)):
        xi_0_rancross.append( xi_l_rancross[i]['xi_0'] )
        xi_2_rancross.append( xi_l_rancross[i]['xi_2'] )
        xi_4_rancross.append( xi_l_rancross[i]['xi_4'] )
        xi_6_rancross.append( xi_l_rancross[i]['xi_6'] )

    #Calculo las medias
    xi_0bar_ran=[]
    xi_2bar_ran=[]
    xi_4bar_ran=[]
    xi_6bar_ran=[]

    xi_0var_ran=[]
    xi_2var_ran=[]
    xi_4var_ran=[]
    xi_6var_ran=[]

    xi_0bar_rancross=[]
    xi_2bar_rancross=[]
    xi_4bar_rancross=[]
    xi_6bar_rancross=[]

    xi_0var_rancross=[]
    xi_2var_rancross=[]
    xi_4var_rancross=[]
    xi_6var_rancross=[]


    nbins_s = 29

    for i in range(nbins_s):
        xi_0bar_ran.append( np.mean([xi[i] for xi in xi_0_ran]) )
        xi_2bar_ran.append( np.mean([xi[i] for xi in xi_2_ran]) )
        xi_4bar_ran.append( np.mean([xi[i] for xi in xi_4_ran]) )
        xi_6bar_ran.append( np.mean([xi[i] for xi in xi_6_ran]) )

        xi_0var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_0_ran],ddof=1)) )
        xi_2var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_2_ran],ddof=1)) )
        xi_4var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_4_ran],ddof=1)) )
        xi_6var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_6_ran],ddof=1)) )
        
        xi_0bar_rancross.append( np.mean([xi[i] for xi in xi_0_rancross]) )
        xi_2bar_rancross.append( np.mean([xi[i] for xi in xi_2_rancross]) )
        xi_4bar_rancross.append( np.mean([xi[i] for xi in xi_4_rancross]) )
        xi_6bar_rancross.append( np.mean([xi[i] for xi in xi_6_rancross]) )

        xi_0var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_0_rancross],ddof=1)) )
        xi_2var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_2_rancross],ddof=1)) )
        xi_4var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_4_rancross],ddof=1)) )
        xi_6var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_6_rancross],ddof=1)) )        
        
        
    if 0 in mps:
        plt.scatter(xi_0var_ran[id_r],nran,label=r'$\sigma_{\xi_0(r)}\, Ran$' if nran==87**3 else "",marker='o',color=seaborn.color_palette()[0])
    if 2 in mps:
        plt.scatter(xi_2var_ran[id_r],nran,label=r'$\sigma_{\xi_2(r)}\, Ran$'if nran==87**3 else "",marker='x',color=seaborn.color_palette()[0])
    if 4 in mps:
        plt.scatter(xi_4var_ran[id_r],nran,label=r'$\sigma_{\xi_4(r)}\, Ran$'if nran==87**3 else "",marker='^',color=seaborn.color_palette()[0])

    if 0 in mps:
        plt.scatter(xi_0var_rancross[id_r],nran,label=r'$\sigma_{\xi_0(r)}\, RanCross$'if nran==87**3 else "",marker='o',color=seaborn.color_palette()[1])
    if 2 in mps:
        plt.scatter(xi_2var_rancross[id_r],nran,label=r'$\sigma_{\xi_2(r)}\, RanCross$'if nran==87**3 else "",marker='x',color=seaborn.color_palette()[1])
    if 4 in mps:
        plt.scatter(xi_4var_rancross[id_r],nran,label=r'$\sigma_{\xi_4(r)}\, RanCross$'if nran==87**3 else "",marker='^',color=seaborn.color_palette()[1])






for nran in [87**3,2*87**3]:
    for i in range(0,200,2):
        xi_l_glass.append( ascii.read('../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

    xi_0_glass=[]
    xi_2_glass=[]
    xi_4_glass=[]
    xi_6_glass=[]
    for i in range(len(xi_l_glass)):
        xi_0_glass.append( xi_l_glass[i]['xi_0'] )
        xi_2_glass.append( xi_l_glass[i]['xi_2'] )
        xi_4_glass.append( xi_l_glass[i]['xi_4'] )
        xi_6_glass.append( xi_l_glass[i]['xi_6'] )
        
    xi_0bar_glass=[]
    xi_2bar_glass=[]
    xi_4bar_glass=[]
    xi_6bar_glass=[]

    xi_0var_glass=[]
    xi_2var_glass=[]
    xi_4var_glass=[]
    xi_6var_glass=[]

    for i in range(nbins_s):
        xi_0bar_glass.append( np.mean([xi[i] for xi in xi_0_glass]) )
        xi_2bar_glass.append( np.mean([xi[i] for xi in xi_2_glass]) )
        xi_4bar_glass.append( np.mean([xi[i] for xi in xi_4_glass]) )
        xi_6bar_glass.append( np.mean([xi[i] for xi in xi_6_glass]) )

        xi_0var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_0_glass],ddof=1)) )
        xi_2var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_2_glass],ddof=1)) )
        xi_4var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_4_glass],ddof=1)) )
        xi_6var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_6_glass],ddof=1)) )




    if 0 in mps:
            plt.scatter(xi_0var_glass[id_r],nran,label=r'$\sigma_{\xi_0(r)}\, Glass$'if nran==87**3 else "",marker='o',color=seaborn.color_palette()[2])
    if 2 in mps:
            plt.scatter(xi_2var_glass[id_r],nran,label=r'$\sigma_{\xi_2(r)}\, Glass$'if nran==87**3 else "",marker='x',color=seaborn.color_palette()[2])
    if 4 in mps:
            plt.scatter(xi_4var_glass[id_r],nran,label=r'$\sigma_{\xi_4(r)}\, Glass$'if nran==87**3 else "",marker='^',color=seaborn.color_palette()[2])
    
    
    
    
plt.title(r'$R = {}Mpc$'.format(nr[19]))
plt.xlabel(r'$\sigma$')
plt.ylabel('Nran')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
