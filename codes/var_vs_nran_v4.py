"""
Variances vs nran for a specific distance 

v2: Split the code into calculation and plotting. No need to recalculate everytime I want to see a plot

v3: -Incluye un punto con las ccvt
    -Para varias escalas
    
v4: -En vez de calcular la media, usar la xi calculada sin randoms. Calcularla varianza respecto de esa xi !!!Parece que quise empezar a hacer esto pero no lo hice. Al final charlamos que era mejor hacerlo como estaba
    -Modifique el codigo para que calcule ZR y Glass al mismo tiempo que el resto de las distribuciones

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy

names=['xi_0','xi_2','xi_4','xi_6']

nr = 2.5+np.linspace(5.,150.,30)[:-1]

mps=[0,2,4]

#Cantidad de archivos que tengo
Nglass = 400
Nrancross = 400
Nran = 200
Nransplit = 400
Nccvt = 60
Nzr = 400

#S o m e   P a r a m e t e r s
fs=20
lfs=14
red=seaborn.color_palette()[3]
blue=seaborn.color_palette()[0]
green=seaborn.color_palette()[2]
other=seaborn.color_palette()[4]


#id_r = np.where(nr==102.5)[0][0]

#for id_r in [np.where(nr==12.5)[0][0],np.where(nr==102.5)[0][0],np.where(nr==142.5)[0][0]]:
#for id_r in [np.where(nr==12.5)[0][0]]:
for id_r in range(len(nr)):

    print('r =',nr[id_r])
    
    xi_l_ran = []
    xi_l_rancross = []
    xi_l_glass = []
    xi_l_ransplit = []
    xi_l_zelrec = []
    xi_l_ccvt = []

    ran0=[]
    ran2=[]
    ran4=[]
    rancross0=[]
    rancross2=[]
    rancross4=[]
    glass0=[]
    glass2=[]
    glass4=[]
    ransplit0=[]
    ransplit2=[]
    ransplit4=[]
    zelrec0=[]
    zelrec2=[]
    zelrec4=[]
    ccvt0=[]
    ccvt2=[]
    ccvt4=[]
	
    print('Ran, RanCross, RanSplit, Glass, and Zelrec')
    for nran in [87**3,2*87**3,4*87**3,8*87**3]:
        for i in range(Nran):
            xi_l_ran.append( ascii.read('../data/out/ran/xi_l_{}_{}.txt'.format(nran,i),names=names))
        for i in range(0,Nrancross,2):
            xi_l_rancross.append( ascii.read('../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
            xi_l_ransplit.append( ascii.read('../data/out/ransplit/xi_l_ransplit_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
        for i in range(0,Nglass,2):
            xi_l_glass.append( ascii.read('../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
        for i in range(0,Nzr,2):
            xi_l_zelrec.append( ascii.read('../data/out/zelrec/xi_l_zelrec_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

            

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

        xi_0_ransplit=[]
        xi_2_ransplit=[]
        xi_4_ransplit=[]
        xi_6_ransplit=[]
        for i in range(len(xi_l_rancross)):
            xi_0_ransplit.append( xi_l_ransplit[i]['xi_0'] )
            xi_2_ransplit.append( xi_l_ransplit[i]['xi_2'] )
            xi_4_ransplit.append( xi_l_ransplit[i]['xi_4'] )
            xi_6_ransplit.append( xi_l_ransplit[i]['xi_6'] )
            
        xi_0_glass=[]
        xi_2_glass=[]
        xi_4_glass=[]
        xi_6_glass=[]
        for i in range(len(xi_l_glass)):
            xi_0_glass.append( xi_l_glass[i]['xi_0'] )
            xi_2_glass.append( xi_l_glass[i]['xi_2'] )
            xi_4_glass.append( xi_l_glass[i]['xi_4'] )
            xi_6_glass.append( xi_l_glass[i]['xi_6'] )

        xi_0_zelrec=[]
        xi_2_zelrec=[]
        xi_4_zelrec=[]
        xi_6_zelrec=[]
        for i in range(len(xi_l_zelrec)):
            xi_0_zelrec.append( xi_l_zelrec[i]['xi_0'] )
            xi_2_zelrec.append( xi_l_zelrec[i]['xi_2'] )
            xi_4_zelrec.append( xi_l_zelrec[i]['xi_4'] )
            xi_6_zelrec.append( xi_l_zelrec[i]['xi_6'] )

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

        xi_0bar_ransplit=[]
        xi_2bar_ransplit=[]
        xi_4bar_ransplit=[]
        xi_6bar_ransplit=[]

        xi_0var_ransplit=[]
        xi_2var_ransplit=[]
        xi_4var_ransplit=[]
        xi_6var_ransplit=[]

        xi_0bar_glass=[]
        xi_2bar_glass=[]
        xi_4bar_glass=[]
        xi_6bar_glass=[]

        xi_0var_glass=[]
        xi_2var_glass=[]
        xi_4var_glass=[]
        xi_6var_glass=[]

        xi_0bar_zelrec=[]
        xi_2bar_zelrec=[]
        xi_4bar_zelrec=[]
        xi_6bar_zelrec=[]

        xi_0var_zelrec=[]
        xi_2var_zelrec=[]
        xi_4var_zelrec=[]
        xi_6var_zelrec=[]

        nbins_s = 29

        for i in range(nbins_s):
        #RAN
            xi_0bar_ran.append( np.mean([xi[i] for xi in xi_0_ran]) )
            xi_2bar_ran.append( np.mean([xi[i] for xi in xi_2_ran]) )
            xi_4bar_ran.append( np.mean([xi[i] for xi in xi_4_ran]) )
            xi_6bar_ran.append( np.mean([xi[i] for xi in xi_6_ran]) )

            xi_0var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_0_ran],ddof=1)) )
            xi_2var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_2_ran],ddof=1)) )
            xi_4var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_4_ran],ddof=1)) )
            xi_6var_ran.append( np.sqrt(np.var([xi[i] for xi in xi_6_ran],ddof=1)) )
            
        #RANCROSS
            xi_0bar_rancross.append( np.mean([xi[i] for xi in xi_0_rancross]) )
            xi_2bar_rancross.append( np.mean([xi[i] for xi in xi_2_rancross]) )
            xi_4bar_rancross.append( np.mean([xi[i] for xi in xi_4_rancross]) )
            xi_6bar_rancross.append( np.mean([xi[i] for xi in xi_6_rancross]) )

            xi_0var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_0_rancross],ddof=1)) )
            xi_2var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_2_rancross],ddof=1)) )
            xi_4var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_4_rancross],ddof=1)) )
            xi_6var_rancross.append( np.sqrt(np.var([xi[i] for xi in xi_6_rancross],ddof=1)) )        

        #RANSPLIT
            xi_0bar_ransplit.append( np.mean([xi[i] for xi in xi_0_ransplit]) )
            xi_2bar_ransplit.append( np.mean([xi[i] for xi in xi_2_ransplit]) )
            xi_4bar_ransplit.append( np.mean([xi[i] for xi in xi_4_ransplit]) )
            xi_6bar_ransplit.append( np.mean([xi[i] for xi in xi_6_ransplit]) )

            xi_0var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_0_ransplit],ddof=1)) )
            xi_2var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_2_ransplit],ddof=1)) )
            xi_4var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_4_ransplit],ddof=1)) )
            xi_6var_ransplit.append( np.sqrt(np.var([xi[i] for xi in xi_6_ransplit],ddof=1)) )        

		#GLASS
		
            xi_0bar_glass.append( np.mean([xi[i] for xi in xi_0_glass]) )
            xi_2bar_glass.append( np.mean([xi[i] for xi in xi_2_glass]) )
            xi_4bar_glass.append( np.mean([xi[i] for xi in xi_4_glass]) )
            xi_6bar_glass.append( np.mean([xi[i] for xi in xi_6_glass]) )

            xi_0var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_0_glass],ddof=1)) )
            xi_2var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_2_glass],ddof=1)) )
            xi_4var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_4_glass],ddof=1)) )
            xi_6var_glass.append( np.sqrt(np.var([xi[i] for xi in xi_6_glass],ddof=1)) )

		#ZELREC
		
            xi_0bar_zelrec.append( np.mean([xi[i] for xi in xi_0_zelrec]) )
            xi_2bar_zelrec.append( np.mean([xi[i] for xi in xi_2_zelrec]) )
            xi_4bar_zelrec.append( np.mean([xi[i] for xi in xi_4_zelrec]) )
            xi_6bar_zelrec.append( np.mean([xi[i] for xi in xi_6_zelrec]) )

            xi_0var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_0_zelrec],ddof=1)) )
            xi_2var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_2_zelrec],ddof=1)) )
            xi_4var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_4_zelrec],ddof=1)) )
            xi_6var_zelrec.append( np.sqrt(np.var([xi[i] for xi in xi_6_zelrec],ddof=1)) )

        ran0.append([nran,xi_0var_ran[id_r]])
        ran2.append([nran,xi_2var_ran[id_r]])
        ran4.append([nran,xi_4var_ran[id_r]])
            
        rancross0.append([nran,xi_0var_rancross[id_r]])
        rancross2.append([nran,xi_2var_rancross[id_r]])
        rancross4.append([nran,xi_4var_rancross[id_r]])

        ransplit0.append([nran,xi_0var_ransplit[id_r]])
        ransplit2.append([nran,xi_2var_ransplit[id_r]])
        ransplit4.append([nran,xi_4var_ransplit[id_r]])

        glass0.append([nran,xi_0var_glass[id_r]])
        glass2.append([nran,xi_2var_glass[id_r]])
        glass4.append([nran,xi_4var_glass[id_r]])
        
        zelrec0.append([nran,xi_0var_zelrec[id_r]])
        zelrec2.append([nran,xi_2var_zelrec[id_r]])
        zelrec4.append([nran,xi_4var_zelrec[id_r]])


    print('CCVT')

    for nran in [87**3]:
        for i in range(1,Nccvt+1,2):
            xi_l_ccvt.append( ascii.read('../data/out/ccvt/xi_l_ccvt_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

        xi_0_ccvt=[]
        xi_2_ccvt=[]
        xi_4_ccvt=[]
        xi_6_ccvt=[]
        for i in range(len(xi_l_ccvt)):
            xi_0_ccvt.append( xi_l_ccvt[i]['xi_0'] )
            xi_2_ccvt.append( xi_l_ccvt[i]['xi_2'] )
            xi_4_ccvt.append( xi_l_ccvt[i]['xi_4'] )
            xi_6_ccvt.append( xi_l_ccvt[i]['xi_6'] )
            
        xi_0bar_ccvt=[]
        xi_2bar_ccvt=[]
        xi_4bar_ccvt=[]
        xi_6bar_ccvt=[]

        xi_0var_ccvt=[]
        xi_2var_ccvt=[]
        xi_4var_ccvt=[]
        xi_6var_ccvt=[]

        for i in range(nbins_s):
            xi_0bar_ccvt.append( np.mean([xi[i] for xi in xi_0_ccvt]) )
            xi_2bar_ccvt.append( np.mean([xi[i] for xi in xi_2_ccvt]) )
            xi_4bar_ccvt.append( np.mean([xi[i] for xi in xi_4_ccvt]) )
            xi_6bar_ccvt.append( np.mean([xi[i] for xi in xi_6_ccvt]) )

            xi_0var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_0_ccvt],ddof=1)) )
            xi_2var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_2_ccvt],ddof=1)) )
            xi_4var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_4_ccvt],ddof=1)) )
            xi_6var_ccvt.append( np.sqrt(np.var([xi[i] for xi in xi_6_ccvt],ddof=1)) )


        ccvt0.append([nran,xi_0var_ccvt[id_r]])
        ccvt2.append([nran,xi_2var_ccvt[id_r]])
        ccvt4.append([nran,xi_4var_ccvt[id_r]])


    ran0 = np.array(ran0)
    ran2 = np.array(ran2)
    ran4 = np.array(ran4)

    rancross0 = np.array(rancross0)
    rancross2 = np.array(rancross2)
    rancross4 = np.array(rancross4)

    glass0 = np.array(glass0)
    glass2 = np.array(glass2)
    glass4 = np.array(glass4)

    ccvt0 = np.array(ccvt0)
    ccvt2 = np.array(ccvt2)
    ccvt4 = np.array(ccvt4)

    zelrec0 = np.array(zelrec0)
    zelrec2 = np.array(zelrec2)
    zelrec4 = np.array(zelrec4)

    ransplit0 = np.array(ransplit0)
    ransplit2 = np.array(ransplit2)
    ransplit4 = np.array(ransplit4)

    
    #P L O T S
    
    f=plt.figure(figsize=(6,10))
    ax1=f.add_subplot(311)#, sharex=ax2)
    ax1.plot(ran0[:,0],ran0[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran\,\, Nfiles={}$'.format(Nran))
    ax1.plot(rancross0[:,0],rancross0[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross\,\, Nfiles={}$'.format(Nrancross))
    ax1.plot(glass0[:,0],glass0[:,1],color=red,marker='o',linestyle='-',label=r'$Glass\,\, Nfiles={}$'.format(Nglass))
    ax1.plot(ransplit0[:,0],ransplit0[:,1],color='k',marker='x',linestyle=':',label=r'$Split\,\, Nfiles={}$'.format(Nransplit))
    ax1.plot(ccvt0[:,0],ccvt0[:,1],color=other,marker='s',linestyle='-.',label=r'$CCVT\,\, Nfiles={}$'.format(Nccvt))
    ax1.plot(zelrec0[:,0],zelrec0[:,1],color='m',marker='X',linestyle='-.',label=r'$Zelrec\,\, Nfiles={}$'.format(Nzr))
    ax1.set_ylabel(r'$\sigma_{\xi_0(r)}$',fontsize=fs)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width, box.height*.8])
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),ncol=3, fancybox=True, shadow=True,fontsize=lfs)
    ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax1.set_ylim(5E-5,9E-3)
    
    ax2 = f.add_subplot(312,sharex=ax1)
    ax2.plot(ran2[:,0],ran2[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran\,\, N={}$'.format(Nran))
    ax2.plot(rancross2[:,0],rancross2[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross Nfiles={}$'.format(Nrancross))
    ax2.plot(glass2[:,0],glass2[:,1],color=red,marker='o',linestyle='-',label=r'$Glass\,\, Nfiles={}$'.format(Nglass))
    ax2.plot(ccvt2[:,0],ccvt2[:,1],color=other,marker='s',linestyle='-.',label=r'$CCVT\,\, Nfiles={}$'.format(Nccvt))
    ax2.plot(ransplit2[:,0],ransplit2[:,1],color='k',marker='x',linestyle=':',label=r'$Split\,\, Nfiles={}$'.format(Nransplit))
    ax2.plot(zelrec2[:,0],zelrec2[:,1],color='m',marker='X',linestyle='-.',label=r'$Zelrec\,\, Nfiles={}$'.format(Nzr))
    ax2.set_ylabel(r'$\sigma_{\xi_2(r)}$',fontsize=fs)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax2.set_ylim(1E-4,7E-3)


    ax3 = f.add_subplot(313, sharex=ax2)
    ax3.plot(ran4[:,0],ran4[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran\,\, Nfiles={}$'.format(Nran))
    ax3.plot(rancross4[:,0],rancross4[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross$')
    ax3.plot(glass4[:,0],glass4[:,1],color=red,marker='o',linestyle='-',label=r'$Glass\,\, Nfiles={}$'.format(Nglass))
    ax3.plot(ransplit4[:,0],ransplit4[:,1],color='k',marker='x',linestyle=':',label=r'$Split\,\, Nfiles={}$'.format(Nransplit))
    ax3.plot(ccvt4[:,0],ccvt4[:,1],color=other,marker='s',linestyle='-.',label=r'$CCVT\,\, Nfiles={}$'.format(Nccvt))
    ax3.plot(zelrec4[:,0],zelrec4[:,1],color='m',marker='X',linestyle='-.',label=r'$Zelrec\,\, Nfiles={}$'.format(Nzr))
    ax3.set_ylabel(r'$\sigma_{\xi_4(r)}$',fontsize=fs)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('$N_{ran}$',fontsize=fs)
    ax3.set_ylim(1E-4,10E-3)

    f.suptitle(r'$R = {}Mpc$'.format(nr[id_r]),fontsize=15)
    f.subplots_adjust(hspace=0)
    f.tight_layout()
    f.savefig('variances_{}.png'.format(id_r),dpi=f.dpi)
    #plt.show()
    plt.close()
    #-----------------------------------------------------


    # S L O P E S

    #glass_m=[]
    #for data in [glass0,glass2,glass4]:
    #    x = data[:,0]
    #    y = data[:,1]
    #    m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
    #    glass_m.append(m)	
    #    #print 'Glass:',m
    #    #ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    

    #ran_m=[]
    #for data in [ran0,ran2,ran4]:
    #    x = data[:,0]
    #    y = data[:,1]
    #    m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
    #    ran_m.append(m)	
    #    #print 'Glass:',m
    #    #ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    

    #rancross_m=[]
    #for data in [rancross0,rancross2,rancross4]:
    #    x = data[:,0]
    #    y = data[:,1]
    #    m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
    #    rancross_m.append(m)	
    #    #print 'Glass:',m
    #    #ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    



    #plt.scatter(range(len(glass_m)),glass_m,marker='^',color=red)
    #plt.scatter(range(len(ran_m)),ran_m,marker='o',color=blue)
    #plt.scatter(range(len(rancross_m)),rancross_m,marker='^',color=green)

    #plt.show()

###
#SigmaGlass / SigmaZelrec
###


###


#-----------------------------------------------------
##SOME TESTS
#nran=1*87**3
#Nglass = 400
#xi_l_glass=[]
#for i in range(0,Nglass,2):
    #xi_l_glass.append( ascii.read('../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

##Agrupo los xi_l
#xi_0_glass=[]
#for i in range(len(xi_l_glass)):
    #xi_0_glass.append( xi_l_glass[i]['xi_0'] )

#for i in range(150,200):
    #plt.plot(nr, xi_l_glass[i]['xi_0'] )
    #plt.yscale('log')
#plt.show()

########################
#Nccvt = 20
#xi_l_ccvt=[]
#for i in range(1,Nccvt+1,2):
    #xi_l_ccvt.append( ascii.read('../data/out/ccvt/xi_l_ccvt_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

#xi_0_ccvt=[]
#for i in range(len(xi_l_ccvt)):
    #plt.loglog(nr, xi_l_ccvt[i]['xi_0'] )
#plt.show()


########################
#Nrancross = 400
#xi_l_rancross=[]
#for i in range(0,Nrancross,2):
    #xi_l_rancross.append( ascii.read('../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

#xi_0_rancross=[]
#for i in range(len(xi_l_rancross)):
    #plt.loglog(nr, xi_l_rancross[i]['xi_0'] )
#plt.show()
