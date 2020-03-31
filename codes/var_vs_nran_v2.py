"""
Variances vs nran for a specific distance 

v2: Split the code into calculation and plotting. No need to recalculate everytime I want to see a plot

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import seaborn
import scipy

names=['xi_0','xi_2','xi_4','xi_6']

nr = 2.5+np.linspace(5.,150.,30)[:-1]

id_r = np.where(nr==102.5)[0][0]

mps=[0,2,4]

#Cantidad de archivos que tengo
Nglass = 200
Nrancross = 400
Nran = 200
Nransplit = 400

xi_l_ran = []
xi_l_rancross = []
xi_l_glass = []
xi_l_ransplit = []

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

for nran in [87**3,2*87**3,4*87**3,8*87**3]:
    for i in range(Nran):
        xi_l_ran.append( ascii.read('../data/out/ran/xi_l_{}_{}.txt'.format(nran,i),names=names))
    for i in range(0,Nrancross,2):
        xi_l_rancross.append( ascii.read('../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(nran,i,i+1),names=names))
        xi_l_ransplit.append( ascii.read('../data/out/ransplit/xi_l_ransplit_{}_{}-{}.txt'.format(nran,i,i+1),names=names))

        

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


    ran0.append([nran,xi_0var_ran[id_r]])
    ran2.append([nran,xi_2var_ran[id_r]])
    ran4.append([nran,xi_4var_ran[id_r]])
        
    rancross0.append([nran,xi_0var_rancross[id_r]])
    rancross2.append([nran,xi_2var_rancross[id_r]])
    rancross4.append([nran,xi_4var_rancross[id_r]])

    ransplit0.append([nran,xi_0var_ransplit[id_r]])
    ransplit2.append([nran,xi_2var_ransplit[id_r]])
    ransplit4.append([nran,xi_4var_ransplit[id_r]])

        



for nran in [87**3,2*87**3,4*87**3,8*87**3]:
    for i in range(0,Nglass,2):
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


    glass0.append([nran,xi_0var_glass[id_r]])
    glass2.append([nran,xi_2var_glass[id_r]])
    glass4.append([nran,xi_4var_glass[id_r]])


ran0 = np.array(ran0)
ran2 = np.array(ran2)
ran4 = np.array(ran4)

rancross0 = np.array(rancross0)
rancross2 = np.array(rancross2)
rancross4 = np.array(rancross4)

glass0 = np.array(glass0)
glass2 = np.array(glass2)
glass4 = np.array(glass4)

ransplit0 = np.array(ransplit0)
ransplit2 = np.array(ransplit2)
ransplit4 = np.array(ransplit4)

#S o m e   P a r a m e t e r s
fs=20
lfs=14
red=seaborn.color_palette()[3]
blue=seaborn.color_palette()[0]
green=seaborn.color_palette()[2]

ax2 = plt.subplot(312)
plt.plot(ran2[:,0],ran2[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran$')
plt.plot(rancross2[:,0],rancross2[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross$')
plt.plot(glass2[:,0],glass2[:,1],color=red,marker='o',linestyle='-',label=r'$Glass$')
plt.plot(ransplit2[:,0],ransplit2[:,1],color='k',marker='x',linestyle=':',label=r'$Split$')
plt.ylabel(r'$\sigma_{\xi_2(r)}$',fontsize=fs)
plt.xscale('log')
plt.yscale('log')
box = ax2.get_position()
#ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=lfs)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)


 
ax1 = plt.subplot(311, sharex=ax2, sharey=ax2)
plt.plot(ran0[:,0],ran0[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran$')
plt.plot(rancross0[:,0],rancross0[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross$')
plt.plot(glass0[:,0],glass0[:,1],color=red,marker='o',linestyle='-',label=r'$Glass$')
plt.plot(ransplit0[:,0],ransplit0[:,1],color='k',marker='x',linestyle=':',label=r'$Split$')
plt.ylabel(r'$\sigma_{\xi_0(r)}$',fontsize=fs)
plt.xscale('log')
plt.yscale('log')
box = ax1.get_position()
#ax1.set_position([box.x0, box.y0, box.width, box.height*.8])
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),ncol=3, fancybox=True, shadow=True,fontsize=lfs)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)

ax3 = plt.subplot(313, sharex=ax2, sharey=ax2)
plt.plot(ran4[:,0],ran4[:,1],color=blue,marker='^',linestyle='--',label=r'$Ran$')
plt.plot(rancross4[:,0],rancross4[:,1],color=green,marker='s',linestyle='-.',label=r'$RanCross$')
plt.plot(glass4[:,0],glass4[:,1],color=red,marker='o',linestyle='-',label=r'$Glass$')
plt.plot(ransplit4[:,0],ransplit4[:,1],color='k',marker='x',linestyle=':',label=r'$Split$')
plt.ylabel(r'$\sigma_{\xi_4(r)}$',fontsize=fs)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$N_{ran}$',fontsize=fs)
box = ax3.get_position()
#ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=lfs)

plt.suptitle(r'$R = {}Mpc$'.format(nr[19]))
plt.subplots_adjust(hspace=0)
plt.show()

#-----------------------------------------------------


# S L O P E S

glass_m=[]
for data in [glass0,glass2,glass4]:
	x = data[:,0]
	y = data[:,1]
	m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
	glass_m.append(m)	
	#print 'Glass:',m
	#ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    

ran_m=[]
for data in [ran0,ran2,ran4]:
	x = data[:,0]
	y = data[:,1]
	m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
	ran_m.append(m)	
	#print 'Glass:',m
	#ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    

rancross_m=[]
for data in [rancross0,rancross2,rancross4]:
	x = data[:,0]
	y = data[:,1]
	m,b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(x),np.log10(y))
	rancross_m.append(m)	
	#print 'Glass:',m
	#ax.plot(x,10.**(np.log10(x)*m+b),color=seaborn.color_palette()[3])    



plt.scatter(range(len(glass_m)),glass_m,marker='^',color=red)
plt.scatter(range(len(ran_m)),ran_m,marker='o',color=blue)
plt.scatter(range(len(rancross_m)),rancross_m,marker='^',color=green)

plt.show()



