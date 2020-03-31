########################################################################################	
xi0_ran = [] # Monopole
xi2_ran = [] # Quadrupole
xi4_ran = [] # Hexadecapole
N = len(glob.glob('../../data/out/ran/*{}*'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_Ran is: {} (should be 200)'.format(N))
for iN in range(N):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_ran = ascii.read('../../data/out/ran/xi_l_{}_{}.txt'.format(Nran,iN),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_ran.append( xil_ran['xi_0'].data )
	xi2_ran.append( xil_ran['xi_2'].data )
	xi4_ran.append( xil_ran['xi_4'].data )

#Calculate Standard Deviation  
std_xi0_ran = []
std_xi2_ran = []	
std_xi4_ran = []		
for ir in range(len(nr)):
	std_xi0_ran.append( np.sqrt( np.var([item[ir] for item in xi0_ran],ddof=1) ))
	std_xi2_ran.append( np.sqrt( np.var([item[ir] for item in xi2_ran],ddof=1) ))
	std_xi4_ran.append( np.sqrt( np.var([item[ir] for item in xi4_ran],ddof=1) ))

std_ran = Table([std_xi0_ran,std_xi2_ran,std_xi4_ran,nr], names=('xi0', 'xi2','xi4','r'))

fname = '../../data/out/std_xi/std_ran_{}.txt'.format(Nran)
ascii.write(std_ran, fname, overwrite=True)
print('Created {}'.format(fname))
#########################################################################################	
	
