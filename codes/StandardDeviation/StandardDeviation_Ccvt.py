########################################################################################	
xi0_c = [] # Monopole
xi2_c = [] # Quadrupole
xi4_c = [] # Hexadecapole
N_c = len(glob.glob('../../data/out/ccvt/*ccvt*{}*'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_ccvt is: {} (should be 60)'.format(N_c))
for iN in range(1,2*N_c,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_c = ascii.read('../../data/out/ccvt/xi_l_ccvt_{}_{}-{}.txt'.format(Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_c.append( xil_c['xi_0'].data )
	xi2_c.append( xil_c['xi_2'].data )
	xi4_c.append( xil_c['xi_4'].data )

#Calculate Standard Deviation  
std_xi0_c = []
std_xi2_c = []	
std_xi4_c = []		
for ir in range(len(nr)):
	std_xi0_c.append( np.sqrt( np.var([item[ir] for item in xi0_c],ddof=1) ))
	std_xi2_c.append( np.sqrt( np.var([item[ir] for item in xi2_c],ddof=1) ))
	std_xi4_c.append( np.sqrt( np.var([item[ir] for item in xi4_c],ddof=1) ))

std_c = Table([std_xi0_c,std_xi2_c,std_xi4_c,nr], names=('xi0', 'xi2','xi4','r'))

fname_c = '../../data/out/std_xi/std_c_{}.txt'.format(Nran)
ascii.write(std_c, fname_c, overwrite=True)
print('Created {}'.format(fname_c))
#########################################################################################	
	
