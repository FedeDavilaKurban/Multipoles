########################################################################################	
zrdir = '/home/fede/mnt/clemente/HomDist/Multipoles/data/out/zelrec/'

xi0_zr = [] # Monopole
xi2_zr = [] # Quadrupole
xi4_zr = [] # Hexadecapole
N_zr = len(glob.glob(zrdir+'*zelrec*{}*_sm{}*'.format(Nran,sm))) # glob.glob works like "ls"
print('Number of xi_zelrec is: {} (should be 200)'.format(N_zr))
for iN in range(0,2*N_zr,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_zr = ascii.read(zrdir+'xi_l_zelrec_{}_{}-{}_sm{}.txt'.format(Nran,iN,iN+1,sm),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_zr.append( xil_zr['xi_0'].data )
	xi2_zr.append( xil_zr['xi_2'].data )
	xi4_zr.append( xil_zr['xi_4'].data )

#Calculate Standard Deviation  
std_xi0_zr = []
std_xi2_zr = []	
std_xi4_zr = []		
for ir in range(len(nr)):
	#print([item[ir] for item in xi0_zr])
	std_xi0_zr.append( np.std([item[ir] for item in xi0_zr],ddof=1) )
	std_xi2_zr.append( np.std([item[ir] for item in xi2_zr],ddof=1) )
	std_xi4_zr.append( np.std([item[ir] for item in xi4_zr],ddof=1) )

std_zr = Table([std_xi0_zr,std_xi2_zr,std_xi4_zr,nr], names=('xi0', 'xi2','xi4','r'))

fname_zr = '/home/fede/Proyectos/Multipoles/data/out/std_xi/std_zr_{}.txt'.format(Nran)
ascii.write(std_zr, fname_zr, overwrite=True)
print('Created {}'.format(fname_zr))
#########################################################################################	
	
