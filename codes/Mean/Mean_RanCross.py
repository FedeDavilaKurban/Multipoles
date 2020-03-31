########################################################################################	
xi0_rc = [] # Monopole
xi2_rc = [] # Quadrupole
xi4_rc = [] # Hexadecapole
N_rc = len(glob.glob('../../data/out/rancross/*rancross*{}*'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_RanCross is: {} (should be 200)'.format(N))
for iN in range(0,2*N_rc,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_rc = ascii.read('../../data/out/rancross/xi_l_rancross_{}_{}-{}.txt'.format(Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_rc.append( xil_rc['xi_0'].data )
	xi2_rc.append( xil_rc['xi_2'].data )
	xi4_rc.append( xil_rc['xi_4'].data )

#Calculate Standard Deviation  
mean_xi0_rc = []
mean_xi2_rc = []	
mean_xi4_rc = []		
for ir in range(len(nr)):
	mean_xi0_rc.append( np.mean([item[ir] for item in xi0_rc]) )
	mean_xi2_rc.append( np.mean([item[ir] for item in xi2_rc]) )
	mean_xi4_rc.append( np.mean([item[ir] for item in xi4_rc]) )

mean_rc = Table([mean_xi0_rc,mean_xi2_rc,mean_xi4_rc,nr], names=('xi0', 'xi2','xi4','r'))


fname_rc = '../../data/out/mean_xi/mean_rc_{}.txt'.format(Nran)
ascii.write(mean_rc, fname_rc, overwrite=True)
print('Created {}'.format(fname_rc))
#########################################################################################	
	
