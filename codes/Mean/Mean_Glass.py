########################################################################################	
xi0_g = [] # Monopole
xi2_g = [] # Quadrupole
xi4_g = [] # Hexadecapole
N_glass = len(glob.glob('../../data/out/glass/*glass*{}*'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_glass is: {} (should be 200)'.format(N_glass))
for iN in range(0,2*N_glass,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_g = ascii.read('../../data/out/glass/xi_l_glass_{}_{}-{}.txt'.format(Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_g.append( xil_g['xi_0'].data )
	xi2_g.append( xil_g['xi_2'].data )
	xi4_g.append( xil_g['xi_4'].data )

#Calculate Standard Deviation  
mean_xi0_g = []
mean_xi2_g = []	
mean_xi4_g = []		
for ir in range(len(nr)):
	mean_xi0_g.append( np.mean([item[ir] for item in xi0_g]) )
	mean_xi2_g.append( np.mean([item[ir] for item in xi2_g]) )
	mean_xi4_g.append( np.mean([item[ir] for item in xi4_g]) )

mean_g = Table([mean_xi0_g,mean_xi2_g,mean_xi4_g,nr], names=('xi0', 'xi2','xi4','r'))

fname_g = '../../data/out/mean_xi/mean_g_{}.txt'.format(Nran)
ascii.write(mean_g, fname_g, overwrite=True)
print('Created {}'.format(fname_g))
#########################################################################################	
	
