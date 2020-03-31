########################################################################################	
xi0_zr = [] # Monopole
xi2_zr = [] # Quadrupole
xi4_zr = [] # Hexadecapole
N_zr = len(glob.glob('../../data/out/zelrec/*zelrec*{}*_clem.txt'.format(Nran))) # glob.glob works like "ls"
print('Number of xi_zelrec is: {} (should be 200)'.format(N_zr))
for iN in range(0,2*N_zr,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_zr = ascii.read('../../data/out/zelrec/xi_l_zelrec_{}_{}-{}_clem.txt'.format(Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_zr.append( xil_zr['xi_0'].data )
	xi2_zr.append( xil_zr['xi_2'].data )
	xi4_zr.append( xil_zr['xi_4'].data )

#Calculate Standard Deviation  
mean_xi0_zr = []
mean_xi2_zr = []	
mean_xi4_zr = []		
for ir in range(len(nr)):
	mean_xi0_zr.append( np.mean([item[ir] for item in xi0_zr]) )
	mean_xi2_zr.append( np.mean([item[ir] for item in xi2_zr]) )
	mean_xi4_zr.append( np.mean([item[ir] for item in xi4_zr]) )

mean_zr = Table([mean_xi0_zr,mean_xi2_zr,mean_xi4_zr,nr], names=('xi0', 'xi2','xi4','r'))

fname_zr = '../../data/out/mean_xi/mean_zr_{}.txt'.format(Nran)
ascii.write(mean_zr, fname_zr, overwrite=True)
print('Created {}'.format(fname_zr))
#########################################################################################	
	
