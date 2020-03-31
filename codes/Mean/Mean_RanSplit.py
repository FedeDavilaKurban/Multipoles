########################################################################################	
xi0_rs = [] # Monopole
xi2_rs = [] # Quadrupole
xi4_rs = [] # Hexadecapole
#N = len(glob.glob('../../data/out/ransplit/*ransplit*{}*'.format(Nran))) # glob.glob works like "ls"
N_rs = 200 # Hay muchos ransplit así que lo fijo a 200 como el resto, debo haber hecho de mas en alguna ocasión
print('Number of xi_ransplit is: {} (should be 200)'.format(N_rs))
for iN in range(0,2*N_rs,2):
	# "xil" meaning $\xi_l$, aka multipoles
	xil_rs = ascii.read('../../data/out/ransplit/xi_l_ransplit_{}_{}-{}.txt'.format(Nran,iN,iN+1),names=['xi_0','xi_2','xi_4','xi_6'])
	xi0_rs.append( xil_rs['xi_0'].data )
	xi2_rs.append( xil_rs['xi_2'].data )
	xi4_rs.append( xil_rs['xi_4'].data )

#Calculate Standard Deviation  
mean_xi0_rs = []
mean_xi2_rs = []	
mean_xi4_rs = []		
for ir in range(len(nr)):
	mean_xi0_rs.append( np.mean([item[ir] for item in xi0_rs]) )
	mean_xi2_rs.append( np.mean([item[ir] for item in xi2_rs]) )
	mean_xi4_rs.append( np.mean([item[ir] for item in xi4_rs]) )

mean_rs = Table([mean_xi0_rs,mean_xi2_rs,mean_xi4_rs,nr], names=('xi0', 'xi2','xi4','r'))

fname_rs = '../../data/out/mean_xi/mean_rs_{}.txt'.format(Nran)
ascii.write(mean_rs, fname_rs, overwrite=True)
print('Created {}'.format(fname_rs))
#########################################################################################	
	
