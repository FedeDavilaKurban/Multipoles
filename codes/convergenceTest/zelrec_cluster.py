"""
Meant to be used inside of a multipole calculation 
"""

class MiniCat:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.we = np.ones(len(x))
        self.size = x.size
        #self.dist = dist
        #self.x = x
        #self.y = y
        #self.z = z
        self.newx = x*1
        self.newy = y*1
        self.newz = z*1 

def get_shift(x, y, z, f_x, f_y, f_z, nbins, binsize):

    xmin = 0.
    ymin = 0.
    zmin = 0.
    binsize1 = binsize
    nbins1 = nbins 
    
    xpos = (x-xmin)/binsize1
    ypos = (y-ymin)/binsize1
    zpos = (z-zmin)/binsize1
    
    i = xpos.astype(int)
    j = ypos.astype(int)
    k = zpos.astype(int)
    #print (k[np.where(k>=63)])
    ddx = xpos-i
    ddy = ypos-j
    ddz = zpos-k

    shift_x = np.zeros(x.size)
    shift_y = np.zeros(x.size)
    shift_z = np.zeros(x.size)

    for ii in range(2):
        for jj in range(2):
            for kk in range(2):

                weight = ( ((1-ddx)+ii*(-1+2*ddx))*\
                            ((1-ddy)+jj*(-1+2*ddy))*\
                            ((1-ddz)+kk*(-1+2*ddz))  )
                iii=i+ii
                jjj=j+jj
                kkk=k+kk
                #print(iii,jjj,kkk)
                
                #Para periodicidad
                iii[iii>=nbins1]-=nbins1
                jjj[jjj>=nbins1]-=nbins1
                kkk[kkk>=nbins1]-=nbins1
                
                pos = (iii, jjj, kkk)
                shift_x += f_x[pos]*weight
                shift_y += f_y[pos]*weight
                shift_z += f_z[pos]*weight

    return shift_x, shift_y, shift_z


box = bs

smooth_mult = 0.5

binsize = box/nbins
#smooth = binsize*smooth_mult

rhom_t = nran/((box)**3)
rhom = rhom_t*(binsize)**3
msepcel = rhom_t**(-1/3) #celda
msepesf = box/((nran*4*np.pi/3)**(1/3)) #esfera
msep=msepcel

smooth = msep*smooth_mult



beta = 0.
bias = 1.
f = 0.
#nthreads = 1#int(os.environ["OMP_NUM_THREADS"])
nthreads = int(os.environ["OMP_NUM_THREADS"])

verbose = 0

#write = 1

dat = MiniCat(data[:,0], data[:,1], data[:,2])

xmin=np.min(dat.x)
ymin=np.min(dat.y)
zmin=np.min(dat.z)

f=plt.figure()
plt.ylim([1E-9,10.])
plt.xlim([7E-3,3.])
x=np.geomspace(5E-3,5E-1,1000)
plt.loglog(x,(x*msepcel/(2*np.pi))**4,dashes=[3,1,2,2],c='b')
#plt.axvline(x=msepcel/(2.*np.pi)**4,label='msepcelda')


if nran==16**3:
    ccvt = ascii.read('../../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
    ccat = ArrayCatalog(ccvt,names=['x','y','z'])      
    ccat['Position'] = transform.StackColumns(ccat['x']*bs, ccat['y']*bs, ccat['z']*bs)
    c_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', Nmesh=512, BoxSize=bs)
    cr = FFTPower(c_mesh, mode='1d',dk=0.7/bs)
    cPk = cr.power
    plt.loglog(cPk['k'], cPk['power'].real*len(ccvt)/bs**3,label='CCVT',dashes=[5,2])
    plt.legend()
if nran==87**3:
    ccvt = ascii.read('../../data/ccvt_particle_87_capacity_10.txt',names=['x','y','z'])
    ccat = ArrayCatalog(ccvt,names=['x','y','z'])      
    ccat['Position'] = transform.StackColumns(ccat['x']*bs, ccat['y']*bs, ccat['z']*bs)
    c_mesh = ccat.to_mesh(compensated=True, resampler='tsc', position='Position', Nmesh=512, BoxSize=bs)
    cr = FFTPower(c_mesh, mode='1d',dk=0.7/bs)
    cPk = cr.power
    plt.loglog(cPk['k'], cPk['power'].real*len(ccvt)/bs**3,label='CCVT',dashes=[5,2])
    plt.legend()


print("Iterating...")
for iloop in range(niter):

    print("Loop %d" % iloop)
    #-- Creating arrays for FFTW
    if iloop==0:
        delta  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        deltak = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        rho    = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        rhok   = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_x  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_y  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')
        psi_z  = pyfftw.empty_aligned((nbins, nbins, nbins), dtype='complex128')

        #-- Initialize FFT objects and load wisdom if available
        wisdom_file = "wisdom."+str(nbins)+"."+str(nthreads)+'.npy'
        if os.path.isfile(wisdom_file) :
            print('Reading wisdom from ', wisdom_file)
            wisd = tuple(np.load(wisdom_file))
            print('Status of importing wisdom', pyfftw.import_wisdom(wisd))
        print('Creating FFTW objects...')
        fft_obj = pyfftw.FFTW(delta, delta, axes=[0, 1, 2], threads=nthreads)
        ifft_obj = pyfftw.FFTW(deltak, psi_x, axes=[0, 1, 2], \
                                threads=nthreads, \
                                direction='FFTW_BACKWARD')
        kr = fftfreq(nbins, d=binsize) * 2 * np.pi * smooth
        norm = np.exp(-0.5 * (  kr[:, None, None] ** 2 \
                                + kr[None, :, None] ** 2 \
                                + kr[None, None, :] ** 2))
    #-- Allocate galaxies and randoms to grid with CIC method
    #-- using new positions
    if verbose:
        print('Allocating galaxies in cells...')
    deltag = np.zeros((nbins, nbins, nbins), dtype='float64')
    fastmodules.allocate_gal_cic(deltag, dat.newx, dat.newy, dat.newz, dat.we, dat.size, xmin, ymin, zmin, box, nbins, 1)
    if verbose:
        print('Smoothing...')
    #deltag = gaussian_filter(deltag, smooth/binsize)
    ##-- Smoothing via FFTs
    rho = deltag + 0.0j
    fft_obj(input_array=rho, output_array=rhok)
    fastmodules.mult_norm(rhok, rhok, norm)
    ifft_obj(input_array=rhok, output_array=rho)
    deltag = rho.real
    #print(deltag)
    if verbose:
        print('Computing density fluctuations, delta...')
    # normalize using the randoms, avoiding possible divide-by-zero errors
    #fastmodules.normalize_delta_survey(delta, deltag, self.alpha, self.ran_min)
    #del(deltag)  # deltag no longer required anywhere
    delta[:]  = deltag/rhom - 1.
    
    del(deltag)
    if verbose:
        print('Fourier transforming delta field...')
    fft_obj(input_array=delta, output_array=delta)
    ## -- delta/k**2
    k = fftfreq(nbins, d=binsize) * 2 * np.pi
    fastmodules.divide_k2(delta, delta, k)

    if verbose:
        print('Inverse Fourier transforming to get psi...')
    fastmodules.mult_kx(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_x)
    fastmodules.mult_ky(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_y)
    fastmodules.mult_kz(deltak, delta, k, bias)
    ifft_obj(input_array=deltak, output_array=psi_z)

    if np.all(psi_x.real[0] == np.zeros(np.shape(psi_x.real[0]))): 
        print('No displacement field. Breaking.')
        break


    # from grid values of Psi_est = IFFT[-i k delta(k)/(b k^2)], compute the values at the galaxy positions
    if verbose:
        print('Calculating shifts...')
    #print('Shape of psi_x:',np.shape(psi_x))
    shift_x, shift_y, shift_z = get_shift(dat.newx,dat.newy,dat.newz,psi_x.real,psi_y.real,psi_z.real,nbins,binsize)
    
    #shiftmult=1.5
    dat.newx += shift_x*shiftmult
    dat.newy += shift_y*shiftmult
    dat.newz += shift_z*shiftmult
    
    for axis in [dat.newx, dat.newy, dat.newz]:
        axis[axis>box] -= box
        axis[axis<0.] += box

    if iloop==0 and not os.path.isfile(wisdom_file):
        wisd=pyfftw.export_wisdom()
        np.save(wisdom_file, wisd)
        print('Wisdom saved at', wisdom_file)
    
    #new = np.column_stack((dat.newx,dat.newy,dat.newz))
    zr1 = Table(np.column_stack((dat.newx,dat.newy,dat.newz)),names=['x','y','z']) 
    zcat1 = ArrayCatalog(zr1)
    zcat1['Position'] = transform.StackColumns(zcat1['x'], zcat1['y'], zcat1['z'])
    real_mesh = zcat1.to_mesh(compensated=True, resampler='tsc', position='Position', BoxSize=[box,box,box],Nmesh=256)
    r = FFTPower(real_mesh, mode='1d',dk=0.7/bs)
    Pkzelrec = r.power
    plt.loglog(Pkzelrec['k'], Pkzelrec['power']*nran/box**3)#,label='ZR1',color='red')
    #if iloop==0: plt.loglog(Pkzelrec['k'],10**8*(Pkzelrec['k']*msep/(2*np.pi))**4,dashes=[3,3,2,2],c='k') 
    f.savefig(pltfolder+'/{}'.format(iloop))

x=np.geomspace(5E-3,5E-1,1000)
plt.loglog(x,(x*msepcel/(2*np.pi))**4,dashes=[3,1,2,2],c='b')
if nran==16**3 or nran==87**3: plt.loglog(cPk['k'], cPk['power'].real*len(ccvt)/bs**3,c='b',dashes=[5,2])
f.savefig(pltfolder+'/{}'.format(iloop))

############################################################################################

##old = np.column_stack((dat.x,dat.y,dat.z))

#new = np.column_stack((dat.newx,dat.newy,dat.newz))

############################################################################################

