#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy
from astropy.io import ascii
import seaborn
from astropy.table import Table
from nbodykit.lab import *


def copybox(pos,box,buf):
    g1 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box and box-buf<i[2]<box])
    g1[:, 0] -= box
    g1[:, 1] -= box
    g1[:, 2] -= box
    g2 = np.asarray([i for i in pos if box-buf<i[1]<box and box-buf<i[2]<box])
    g2[:, 1] -= box
    g2[:, 2] -= box
    g3 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box and box-buf<i[2]<box])
    g3[:, 0] += box
    g3[:, 1] -= box
    g3[:, 2] -= box
    g4 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[2]<box])
    g4[:, 0] -= box
    g4[:, 2] -= box
    g5 = np.asarray([i for i in pos if box-buf<i[2]<box])
    g5[:, 2] -= box
    g6 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[2]<box])
    g6[:, 0] += box
    g6[:, 2] -= box
    g7 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf and box-buf<i[2]<box])
    g7[:, 0] -= box
    g7[:, 1] += box
    g7[:, 2] -= box

    g8 = np.asarray([i for i in pos if 0<i[1]<buf and box-buf<i[2]<box])
    g8[:, 1] += box
    g8[:, 2] -= box
    g9 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf and box-buf<i[2]<box])
    g9[:, 0] += box
    g9[:, 1] += box
    g9[:, 2] -= box
    g10 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box])
    g10[:, 0] -= box
    g10[:, 1] -= box
    g11 = np.asarray([i for i in pos if box-buf<i[1]<box])
    g11[:, 1] -= box
    g12 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box])
    g12[:, 0] += box
    g12[:, 1] -= box
    g13 = np.asarray([i for i in pos if box-buf<i[0]<box])
    g13[:, 0] -= box
    g14 = np.asarray([i for i in pos if 0<i[0]<buf])
    g14[:, 0] += box
    g15 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf])
    g15[:, 0] -= box
    g15[:, 1] += box
    g16 = np.asarray([i for i in pos if 0<i[1]<buf])
    g16[:, 1] += box
    g17 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf])
    g17[:, 0] += box
    g17[:, 1] += box
    g18 = np.asarray([i for i in pos if box-buf<i[0]<box and box-buf<i[1]<box and 0<i[2]<buf])
    g18[:, 0] -= box
    g18[:, 1] -= box
    g18[:, 2] += box
    g19 = np.asarray([i for i in pos if box-buf<i[1]<box and 0<i[2]<buf])
    g19[:, 1] -= box
    g19[:, 2] += box
    g20 = np.asarray([i for i in pos if 0<i[0]<buf and box-buf<i[1]<box and 0<i[2]<buf])
    g20[:, 0] += box
    g20[:, 1] -= box
    g20[:, 2] += box
    g21 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[2]<buf])
    g21[:, 0] -= box
    g21[:, 2] += box
    g22 = np.asarray([i for i in pos if 0<i[2]<buf])
    g22[:, 2] += box
    g23 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[2]<buf])
    g23[:, 0] += box
    g23[:, 2] += box
    g24 = np.asarray([i for i in pos if box-buf<i[0]<box and 0<i[1]<buf and 0<i[2]<buf])
    g24[:, 0] -= box
    g24[:, 1] += box
    g24[:, 2] += box
    g25 = np.asarray([i for i in pos if 0<i[1]<buf and 0<i[2]<buf])
    g25[:, 1] += box
    g25[:, 2] += box
    g26 = np.asarray([i for i in pos if 0<i[0]<buf and 0<i[1]<buf and 0<i[2]<buf])
    g26[:, 0] += box
    g26[:, 1] += box
    g26[:, 2] += box

    pos = np.vstack([pos, g1, g2, g3, g4, g5, g6, g7, g8, g9,
                    g10, g11, g12, g13, g14, g15, g16, g17,
                    g18, g19, g20, g21, g22, g23, g24, g25, g26])

    #plt.scatter(pos[:,0],pos[:,1])
    #plt.show()
    pos += buf
    return pos

glass = ascii.read('../data/glass_16.dat',names=['x','y','z'])
ccvt = ascii.read('../data/ccvt_particle_16_capacity_40.txt',names=['x','y','z'])
ran = Table(np.random.random((len(glass),3)),names=['x','y','z'])
zr = ascii.read('/home/fede/Proyectos/Voro/codes/zelrecon/zelrec_final.dat',names=['x','y','z'])

#%%
box = 1.
buf = 1.

pos = np.column_stack((glass['x'].data,glass['y'].data,glass['z'].data))
glassdata = copybox(pos,box,buf)

pos = np.column_stack((ccvt['x'].data,ccvt['y'].data,ccvt['z'].data))
ccvtdata = copybox(pos,box,buf)

pos = np.column_stack((ran['x'].data,ran['y'].data,ran['z'].data))
randata = copybox(pos,box,buf)

pos = np.column_stack((zr['x'].data,zr['y'].data,zr['z'].data))
zrdata = copybox(pos,box,buf)

#%%
boxsize = box+2*buf
g_tree = scipy.spatial.cKDTree(glassdata, boxsize=boxsize)
c_tree = scipy.spatial.cKDTree(ccvtdata, boxsize=boxsize)
r_tree = scipy.spatial.cKDTree(randata, boxsize=boxsize)
z_tree = scipy.spatial.cKDTree(zrdata, boxsize=boxsize)
#%%
g_s2 = []
c_s2 = []
r_s2 = []
z_s2 = []

np.random.seed(100000)
a, b = [buf,box+buf] #Maximum value for "r" will be 0.38 so 1-0.38=.62
x = (b-a)*np.random.random((1000,3))+a

meansep = (4.*np.pi*4096/3)**(-1./3)

rs = np.geomspace(meansep/10,2.,100)

for r in rs:
    idx = g_tree.query_ball_point(x,r)
    g_s2.append( np.var([len(i) for i in idx],ddof=1) )

    idx = c_tree.query_ball_point(x,r)
    c_s2.append( np.var([len(i) for i in idx],ddof=1) )

    idx = r_tree.query_ball_point(x,r)
    r_s2.append( np.var([len(i) for i in idx],ddof=1) )

    idx = z_tree.query_ball_point(x,r)
    z_s2.append( np.var([len(i) for i in idx],ddof=1) )

#%%
# g_s2 = []
# c_s2 = []
# r_s2 = []
# z_s2 = []

# np.random.seed(100000)
# a, b = [0.,.62] #Maximum value for "r" will be 0.38 so 1-0.38=.62
# x = (b-a)*np.random.random((1000,3))+a

# meansep = (4.*np.pi*4096/3)**(-1./3)

# rs = np.geomspace(meansep/10,10*meansep,100)

# for r in rs:
#     idx = g_tree.query_ball_point(x,r)
#     g_s2.append( np.var([len(i) for i in idx],ddof=1) )

#     idx = c_tree.query_ball_point(x,r)
#     c_s2.append( np.var([len(i) for i in idx],ddof=1) )

#     idx = r_tree.query_ball_point(x,r)
#     r_s2.append( np.var([len(i) for i in idx],ddof=1) )

#     idx = z_tree.query_ball_point(x,r)
#     z_s2.append( np.var([len(i) for i in idx],ddof=1) )

# #%%
# numPoints = 16**3

# meansep = (4.*np.pi*numPoints/3)**(-1./3)

# np.random.seed(100000)
# x = 3*np.random.random((100,3))

# ns = np.arange(1,6000,100)

# rs = []  #pre-allocate an array to store distances
# n = []
# for i in range(len(x)):
#     if i%10==0: print(i)
#     for k in ns:
#         point = x[i]
#         r, idx = g_tree.query(point,k=[k])
#         rs.append(r)
#         n.append(k)

# plt.xlabel('r')
# plt.ylabel('N(r)')
# plt.scatter(rs,n)

#%%
plt.figure(figsize=(8,6))
plt.loglog(rs,r_s2,label='Random')
plt.loglog(rs,g_s2,label='Glass')
plt.loglog(rs,c_s2,label='CCVT')
plt.loglog(rs,z_s2,label='Zeldovich')

plt.vlines(meansep/2,5E-4,10E3,ls='-.',label='Half MeanSep')

plt.loglog(rs,rs**(3)*2E4,ls=':',c='k',label=r'$R^{3}$')
#plt.loglog(rs,rs**(2)*1E3,ls='--',c='k',label=r'$R^{-4}$')

plt.legend()
plt.xlabel('R')
plt.ylabel(r'$\sigma^{2}$')

#plt.savefig('sigma_vs_r.png')
# %%
