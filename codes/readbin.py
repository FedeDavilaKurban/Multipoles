f = open('glass_001.bin','rb')

for i in range(3):
    dtype=np.dtype( [('x','f8'),('y','f8'),('z','f8')  ] )
    f.seek(i)
    data=np.fromfile(f,dtype=dtype)
    print data[0],data[-1]
    print ' '

f = open('glass_001.bin','rb')

#dtype=np.dtype( [('N','i4')] )
#dtype=np.dtype( 'f4' )
dtype=np.dtype( [('x','f4'),('y','f4'),('z','f4')  ] )
f.seek(4)
data=np.fromfile(f,dtype=dtype)
print data[0]
print ' '
np.reshape(data,(87**3,3))
