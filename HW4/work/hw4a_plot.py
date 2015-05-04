from pylab import *
from numpy import ma

x,y,u,v = np.genfromtxt('/home/student/Homeworks/HW4/work/velocity.out', delimiter=',', unpack=True)

print "size of x is: " + str(len(x)) 
print "size of y is: " + str(len(y)) 
print "size of u is: " + str(len(u)) 
print "size of v is: " + str(len(v)) 

Mag = sqrt(pow(u, 2) + pow(v, 2))

#u = u/Mag
#v = v/Mag

#x10 = x[::4]
#y10 = y[::4]
#u10 = u[::4]
#v10 = v[::4]

figure()
#subplot(111)
#for 
#arrow( x, y, u, v )

scatter(x,y,s=80, c=Mag, marker=(5,0))

title('particle velocity magnitude by color')

show()