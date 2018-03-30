from pylab import *

testbeam = loadtxt("testbeam.dat")


close('all')
figure()
imshow(testbeam,interpolation='nearest')

show()
