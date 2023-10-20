import numpy as np

with open('monitor0001.dat','rb') as f:
    for k in xrange(4):
        data = np.fromfile(f, dtype=np.float32, count = 2*3)
        print np.reshape(data,(2,3))