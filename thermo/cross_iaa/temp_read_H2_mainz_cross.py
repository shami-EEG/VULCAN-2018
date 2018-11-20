import numpy as np
import scipy
#import vulcan_cfg
    
ost = ''

h2 = np.genfromtxt('H2_mainz.txt',dtype=float, skip_header=1, names = ['lambda','cross'])

for n in range(len(h2['lambda'])):
    ost +=  str(float(h2['lambda'][n])) + ', ' + "{:.4e}".format( float(h2['cross'][n]) ) + '\n'
    
with open('H2_cross.txt', 'w+') as f: 
    f.write(ost)
        