import numpy as np
import scipy
#import vulcan_cfg
    
ost = ''
with open('H2 + hv -- H + H cross section.dat', 'r') as f:
    for line in f.readlines():
        li = line.split()
        for n in range(0, len(li), 2):
            ost +=  str(float(li[n])*0.1) + ', ' + "{:.4e}".format( float(li[n+1]) ) + '\n'
        
    
with open('h2_cross.txt', 'w+') as f: f.write(ost)
        