import numpy as np
import scipy
#import vulcan_cfg
    
ost = ''

ch3 = np.genfromtxt('wav_sigmaCH3.txt',dtype=float, skip_header=1, names = ['lambda','cross'])

for n in range(len(ch3['lambda'])):
    ost +=  str(float(ch3['lambda'][n]*0.1)) + ', ' + "{:.4e}".format( float(ch3['cross'][n]) ) + '\n'
    if ch3['cross'][n]>0:
        print (str(ch3['lambda'][n]) + '\t' + str(ch3['cross'][n]) )

with open('ch3_cross.txt', 'w+') as f: 
    f.write(ost)
        