import numpy as np
import scipy
#import vulcan_cfg
from scipy import interpolate

import matplotlib.pyplot as plt
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False    
            
sp = 'H2'
plot_name = sp + '_cross'


bins = np.arange(20.,400.1, 0.1) 
cross_raw = np.genfromtxt(sp+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross'])

# for values outside the boundary => fill_value = 0
inter_cross = interpolate.interp1d(cross_raw['lambda'], cross_raw['cross'], bounds_error=False, fill_value=0)

cross = np.empty(len(bins))
# for n, ld in enumerate(bins):
#     cross[n] = inter_cross(ld)

plt.figure()
plt.plot( cross_raw['lambda']*10., cross_raw['cross'] )
plt.gca().set_xscale('log')       
plt.gca().set_yscale('log')
plt.xlim((10.,2000.)) 
plt.ylim((1.E-35,1.E-18))
plt.ylabel("cross section (cm^2)")
plt.xlabel("wavelength (A)")
plt.title(sp)
plt.savefig( plot_name + '.png')
#plt.savefig( plot_name + '.eps')

plot = Image.open(plot_name + '.png')
plot.show()

