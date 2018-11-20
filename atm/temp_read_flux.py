import numpy as np
import scipy
#import vulcan_cfg
    
ost = ''

hd189 = np.genfromtxt('Kurucz_HD189.txt', names=True, dtype=None, skip_header=3)

au = 1.4959787E14  # cm
r_sun = 6.957E10 # cm
c_light = 2.99792E-2 # cm/s 

with open('HD189_flux_Kurucz.txt', 'w+') as f: 
    
    for n, wl in enumerate(hd189['WAVELENGTH']):
        # solar radius = 6.957E10 cm
        flux = hd189['FLUX'][n] * 1.e-4 * np.pi * (0.752*r_sun/(0.03142*au))**2 * (c_light/((wl*100.)**2))
        ost += "{:.6e}".format( wl * 1.e9) + '\t' +  "{:.6e}".format( flux ) + '\n'

    f.write(ost)
        