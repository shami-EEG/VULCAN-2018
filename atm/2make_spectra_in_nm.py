import numpy as np
import scipy

au = 1.4959787E14  # cm
r_sun = 6.957E10 # cm

new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'


with open('HD4539_coolcat.txt') as f:
    for line in f.readlines():
        if not line.startswith("#") and line.split(): 
            li = line.split()
            wl = float(li[0])*0.1
            flux = float(li[1])*10. *(663770.35*au/r_sun*0.75)**2
            
            if flux > 0:
                new_str += '{:<12}'.format(wl) + '\t' + "{:.2E}".format(flux) + '\n'

    
with open('HD4539.txt', 'w+') as f: f.write(new_str)   
