import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle
import numpy as np

# Setting the 2nd input argument as the filename of vulcan output   
vul_data = sys.argv[1]
# # Setting the 3rd input argument as the species names to be plotted (separated by ,)
# plot_spec = sys.argv[2]
# # Setting the 4th input argument as the output eps filename
# plot_name = sys.argv[3]
#
# plot_dir = vulcan_cfg.plot_dir
#
# # taking user input species and splitting into separate strings and then converting the list to a tuple
# plot_spec = tuple(plot_spec.split(','))
# nspec = len(plot_spec)

fc = np.genfromtxt('output/chem_output_kelt9b_inversion.dat', names=True, dtype=None, skip_header=0)


with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
  
vulcan_spec = data['variable']['species']

        
#sp_list = vulcan_spec
sp_list = ['H', 'H2', 'C', 'N', 'C2', 'O', 'H2O', 'OH', 'CO', 'CO2', 'CH4', 'C2H2', 'NH3', 'N2', 'HCN']

output_file = 'output/kelt9b_inversion_EQ.txt'

ost = 'T(K)   P(dyne/cm2)'
for sp in sp_list: ost += "{:>11s}".format(sp)
ost += "{:>11s}".format('Fe')
ost += "{:>11s}".format('Fe+')
ost += '\n'
#ost += 'mu\n'
for lev, tt in enumerate(data['atm']['Tco']):
    ost += "{:4.1f}".format(tt) + "{:>11.3E}".format(data['atm']['pco'][lev])
    for sp in sp_list:
        if sp=='CH4': ost += "{:>12.3E}".format(fc['C1H4'][lev])
        else:
            ost += "{:>12.3E}".format(data['variable']['y_ini'][lev, vulcan_spec.index(sp)]/data['atm']['n_0'][lev]) 
        #ost += "{:>12.3E}".format(data['variable']['ymix'][lev, vulcan_spec.index(sp)])
        
        #if sp=='CH4': print (data['variable']['y_ini'][lev, vulcan_spec.index(sp)]/data['atm']['n_0'][lev] / fc['C1H4'][lev])
        #if sp=='NH3': print (data['variable']['y_ini'][lev, vulcan_spec.index(sp)]/data['atm']['n_0'][lev] / fc['H3N1'][lev])
         
    ost += "{:>12.3E}".format(fc['Fe'][lev]) + "{:>12.3E}".format(fc['Fe1p'][lev])         
    ost += '\n'

ost = ost[:-1]
with open(output_file, "w") as of:
    of.write(ost)






