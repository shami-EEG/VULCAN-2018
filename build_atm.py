import numpy as np
from numpy import polynomial
import scipy
from scipy import interpolate
import scipy.optimize as sop
import subprocess
import pickle

import vulcan_cfg
from phy_const import kb, Navo, r_sun, au
from vulcan_cfg import nz
import chem_funs
from chem_funs import ni, nr  # number of species and reactions in the network
species = chem_funs.spec_list

###
with open(vulcan_cfg.com_file, 'r') as f:
    columns = f.readline() # reading in the first line
    num_ele = len(columns.split())-2 # number of elements (-2 for removing "species" and "mass") 
type_list = ['int' for i in range(num_ele)]
type_list.insert(0,'U20'); type_list.append('float')
compo = np.genfromtxt(vulcan_cfg.com_file,names=True,dtype=type_list)
# dtype=None in python 2.X but Sx -> Ux in python3
compo_row = list(compo['species'])
###


class InitialAbun(object):
    """
    Calculating the appropriate initial mixing ratios with the assigned elemental abundance
    """
    
    def __init__(self):
        self.ini_m = [0.9,0.1,0.,0.,0] # initial guess
        #self.EQ_ini_file = vulcan_cfg.EQ_ini_file
        
        self.atom_list = vulcan_cfg.atom_list

    def abun_lowT(self, x):
        """
        calculating the initial mixing ratios of the following 5 molecules (with CH4) 
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CH4 x4:He x5:NH3
        """
        O_H, C_H, He_H, N_H = vulcan_cfg.O_H, vulcan_cfg.C_H, vulcan_cfg.He_H, vulcan_cfg.N_H
        x1,x2,x3,x4,x5 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2 - (2*x1+2*x2+4*x3+3*x5)*O_H
        f3 = x3 - (2*x1+2*x2+4*x3+3*x5)*C_H
        f4 = x4 - (2*x1+2*x2+4*x3+3*x5)*He_H
        f5 = x5 - (2*x1+2*x2+4*x3+3*x5)*N_H
        return f1,f2,f3,f4,f5
        
    def abun_highT(self, x):
        """
        calculating the initial mixing ratios of the following 4 molecules (with CO)
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CO x4:He x5:N2
        """
        O_H, C_H, He_H, N_H = vulcan_cfg.O_H, vulcan_cfg.C_H, vulcan_cfg.He_H, vulcan_cfg.N_H
        x1,x2,x3,x4,x5 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2+x3 - (2*x1+2*x2)*O_H
        f3 = x3 - (2*x1+2*x2)*C_H
        f4 = x4 - (2*x1+2*x2)*He_H
        f5 = x5*2 - (2*x1+2*x2)*N_H
        return f1,f2,f3,f4,f5
    
    
    def ini_mol(self):
        if vulcan_cfg.ini_mix == 'const_lowT':
            return np.array(sop.fsolve(self.abun_lowT, self.ini_m))
        # somehow const_highT is not stable at high P...
        # elif vulcan_cfg.ini_mix == 'const_highT':
            # return np.array(sop.fsolve(self.abun_highT, self.ini_m))
            
    def ini_fc(self, data_var, data_atm):
        # reading-in the default elemental abundances from Lodders 2009
        with open('fastchem_vulcan/chemistry/elements/element_abundances_lodders.dat' ,'r') as f:
            new_str = ""
            ele_list = list(vulcan_cfg.atom_list)
            ele_list.remove('H')
            
            if vulcan_cfg.use_solar == True: 
                new_str = f.read() # read in as a string
                print ("Initializing with the default solar abundance.")
                
            else: # using costomized elemental abundances
                print ("Initializing with the customized elemental abundance:")
                print ("{:4}".format('H') + str('1.'))
                for line in f.readlines():   
                        li = line.split()
                        if li[0] in ele_list:
                            sp = li[0].strip()
                            # read-in vulcan_cfg.sp_H
                            sp_abun = getattr(vulcan_cfg, sp+'_H')
                            fc_abun = 12. + np.log10(sp_abun)
                            line = sp + '\t' + "{0:.4f}".format(fc_abun) + '\n'
                            print ("{:4}".format(sp) + "{0:.4E}".format(sp_abun))
                        new_str += line
            
            # make the new elemental abundance file
            with open('fastchem_vulcan/chemistry/elements/element_abundances_vulcan.dat', 'w') as f: f.write(new_str)
            
        # write a T-P text file for fast_chem
        with open('fastchem_vulcan/input/vulcan_TP/vulcan_TP.dat' ,'w') as f:
            ost = ''   
            for n, p in enumerate(data_atm.pco): # p in bar in fast_chem
                ost +=  '{:.3e}'.format(p/1.e6) + '\t' + '{:.1f}'.format(data_atm.Tco[n])  + '\n'
            ost = ost[:-1]
            f.write(ost)
        
        subprocess.call(["./fastchem input/config.input"], shell=True, cwd='fastchem_vulcan/')
           
    def ini_y(self, data_var, data_atm): 
        # initial mixing ratios of the molecules
        
        ini_mol = self.ini_mol()  
        ini = np.zeros(ni)
        y_ini = data_var.y
        gas_tot = data_atm.M
      
        if vulcan_cfg.ini_mix == 'fc':
        
            self.ini_fc(data_var, data_atm)
            fc = np.genfromtxt('fastchem_vulcan/output/vulcan_EQ.dat', names=True, dtype=None, skip_header=0) 
            neutral_sp = [sp for sp in species if sp not in vulcan_cfg.excit_sp]

            for sp in neutral_sp:
                if sp in fc.dtype.names:
                    y_ini[:,species.index(sp)] = fc[sp]*gas_tot
                else: print (sp + ' not included in fastchem.')
            
            # remove the fc output
            subprocess.call(["rm vulcan_EQ.dat"], shell=True, cwd='fastchem_vulcan/output/')
        
        elif vulcan_cfg.ini_mix == 'fc_precal':
            
            pre_fc = 'fastchem_vulcan/output/vulcan_EQ_pre.dat'
            print ('\n Using the precalculated fastchem output: '+ pre_fc)
            
            fc = np.genfromtxt(pre_fc, names=True, dtype=None, skip_header=0)   
            for sp in species:
                y_ini[:,species.index(sp)] = fc[sp]*gas_tot
        
        elif vulcan_cfg.ini_mix == 'vulcan_ini':
            with open(vulcan_cfg.vul_ini, 'rb') as handle:
              vul_data = pickle.load(handle) 
            
            y_ini = np.copy(vul_data['variable']['y'])
            data_var.y = np.copy(y_ini)
        
        elif vulcan_cfg.ini_mix == 'const_mix':
            
            y_ini = np.zeros((nz,ni))
            for sp in vulcan_cfg.const_mix.keys():
                y_ini[:,species.index(sp)] = gas_tot* vulcan_cfg.const_mix[sp]
            data_var.y = np.copy(y_ini)
            
        else:
            
            for i in range(nz):
                
                if vulcan_cfg.ini_mix == 'const_lowT':
                    y_ini[i,:] = ini
                    y_ini[i,species.index('H2')] = ini_mol[0]*gas_tot[i]; y_ini[i,species.index('H2O')] = ini_mol[1]*gas_tot[i]; y_ini[i,species.index('CH4')] = ini_mol[2]*gas_tot[i]
                    y_ini[i,species.index('NH3')] = ini_mol[4]*gas_tot[i]
                    # assign rest of the particles to He
                    y_ini[i,species.index('He')] = gas_tot[i] - np.sum(y_ini[i,:])

                else:
                    raise IOError ('\nInitial mixing ratios unknow. Check the setting in vulcan_cfg.py.')
        
        ysum = np.sum(y_ini, axis=1).reshape((-1,1))
        # storing ymix
        data_var.ymix = y_ini/ysum
        
        data_var.y_ini = np.copy(y_ini)
        return data_var
        


    def ele_sum(self, data_var):
        
        for atom in self.atom_list:
            data_var.atom_ini[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = 0.
            data_var.atom_conden[atom] = 0.
            
        return data_var


class Atm(object):
    
    def __init__(self):
        self.g = vulcan_cfg.g # gravity
        self.P_b = vulcan_cfg.P_b
        self.P_t = vulcan_cfg.P_t
        self.type = vulcan_cfg.atm_type
        self.use_Kzz = vulcan_cfg.use_Kzz
        self.Kzz_prof = vulcan_cfg.Kzz_prof
        self.const_Kzz = vulcan_cfg.const_Kzz
        self.use_vz = vulcan_cfg.use_vz
        self.vz_prof = vulcan_cfg.vz_prof
        self.const_vz = vulcan_cfg.const_vz        
        
    def f_pico(self, data_atm):
        '''calculating the pressure at interface'''
        
        pco = data_atm.pco
        
        # construct pico
        pco_up1 = np.roll(pco,1)
        pi = (pco * pco_up1)**0.5
        pi[0] = pco[0]**1.5 * pco[1]**(-0.5)
        pi = np.append(pi,pco[-1]**1.5 * pco[-2]**(-0.5))
        
        # store pico
        data_atm.pico = pi
        #data_atm.pco = pco
             
        return data_atm
    
    
    def load_TPK(self, data_atm, output):
        
        PTK_fun = {}
        
        if self.type == 'isothermal': 
            data_atm.Tco = np.repeat(vulcan_cfg.Tiso,nz)
            data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            data_atm.vz = np.repeat(self.const_vz,nz-1)
            
        elif self.type == 'analytical': 
            
            # plotting T-P on the fly                               
            para_atm = vulcan_cfg.para_anaTP 
            
            # return the P-T function
            PTK_fun['pT'] = lambda pressure: self.TP_H14(pressure, *para_atm)        
            data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            data_atm.vz = np.repeat(self.const_vz,nz-1)
            
        elif self.type == 'file':
            
            if self.Kzz_prof == 'const':     
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file = atm_table['Pressure'], atm_table['Temp']
            
            elif self.Kzz_prof == 'file':
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file, Kzz_file = atm_table['Pressure'], atm_table['Temp'], atm_table['Kzz']
            
            else: raise IOError ('\n"Kzz_prof" (the type of Kzz profile) cannot be recongized.\nPlease assign it as "file" or "const" in vulcan_cfg.')

            if self.vz_prof == 'const': data_atm.vz = np.repeat(self.const_vz,nz-1)
            elif self.vz_prof == 'file': vz_file =  atm_table['vz']


            if max(p_file) < data_atm.pco[0] or min(p_file) > data_atm.pco[-1]:
                print ('Warning! P_b and P_t assgined in vulcan.cfg are out of range of the file.\nConstant extension will be used.')
            
            PTK_fun['pT'] = interpolate.interp1d(p_file, T_file, assume_sorted = False, bounds_error=False,\
             fill_value=(T_file[np.argmin(p_file)], T_file[np.argmax(p_file)] )  )
            # store Tco in data_atm
            try:
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            
            # for SciPy earlier than v0.18.0
            except ValueError:
                PTK_fun['pT'] = interpolate.interp1d(p_file, T_file, assume_sorted = False, bounds_error=False, fill_value=T_file[np.argmin(p_file)] )  
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            
            
            if self.use_Kzz == True and self.Kzz_prof == 'file': 
                
                PTK_fun['pK'] = interpolate.interp1d(p_file, Kzz_file, assume_sorted = False, bounds_error=False, fill_value=(Kzz_file[np.argmin(p_file)], Kzz_file[np.argmax(p_file)]) )
                # store Kzz in data_atm
                try:
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])
                # for SciPy earlier than v0.18.0
                except ValueError:
                    PTK_fun['pK'] = interpolate.interp1d(p_file, Kzz_file, assume_sorted = False, bounds_error=False, fill_value=Kzz_file[np.argmin(p_file)] )
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])
            
            elif self.Kzz_prof == 'const': data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
        
        else: raise IOError ('\n"atm_type" cannot be recongized.\nPlease trassign it in vulcan_cfg.')
                        
        if self.use_Kzz == False:
            # store Kzz in data_atm
            data_atm.Kzz = np.zeros(nz-1)
        if self.use_vz == False: 
            data_atm.vz = np.zeros(nz-1)   
               
        # calculating and storing M(the third body)
        data_atm.M = data_atm.pco/(kb*data_atm.Tco)
        data_atm.n_0 = data_atm.M.copy()
        
        # plot T-P profile
        if vulcan_cfg.plot_TP == True: output.plot_TP(data_atm)
        
        # print warning when T exceeds the valid range of Gibbs free energy (NASA polynomials)
        if np.any(np.logical_or(data_atm.Tco < 200, data_atm.Tco > 6000)): print ('Temperatures exceed the valid range of Gibbs free energy.\n')
        
        return data_atm
        
        
    # T(P) profile in Heng et al. 2014 (126)
    def TP_H14(self, pco, *args_analytical):
        
        # convert args_analytical tuple to a list so we can modify it
        T_int, T_irr, ka_0, ka_s, beta_s, beta_l = list(args_analytical) 
        
        g = vulcan_cfg.g
        P_b = vulcan_cfg.P_b 
     
        # albedo(beta_s) also affects T_irr
        albedo = (1.0-beta_s)/(1.0+beta_s)
        T_irr *= (1-albedo)**0.25    
        eps_L = 3./8; eps_L3=1./3; ka_CIA=0
        m = pco/g; m_0 = P_b/g; ka_l = ka_0 + ka_CIA*m/m_0
        term1 = T_int**4/4*(1/eps_L + m/(eps_L3*beta_l**2)*(ka_0 + ka_CIA*m/(2*m_0) ) )
        term2 = (1/(2*eps_L) + scipy.special.expn(2,ka_s*m/beta_s)*(ka_s/(ka_l*beta_s)- (ka_CIA)*m*beta_s/(eps_L3*ka_s*m_0*beta_l**2) ) )
        term3 = ka_0*beta_s/(eps_L3*ka_s*beta_l**2)*(1./3 - scipy.special.expn(4,ka_s*m/beta_s))
        term4 = 0. #related to CIA
        T = (term1 + T_irr**4/8*(term2 + term3 + term4) )**0.25
    
        return T
    
        
    def mol_mass(self, sp):
        ''' calculating the molar mass of each species?'''
        return compo['mass'][compo_row.index(sp)]

    def mean_mass(self, var, atm, ni):
        mu = np.zeros(nz)
        for i in range(ni):
            mu += self.mol_mass(species[i]) * var.ymix[:,i]
        atm.mu = mu
        return atm        
        
    def f_mu_dz(self, data_var, data_atm):  
            
        dz, zco = np.empty(nz), np.zeros(nz+1) # pressure defined at interfaces
        Tco, pico = data_atm.Tco.copy(), data_atm.pico.copy()
        g = self.g
        
        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)
        # updating and storing Hp
        data_atm.Hp = kb*Tco/(data_atm.mu/Navo*g) 

        for i in range(nz):
            dz[i] = data_atm.Hp[i] * np.log(pico[i]/pico[i+1])
            zco[i+1] = zco[i] + dz[i]

        zmco = 0.5*(zco + np.roll(zco,-1))
        zmco = zmco[:-1]
        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        # for the j grid, dzi[j] from the grid above and dz[j-1] from the grid below
        
        # updating and storing dz and dzi
        data_atm.dz = dz
        data_atm.dzi = dzi
        data_atm.zmco = zmco
        
        return data_atm
        
    def read_sflux(self, var, atm):
        
        atm.sflux_raw = np.genfromtxt(vulcan_cfg.sflux_file, dtype=float, skip_header=1, names = ['lambda','flux'])
        
        # for values outside the boundary => fill_value = 0
        inter_sflux = interpolate.interp1d(atm.sflux_raw['lambda'], atm.sflux_raw['flux'], bounds_error=False, fill_value=0)
        for n, ld in enumerate(var.bins):
            var.sflux_top[n] = inter_sflux(ld) * (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2 # 1/Pi for a half-hemisphere
            # not converting to actinic flux yet *1/(hc/ld)  
            
            # for TOA 1AU Earth
            #var.sflux_top[n] = inter_sflux(ld) * (1./vulcan_cfg.orbit_radius)**2 /(hc/ld)
    
    def mol_diff(self, atm):
        '''
        choosing the formulea of molecular diffusion for each species
        then constucting Dzz(z) 
        '''
        Tco = atm.Tco
        n_0 = atm.n_0 
        
        # using the value defined on the interface
        Tco_i = np.delete((Tco + np.roll(Tco,1))*0.5, 0)
        n0_i = np.delete((n_0 + np.roll(n_0,1))*0.5, 0)
        
        if vulcan_cfg.use_moldiff == False:
            for i in range(len(species)):
                # this is required even without molecular weight
                atm.ms[i] = compo[compo_row.index(species[i])][-1]
            return
        
        if vulcan_cfg.atm_base == 'H2':
            Dzz_gen = lambda T, n_tot, mi: 2.3E17*T**0.765/n_tot *( 16.04/mi*(mi+2.016)/18.059 )**0.5
        else: raise IOError ('\n Unknow atm_base!')
        
        for i in range(len(species)):
            # input should be float or in the form of nz-long 1D array
            atm.Dzz[:,i] = Dzz_gen(Tco_i, n0_i, self.mol_mass(species[i]))
            
            # constructing the molecular weight for every species
            # this is required even without molecular weight
            atm.ms[i] = compo[compo_row.index(species[i])][-1]
        
        # no exception needed!?    
        # exceptions
        # if vulcan_cfg.atm_base == 'H2':
#             atm.Dzz[:,species.index('H2')] = np.zeros(nz-1)
#         else: raise IOError ('\n Unknow atm_base!')
        
        
        # # thermal diffusion for H and H2
        # atm.alpha[species.index('H')] = -0.25
        # atm.alpha[species.index('H2')] = -0.25
    
    def BC_flux(self, atm):
        '''
        Reading-in the boundary conditions of constant flux (cm^-2 s^-1) at top/bottom
        '''
        # read in the const top BC
        if vulcan_cfg.use_topflux == True: 
            print ("Using the prescribed constant top flux.")
            with open (vulcan_cfg.top_BC_mix_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.top_flux[species.index(li[0])] = li[1]
        
        # read in the const bottom BC
        if vulcan_cfg.use_botflux == True: 
            print ("Using the prescribed constant bottom flux.")
            with open (vulcan_cfg.bot_BC_mix_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.bot_flux[species.index(li[0])] = li[1]
                        atm.bot_vdep[species.index(li[0])] = li[2]
                        
        # using fixed-mixing-ratio BC          
        if vulcan_cfg.use_fix_bot == True: 
            print ("Using the prescribed fixed bottom mixing ratios.")
            with open (vulcan_cfg.bot_BC_mix_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.bot_fix_sp[species.index(li[0])] = li[3]
         
        

if __name__ == "__main__":
    print("This module stores classes for constructing atmospheric structure \
    and initializing its chemical composition from the desinated elemental abudance.")