import numpy as np
from scipy import interpolate
import timeit
import chem_funs
from phy_const import kb
import matplotlib.pyplot as plt
import matplotlib.legend as lg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False

# make a grid to compute EQ
# i.e. a single T-P convering from 1000K - 3000K and 100 bar to 1mbar  (size:20 x 50)
# and only need 2 or 3 different files for different C/O and metallicity

# read in fast_chem output as a table file    
fc = np.genfromtxt('output/TP_grid_solar.dat', names=True, dtype=None, skip_header=0)

p_fc = np.unique(fc['P'])
t_fc = np.unique(fc['T'])

def fEQ(sp,pres,temp): # faster
    # IMPORTANT: x:p, y:T because z is in shape of (y,x)
    '''
    pres in bar
    using log pres to interpolate 
    '''
    fEQ = interpolate.interp2d(np.log10(p_fc), t_fc, fc[sp].reshape((t_fc.size,p_fc.size)), kind='linear')  
    # IMPORTANT: Z must be in shape of (Y,X)
    
    return fEQ(np.log10(pres), temp).flatten()

def fEQ2(sp,temp,pres):
    '''
    pres in bar
    using log pres to interpolate 
    '''
    
    fEQ = interpolate.interpolate.RectBivariateSpline(t_fc, np.log10(p_fc), fc[sp].reshape(t_fc.size, p_fc.size), s=0)  # kind='cubic'
    
    return fEQ(temp, np.log10(pres)).reshape(-1) #.flatten()
    
def fEQ3(sp,temp,pres, use_2D_grids = False):
    # IMPORTANT: first T than P
    '''
    pres in bar
    using log pres to interpolate 
    '''
    
    fEQ = interpolate.interp2d(t_fc, np.log10(p_fc), fc[sp], kind='linear')  # IMPORTANT: first T than P
    
    if use_2D_grids == True:
        return fEQ(temp, np.log10(pres))
    else:
        return np.diag(fEQ(temp, np.log10(pres))).flatten()

# 3D interpolation for C/O or metallicity
# for sp in species:
#     fEQ_interp[sp] = RegularGridInterpolator((temp, pres, CtoO), data)  # size of x*y*z


###
# temp = np.arange(1000,2501,100)
# pres = np.logspace(2,-4,20)
# CtoO = 




# R1 H + H2O -> OH + H2 
def k1(T):
    return 7.50E-16 *T**1.600 *np.exp(-9720.0/T)   	

# R2 H + H2O <- OH + H2 
def k2(T):
    return k1(T)/ chem_funs.Gibbs(1,T)


# R3 O + H2 -> OH + H 
def k3(T):
    return 8.52E-20 *T**2.670 *np.exp(-3160.0/T )   	 

def k4(T):
    return k3(T)/ chem_funs.Gibbs(3,T)


# the reverse reaction of (9) in Cooper & Showman 2006 
# CH3O + M -> H + H2CO + M
def k9_inf(T):
    k = 1.5E11 *T * np.exp(-12880./T)
    return k
def k9_0(T):
    k = 1.4E-6 *T**-1.2 *np.exp(-7800./T)
    return k
    
def k9(p,T):
    M = p*1.E6/(kb*T)
    k = k9_0(T)*k9_inf(T) / (k9_0(T)*M + k9_inf(T) )
    # M is n in Cooper & Showman 2006
    return k

# R13 H + CH4 -> CH3 + H2
def k13(T):
    return 2.20E-20 *T**3.000 *np.exp(-4040.0/T)   	

def k14(T):
    return k13(T)/ chem_funs.Gibbs(13,T)

# R107 O + CH -> H + CO    
def k107(T):
    return 6.59E-11    

def k108(T):
    return k107(T)/ chem_funs.Gibbs(107,T)

# CH2 + OH -> H2CO + H    
def k153(T):
    return 3.01E-11
    
def k154(T):
    return k153(T)/ chem_funs.Gibbs(153,T)
    
    
def K17(T):
    K = np.exp( -( chem_funs.gibbs_sp('CO',T) + 3./2* chem_funs.gibbs_sp('H2',T) -  chem_funs.gibbs_sp('CH3O',T) ) )
    # * (kb*T/1.e6) ** (-3./2)
    return K
    
def tau_co(p,T):
    t = K17(T)/k9(p,T) *(p*1.E6)**-1 *(kb*T)**-0.5 *(1.E6)**1.5  *(kb*T/1.E6)**1.5
    return t
def tau_co_cs_old(p,T):
    t = K17(T)*kb*T /( k9(p,T)*p**1.5) /(p*1.E6)
    return t

# This reaction is rarely important!
# R293: OH + CH3 + M -> CH3OH + M
def k293_0(T):
    return 1.932E3*T**-9.88 *np.exp(-7544./T) + 5.109E-11*T**-6.25 *np.exp(-1433./T)
def k293_inf(T):
    return 1.031E-10 * T**-0.018 *np.exp(16.74/T)
def k293(p,T):
    M = p*1.E6/(kb*T)
    return k293_0(T)*k293_inf(T) / (k293_0(T)*M + k293_inf(T) )


def k294_0(T):
    return 3.32E-07 * np.exp(-34400./T)

def k294_inf(T):
    return 1.90E+16 * np.exp(-46200./T)

def k294_exp(p,T):
    M = p*1.E6/(kb*T)
    return k294_0(T)*k294_inf(T) / (k294_0(T)*M + k294_inf(T) )

# R294: reverse of R293
def k294(p,T):
    return k293(p,T)/ chem_funs.Gibbs(293,T)

# H + H + M -> H2 + M  
def k233_0(T):
    return 2.70E-31 *T**-0.600
def k233_inf(T):
    return 3.31E-06 *T**-1.00
    
def k233(p,T):
    M = p*1.E6/(kb*T)
    return k233_0(T) / (k233_0(T)*M/k233_inf(T) + 1. ) 

def k234(p,T):
    return k233(p,T)/ chem_funs.Gibbs(233,T)

# R125: CH2OH + H -> OH + CH3    
def k125(T):
    return 1.60E-10
    
# R63: H2 + CO2 -> CO + H2O
def k63(T):
    return 1.66E-15*T**0.5 *np.exp(-7550.0/T)
    
def k128(T):
    return k125(T)/ chem_funs.Gibbs(125,T)

# R115 	[  OH + C -> CO + H  ]
def k115(T):
    return 1.05E-12	* T**0.500

def k116(T):
    return k115(T)/ chem_funs.Gibbs(115,T)
    
# R125 	[  CH2OH + H -> OH + CH3  ]
def k125(T):
    return 1.60E-10
    
def k126(T):
    return k125(T)/ chem_funs.Gibbs(125,T)

# CH3OH + H -> CH3O + H2   
def k147(T):
    return 6.82E-20 *T**2.685 *np.exp(-4643./T)
    
def k148(T):
    return k147(T)/ chem_funs.Gibbs(147,T)

# CH2OH + M -> H + H2CO + M  
def k269(p,T):
    M = p*1.E6/(kb*T)
    k0 = 1.66E-10 *np.exp(-12630./T)
    kinf = 3.00E+09	*np.exp(-14600./T)
    return k0/(1+k0*M/kinf)

def k270(p,T):
    return k269(p,T)/ chem_funs.Gibbs(269,T)


# R141: CH3 + O -> H2CO + H 
def k141(T):
    return 1.40E-10
   
def k142(T):
    return k141(T)/ chem_funs.Gibbs(141,T)	

# CH2 + O -> CO + H + H     
def k151(T):    
    return 1.33E-10
    
def k152(T):
    return k151(T)/ chem_funs.Gibbs(151,T)	

# Moses: H2O + CH -> H2CO + H
def kM740(T):
    return 9.5E-12 *np.exp(380./T)



def Ka(T):
    K = np.exp( -( chem_funs.gibbs_sp('CH3OH',T) - 2* chem_funs.gibbs_sp('H2',T) -  chem_funs.gibbs_sp('CO',T) ) )
    return K
    
def Kb(T):
    K = np.exp( -( chem_funs.gibbs_sp('CH2OH',T) +  chem_funs.gibbs_sp('H',T) - 2* chem_funs.gibbs_sp('H2',T) -  chem_funs.gibbs_sp('CO',T) ) )
    return K

##################################

# # R141: CH3 + O -> H2CO + H
# def tau_co_tsai(p,T):
#     n = p*1.E6/(kb*T)
#     a = fEQ('CO',T,p)*n
#     b = k142(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2 + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2
#     #b = np.maximum(k142(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2, k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2)
#     return a/b
    
def tau_co_viss(p,T):
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    # fEQ is the mixing ratio
    a = fEQ('CO',T,p)*n
    b = k294(p,T)*fEQ('CH3OH',T,p)*n*M + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2
    #b = np.maximum(k294(p,T)*fEQ('CH3OH',T,p)*n*M , k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n)
    
    return a/b
    
def tau_co_CS(p,T):
    M = p*1.E6/(kb*T)
    a = kb*T/(p*1.e6)
    #b = np.maximum(k294(p,T)*Ka(T)*(p)**2 , k125(T)*Kb(T) *p )
    b =  k294(p,T)*Ka(T)*(p)**2 + k125(T)*Kb(T) *p  #
    return a/b

# maximum of Fig5 in Zahnle & Marley 2014
def tau_co_zah(p,T):
    t1 = 1/p *3.E-6 * np.exp(42000./T) 
    t2 = 1/p**2 *40. * np.exp(25000./T) 
    t3 = 1/p**2 * 3.E-4 * np.exp(36000./T) 
    return t3 
    
def tau_co_R14(p,T):
    n = p*1.E6/(kb*T)
    M = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k14(T)*fEQ('H2',T,p)*fEQ('CH3',T,p)*n**2 + k142(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2 + 0*(k154(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2 + k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2 )
    return a/b

def tau_co_R142_R154(p,T):
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k142(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2 + k154(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2
    return a/b

def tau_co_R270(p,T):
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2
    #b = k294(p,T)*fEQ('CH3OH',T,p)*M*n + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p) *n*n   
    return a/b
        
def tau_co_R125_R270(p,T):
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    a = fEQ('CH4',T,p)*n
    b = k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2 + k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2
    #b = k294(p,T)*fEQ('CH3OH',T,p)*M*n + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p) *n*n   
    return a/b
    
def tau_co_R142(p,T):
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k142(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2
    #b = k294(p,T)*fEQ('CH3OH',T,p)*M*n + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p) *n*n   
    return a/b
    
def tau_co_R148(p,T):
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k148(T)*fEQ('CH3O',T,p)*fEQ('H2',T,p)*n**2
    #b = k294(p,T)*fEQ('CH3OH',T,p)*M*n + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p) *n*n   
    return a/b
    
def tau_co_R125(p,T):
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    b = k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2
    return a/b
        
def tau_co_R125_R294(p,T):
    n = p*1.E6/(kb*T)
    M = n
    a = fEQ('CO',T,p)*n
    b = k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2 + k294(p,T)*fEQ('CH3OH',T,p)*M*n
    return a/b
    
def tau_co_R294(p,T):
    n = p*1.E6/(kb*T)
    M = n
    a = fEQ('CO',T,p)*n
    b = k294(p,T)*fEQ('CH3OH',T,p)*M*n #+ k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2 *0
    return a/b
    
def tau_co_R294_R125(p,T):
    n = p*1.E6/(kb*T)
    M = n
    a = fEQ('CO',T,p)*n
    b = k294(p,T)*fEQ('CH3OH',T,p)*M*n + k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2
    return a/b

def tau_co_R294_R125_270(p,T):
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    #b = k293(p,T)*fEQ('OH',T,p)*fEQ('CH3',T,p)*M*n**2  + k126(T)*fEQ('CH3',T,p)*fEQ('OH',T,p)*n**2
    b = k294(p,T)*fEQ('CH3OH',T,p)*M*n  + np.minimum( k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2, k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2)  
    return a/b

def tau_co_R142_R294_R125_270_R14_116(p,T):
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    a = fEQ('CO',T,p)*n
    #b = k293(p,T)*fEQ('OH',T,p)*fEQ('CH3',T,p)*M*n**2  + k126(T)*fEQ('CH3',T,p)*fEQ('OH',T,p)*n**2
    b = k294(p,T)*fEQ('CH3OH',T,p)*M*n  + np.minimum( k125(T)*fEQ('CH2OH',T,p)*fEQ('H',T,p)*n**2, k270(p,T)*fEQ('H',T,p)*fEQ('H2CO',T,p)*M*n**2)\
    + k142(T)*fEQ('H2CO',T,p)*fEQ('H',T,p)*n**2 + np.minimum( k14(T)*fEQ('H2',T,p)*fEQ('CH3',T,p)*n**2,  k116(T)*fEQ('CO',T,p)*fEQ('H',T,p)*n**2)
    return a/b

def tau_ch4(T, p):
    # CH4 timescale using R141_R293_R126_269_R13_115
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    a = fEQ('CH4',p,T)*n
    #b = k293(p,T)*fEQ('OH',p,T)*fEQ('CH3',p,T)*M*n**2  + k126(T)*fEQ('CH3',p,T)*fEQ('OH',p,T)*n**2
    b = k293(p,T)*fEQ('OH',p,T)*fEQ('CH3',p,T)*M*n**2  + np.minimum( k126(T)*fEQ('CH3',p,T)*fEQ('OH',p,T)*n**2, k269(p,T)*fEQ('CH2OH',p,T)*M*n)\
    + k141(T)*fEQ('CH3',p,T)*fEQ('O',p,T)*n**2 + np.minimum( k13(T)*fEQ('H',p,T)*fEQ('CH4',p,T)*n**2, k115(T)*fEQ('OH',p,T)*fEQ('C',p,T)*n**2)
    return a/b

def tau_ch4_contour(T, p):
    # CH4 timescale using R141_R293_R126_269_R13_115
   
    M = p*1.E6/(kb*T)
    n = p*1.E6/(kb*T)
    
    a = fEQ3('CH4',T,p)*n
    b = k293(p,T)*fEQ3('OH',T,p)*fEQ3('CH3',T,p)*M*n**2  + np.minimum( k126(T)*fEQ3('CH3',T,p)*fEQ3('OH',T,p)*n**2, k269(p,T)*fEQ3('CH2OH',T,p)*M*n)\
    + k141(T)*fEQ3('CH3',T,p)*fEQ3('O',T,p)*n**2 + np.minimum( k13(T)*fEQ3('H',T,p)*fEQ3('CH4',T,p)*n**2, k115(T)*fEQ3('OH',T,p)*fEQ3('C',T,p)*n**2)
    return a/b


plt.figure(0)
x = np.arange(500,2501,50)
y = np.logspace(-4,3,100)
Z_co, Z_ch4, Z_h = np.empty( (len(x),len(y) ) ), np.empty( (len(x),len(y) ) ), np.empty( (len(x),len(y) ) )
#X, Y = np.meshgrid(x, y)
#Z_co = tau_co_R142_R294_R125_270_R14_116(y,x)[0]
for n, tt in enumerate(x):
    Z_ch4[n,:] = tau_ch4(tt, y) 
    
    #Z_h[:,n] = tau_h_R233(y,tt) *0.1
    
    ###
    # *0.1 to fit data!!!
    ###
    
    #Z_co[:,n] = np.maximum(Z_co[:,n], Z_h[:,n])

CS = plt.contour(x,y,np.transpose(np.log10(Z_ch4)), np.arange(-2,21,2), colors='k')
#CS1 = plt.contour(x,y,np.transpose(np.log10(Z_ch4)), np.arange(21,24), colors='k')
plt.clabel(CS, [int(cs) for cs in CS.levels], fontsize=10, inline=True, inline_spacing=0, fmt='%i')

plt.ylim(ymin=1.E-4)
plt.title(r'log$_{10}$ $\tau$ $_{\mathrm{CH4}}$')
plt.gca().set_yscale('log')
plt.gca().invert_yaxis()
plt.xlabel('Temperature (K)' )
plt.ylabel('Pressure (bar)' )
plt.savefig('../../Plots/timescale_contour.png')
plot = Image.open("../../Plots/timescale_contour.png")
plot.show()
