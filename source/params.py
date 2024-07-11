# all physical and numerical parameters are set here.

# main parameters to set:
k0 = 6.0e-14                                  # permeability pre-factor [m^2]
Tz0 = 25/1e3                                  # temperature gradient [K/m]

# porosity-stress relation empirical parameters
# (normal consolidation line)
phi0 = 0.445                                  # reference porosity at 1 kPa
e0 = phi0/(1-phi0)                            # reference void ratio at 1 kPa
d0 = 0.15                                     # stress parameter  

# physical parameters:
g = 9.81                                      # gravitational acceleration [m/s^2]
rho_i = 917                                   # ice density [kg/m^3]
rho_w = 1000                                  # density of water [kg/m^3]
rho_s = 2500                                  # density of sediment [kg/m^3]
mu = 1.8e-3                                   # water viscosity [Pa s]
a,b = 6,0.5                                   # saturation and permeability exponents
Tm = 273.15                                   # reference temperature [K]
Lh = 3.34e5                                   # latent heat [J/kg]
g_iw = 0.034                                  # interfacial energy [J/m^2]
r_p = 1e-6                                    # pore throat radius [m]
Tf = (1-2*g_iw/(r_p*rho_i*Lh))*Tm             # infiltration temperature [K]

# scales 
N_sc = 2*g_iw/r_p                             # effective stress scale
T_sc = Tm*N_sc/(rho_i*Lh)                     # temperature scale
z_sc = N_sc/((rho_s-rho_w)*g)                 # depth scale
v_sc = k0*N_sc/(mu*z_sc)                      # *** velocity scale
t_sc = z_sc/v_sc                              # *** time scale 

# nondimensional parameters
alpha = 1e3/N_sc                              # empirical stress scale N_0/[N]
Tz_sc = Tz0*z_sc/T_sc                         # dimensionless temperature gradient

# domain parameters:
nz = 1000                                     # Number of elements in z direction



## print some scales
## (example values for gamma, zeta, compaction length) 
# import numpy as np
# gamma = 4**2
# zeta = gamma*z_sc*N_sc/v_sc
# compact_length = ((k0/mu)*zeta)**0.5
# print('compaction length = '+str(np.round(compact_length,2))+' m')
# print('\n')
# print('With [k] = 10^'+str(str(np.round(np.log10(k0),1)))+' m^2:')
# print('[z] = '+str(z_sc)+' m')
# print('[t] = '+str(np.round(t_sc/3.154e7,3))+' yr')
# print('[v] = '+str(np.round(v_sc*3.154e7,3))+' m/yr')
# print('\n')
# print('example calculations for FIG 6:')
# print('v_i paper ~ 0.038*[v] --> '+str(np.round(3.8e-2*v_sc*3.154e7,3))+' m/yr')
# print('v_* paper ~ 0.2*v_i --> '+str(np.round(0.2*3.8e-2*v_sc*3.154e7,3))+' m/yr')
# print('t_max ~ 3000*[t] --> '+str(np.round(3e3*t_sc/3.154e7,2))+' yr')
# print('t_lens ~ t_max/18 --> '+str(np.round((1/18.)*3e3*t_sc/3.154e7,2))+' yr')
# print('max zeta --> '+'{:.2e}'.format(zeta)+' Pa s')
# print('\n')

# p_sc = 1e4
# z_l = 10
# eta_i = 5e14
# nu = (p_sc*z_sc / (2*eta_i*v_sc)) 
# print('nu = '+'{:.2e}'.format(nu))
