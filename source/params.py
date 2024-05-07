# all physical and numerical parameters are set here.

# main parameters to set:
k0 = 1.0e-15                                  # permeability pre-factor [m^2]
Ts = 272                                      # temperature at surface [K]
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
alpha = 1e3/N_sc                          # empirical stress scale N_0/[N]
Tz_sc = Tz0*z_sc/T_sc                     # dimensionless temperature gradient

# domain parameters:
nz = 1000                                 # Number of elements in z direction

# import numpy as np
# gamma = 1
# zeta = gamma*z_sc*N_sc/v_sc
# compact_length = np.sqrt((k0/mu)*zeta)

# print(compact_length)
# print(z_sc)
# print(np.log10(zeta))

