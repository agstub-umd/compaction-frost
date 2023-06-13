# all physical and numerical parameters are set here.
import numpy as np

# main parameters to set:
phi0 = 0.5                                    # initial porosity 
eta = 1e10                                    # Newtonian viscosity (Pa s)
k0 = 1e-12                                    # permeability pre-factor (m^2)
T0 = 270                                      # temperature at surface (K)
G = 25/1e3                                    # temperature gradient (K/m)

# physical parameters:
g = 9.81                                      # gravitational acceleration (m/s^2)
rho_i = 917                                   # ice density (kg/m^3)
rho_w = 1000                                  # density of water (kg/m^3)
rho_s = 2500                                  # density of sediment (kg/m^3)
mu = 1.8e-3                                   # water viscosity (Pa s)
zeta = eta/phi0                               # initial bulk viscosity (Pa s)
a,b = 6,0.5                                   # permeability exponents
delta = np.sqrt((k0/mu)*((4./3.)*eta + zeta)) # compaction length (m)
Tm = 273.15                                   # reference temperature (K)

Lh = 3.34e5                                   # latent heat (J/kg)
g_iw = 0.034                                  # interfacial energy (J/m^2)
r_p = 1e-6                                    # pore throat radius (m)

Tf = (1-2*g_iw/(r_p*rho_i*Lh))*Tm             # infiltration temperature     

Ts = (Tf - T0)/(Tm-Tf)                        # scaled surface temp

# some scales for reference
w_scale = k0*rho_i*Lh*(1-Tf/Tm)/(mu*delta)
t_scale = delta/w_scale
stress_scale = w_scale*((4./3.)*eta + zeta)/delta 

# domain parameters:
L = 50                                        # length of domain in compaction lengths
nz = 1000                                     # Number of elements in z direction
dz = L/(nz+1)

# time-stepping parameters:
t_f = 10                                      # Final time (relative to time scale)
nt = 1000                                     # Number of time steps
dt = t_f/nt                                   # Timestep size
theta = 0.5                                   # time-integration parameter 
                                              # (0=forward euler, 0.5 = trapezoid, 1 = backward euler)

# misc:
phi_min = 1e-2                                # minimum porosity constraint

