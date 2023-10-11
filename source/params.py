# all physical and numerical parameters are set here.

# main parameters to set:
eta = 1e12                                    # Newtonian viscosity [Pa s]
k0 = 1e-17                                    # permeability pre-factor [m^2]
Ts = 272                                      # temperature at surface [K]
Tz0 = 25/1e3                                  # temperature gradient [K/m]
eta_i = 1e16                                  # ice viscosity [Pa s]
phi0 = 0.3                                    # initial porosity

discrete_lens = True

# physical parameters:
g = 9.81                                      # gravitational acceleration [m/s^2]
rho_i = 917                                   # ice density [kg/m^3]
rho_w = 1000                                  # density of water [kg/m^3]
rho_s = 2500                                  # density of sediment [kg/m^3]
mu = 1.8e-3                                   # water viscosity [Pa s]
zeta = eta/phi0                               # initial bulk viscosity [Pa s]
a,b = 6,0.5                                   # saturation and permeability exponents
Tm = 273.15                                   # reference temperature [K]
Lh = 3.34e5                                   # latent heat [J/kg]
g_iw = 0.034                                  # interfacial energy [J/m^2]
r_p = 1e-6                                    # pore throat radius [m]
Tf = (1-2*g_iw/(r_p*rho_i*Lh))*Tm             # infiltration temperature [K]
q_in = 0.07                                   # heat flux below fringe [W/m^2]
Ki = 2.1                                      # thermal conductivity of ice [kg/(m s^3 K)]
Kw = 0.56                                     # thermal conductivity of water [kg/(m s^3 K)]
Ks = 4.0                                      # thermal conductivity of sediment [kg/(m s^3 K)]

# scales 
N_sc = 2*g_iw/r_p
T_sc = Tm*N_sc/(rho_i*Lh)
z_sc = Ki*T_sc/q_in
w_sc = k0*N_sc/(mu*z_sc)
t_sc = z_sc/w_sc
stress_sc = ((4./3.)*eta + zeta)*w_sc/z_sc

# nondimensional parameters
G = rho_w*g*z_sc/N_sc
nu = rho_s/rho_w
gamma = 1e-1 #((4./3.)*eta + zeta)*w_sc/(N_sc*z_sc) # [sediment stess] / [N] 
alpha = 1e-1  #1e4/N_sc                             # plasticity [Pi]/[N]
beta = (4./3.)*eta_i*(w_sc/z_sc)/N_sc               # [ice stress] / [N]
Tz_sc = Tz0*z_sc/(Tm-Tf)                            # dimensionless temp. gradient

# domain parameters:
nz = 500                                        # Number of elements in z direction
zf = (Ts-Tf)/(Tz0*z_sc)                         # base of fringe position

# time-stepping parameters:
theta = 0.5                                     # time-integration parameter 
                                                # (0=forward euler, 0.5 = trapezoid, 1 = backward euler)

CFL = 0.25                                      # CFL constant for timestepping

# misc:
phi_min = 0.01                                  # minimum porosity constraint
phi_max = 0.99                                  # maximum porosity constraint

# print(gamma)
# print(beta)
# print(1e4/N_sc)