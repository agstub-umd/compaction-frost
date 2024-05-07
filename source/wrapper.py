from solvers import time_stepping, initialize
from dolfinx.mesh import create_interval
from mpi4py import MPI
from params import nz
import numpy as np

z_0 = np.arange(3.9,12.1,0.1)  
timesteps = np.linspace(0,1e3,1000)


def wrapper(i,v_i):
    N_f = 2.0  # effective stress at base of fringe 
    z_l0 = z_0[i]
    z_b = 1e-3
    gamma = 0
    N_c = 0.05
    const_phi = False
    domain = create_interval(MPI.COMM_WORLD,nz,[z_b,z_l0])
    initial = initialize(domain,N_f,gamma,const_phi,eps_min=1e-10)
    N, heave, visc, z, new_lens, converged = time_stepping(domain,initial,N_f,N_c,gamma,const_phi,v_i,timesteps,eps=1e-10)
    z_l = z[:,-1]
    return([z_l,new_lens,[converged]])