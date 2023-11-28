# This file contains the functions needed for solving the compaction problem.
import numpy as np
from constitutive import D, Phi, dPhi, get_fields
from dolfinx.fem import (Constant, Function, FunctionSpace, dirichletbc,
                         locate_dofs_topological)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import create_interval, locate_entities_boundary
from dolfinx.nls.petsc import NewtonSolver
from mpi4py import MPI
from params import G, nu, nz
from petsc4py import PETSc
from post_process import interp
from ufl import Dx, Measure, SpatialCoordinate, TestFunction, ds


def weak_form(N,N_t,N_prev,v_i,domain,dt,eps,penalty,steady):
    # Weak form of the residual for the compaction problem
    x = SpatialCoordinate(domain)
    T,S,k = get_fields(x[0])
    dz = Measure("dx", metadata={"quadrature_degree": 10})

    phi = Phi(N) 

    v_i = Constant(domain, PETSc.ScalarType(v_i)) 
    dt = Constant(domain, PETSc.ScalarType(dt))
    pen = Constant(domain, PETSc.ScalarType(1e6))

    # define forcing function
    f = (1-phi)*(G*(1-phi)*(nu-1) + (1+T)*Dx(phi*S,0))*k/((1-phi*S)**2)

    # define effective stress below ice lens
    N_l = (1-phi*S)*(1+T)
    
    # define residual of weak form 
    F = dt*dPhi(N)*v_i*Dx(N,0)*N_t*dz
    F += dt*D(phi,S,eps)*Dx(N,0)*Dx(N_t,0)*dz - dt*Dx(f,0)*N_t*dz
  
    if steady == False:
        F += dPhi(N)*(N-N_prev)*N_t*dz 

    if penalty == True:
        F += pen*(N-N_l)*N_t*ds  

    return F

def solve_pde(domain,N_prev,N_f,v_i,dt,eps=1e-10,penalty = False,steady=False):
        # solves the momentum balance for a fixed porosity field
     
        # Define function space
        V = FunctionSpace(domain, ("CG", 1))
        N = Function(V)
        N_t = TestFunction(V)

        # define weak form
        F = weak_form(N,N_t,N_prev,v_i,domain,dt,eps,penalty=penalty,steady=steady)  

        x = SpatialCoordinate(domain)
        T,S,k = get_fields(x[0])
        z_,N_l = interp((1-Phi(N_prev)*S)*(1+T),domain)
        N_l = N_l[-1]

        # # effective stress beneath lens
        zl = domain.geometry.x[:,0].max()
        facets_l = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zl))
        dofs_l = locate_dofs_topological(V, domain.topology.dim-1, facets_l)
        bc_l = dirichletbc(PETSc.ScalarType(N_l), dofs_l,V)  
        
        # # effective stress at base 
        zb = domain.geometry.x[:,0].min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs_b = locate_dofs_topological(V, domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(N_f), dofs_b,V)      

        if penalty == False:
            bcs = [bc_b,bc_l]
     
        elif penalty == True:
            bcs = [bc_b]

        N.interpolate(N_prev)

        # Solve for phi
        problem = NonlinearProblem(F, N, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)

        solver.solve(N)
      
        return N

def time_stepping(domain,initial,N_f,v_i,timesteps,eps=1e-10,penalty=False):
    # solve the compaction problem given:
    # domain: the computational domain
    # initial: initial conditions 
    # N_f: effective stress at base of fringe
    # v_i: pulling speed 
    # timesteps: array of time steps
    # eps: diffusivity regularization parameter
    # penalty: flag for whether to enforce BC with penalty method
    # *see example.ipynb for an example of how to set these
    #
    # The solver returns:
    # N: effective stress (solution)
    # z: domain coordinates (change over time due ice lens motion)
    # new_lens: array of time steps when a new lens forms (new lens == 1, otherwise 0)

    dt = timesteps[1]-timesteps[0]
    nt = np.size(timesteps)

    # create arrays for saving solution
    N = np.zeros((nt,nz+1))
    z = np.zeros((nt,nz+1))
    new_lens = np.zeros(nt)

    V = FunctionSpace(domain, ("CG", 1)) 
    sol_n = Function(V)
    sol_n.interpolate(initial)

    # time-stepping loop
    for i in range(nt):

        if i>0:
            print('time step '+str(i+1)+' out of '+str(nt)+', N_min = '+'{:.2f}'.format(N_min)
                  +', ('+str(np.size(np.where(new_lens==1)))+' ice lenses)'+' \r',end='')
           
        # solve the compaction problem for sol = N
        sol = solve_pde(domain,sol_n,N_f,v_i,dt,eps,steady=False,penalty=penalty)

        # save the solution as numpy arrays
        z_i,N_i = interp(sol,domain)   
        N[i,:] = N_i
        z[i,:] = z_i

        # compute minimum effective stress
        N_min = N_i.min()

        if N_min <= 0.05: 
        #  initiate new lens position where effective stress reaches threshold value
           
            # set flag to 1 for lens initiation times 
            new_lens[i] = 1 

            # create new domain [z_b, z_n] where z_n corresponds to minimum
            # where the effective stress N reaches threshold          
            l = np.argmin(N_i)
            z_n = z_i[l]
            z_b = z_i.min()
            domain = create_interval(MPI.COMM_WORLD,nz,[z_b,z_n])

            # initialize solution on new domain        
            sol_n = initialize(domain,N_f)
    
            z_i,N_i = interp(sol_n,domain)     
            N[i,:] = N_i
            z[i,:] = z_i
         
        else:
            # calculate sediment velocity below ice lens
            v_s = get_vel(sol,v_i,domain)[-1]

            # evolve domain geometry according to velocity
            Z = domain.geometry.x
            z_l = Z[:,0].max()
            z_b = Z[:,0].min()
            Z[:,0] += dt*v_s*(Z[:,0]-z_b)/(z_l-z_b)

            V = FunctionSpace(domain, ("CG", 1)) 
            sol_n = Function(V)
            sol_n.interpolate(sol)
            sol_l = sol.x.array[-1]
            sol_n.x.array[-1] = sol_l

    return N, z, new_lens

def get_vel(N,v_i,domain):
    # calculate sediment velocity below ice lens
    phi = Phi(N)
    x = SpatialCoordinate(domain)
    T,S,k = get_fields(x[0])
    f = G*(1-phi)*(nu-1) + (1+T)*Dx(phi*S,0)
    v_s_expr = v_i - (Dx(N,0)+ f)*k/((1-phi*S)**2)
    z,v_s = interp(v_s_expr,domain)  
    return v_s

def initialize(domain,N_f):
    # initialize effective stress solution given a domain
    # and the boundary condition at the base
    # this assumes that solution is rigid (dN/dt + v_i dN/dz = 0)
    # so that N solves an elliptic pde 
    V = FunctionSpace(domain, ("CG", 1))
    z_l = domain.geometry.x[:,0].max()
    z_b = domain.geometry.x[:,0].min()
    sol_n = Function(V)
    sol_n.interpolate(lambda x: 0.001*(x[0]-z_b)/(z_l-z_b) + N_f)
    sol_n = solve_pde(domain,sol_n,N_f,0,1,eps=1e-2,penalty=True,steady=True)
    
    # decrease regularization parameter 
    eps_ = np.flipud(np.logspace(-5,-2,10))
    for j in range(eps_.size):
        sol_n = solve_pde(domain,sol_n,N_f,0,1,eps=eps_[j],penalty=False,steady=True)
    return sol_n
