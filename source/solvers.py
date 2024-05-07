# This file contains the functions needed for solving the compaction problem.
import numpy as np
from constitutive import D, Phi, dPhi, get_fields,N_visc
from dolfinx.fem import (Constant,Expression, Function, FunctionSpace, dirichletbc,
                         locate_dofs_topological)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import create_interval, locate_entities_boundary
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.log import LogLevel, set_log_level
from mpi4py import MPI
from params import nz
from petsc4py import PETSc
from post_process import interp
from ufl import Dx, Measure, FiniteElement,SpatialCoordinate, TestFunction, ds,split,TestFunctions



def weak_form_N(N,N_t,N_n,v,gamma,const_phi,v_i,domain,dt,eps,penalty,elliptic):
    # Weak form of the residual for the plastic effective stress equation
    x = SpatialCoordinate(domain)
    T,S,k = get_fields(x[0])
    dz = Measure("dx", metadata={"quadrature_degree": 10})

    # wrap constants in dolfinx
    v_i = Constant(domain, PETSc.ScalarType(v_i)) 
    dt = Constant(domain, PETSc.ScalarType(dt))
    pen = Constant(domain, PETSc.ScalarType(1e10))

    # define porosity
    phi = Phi(N,const_phi=const_phi) 

    # define viscous compaction stress
    N_v = N_visc(v,phi,gamma) 

    # define forcing function   
    f = ((1-phi) + (1+T)*Dx(phi*S,0))*D(phi,S,eps)

    # define effective stress below ice lens
    N_l = (1-phi*S)*(1+T)
    
    # define residual of weak form 
    F = dt*dPhi(N,const_phi=const_phi)*v_i*Dx(N,0)*N_t*dz
    F += dt*D(phi,S,eps)*Dx(N,0)*Dx(N_t,0)*dz - dt*Dx(f,0)*N_t*dz

    if gamma>0:
        # append heave-dependent forcing function
        F += dt*D(phi,S,eps)*Dx(N_v,0)*Dx(N_t,0)*dz
  
    if elliptic == False:
        # add evolution term for dynamic simulations
        F += dPhi(N,const_phi=const_phi)*(N-N_n)*N_t*dz 

    if penalty == True:
        # penalty method used in initialize() function
        F += pen*(N-N_l)*N_t*ds  

    return F


def weak_form_v(v,v_t,N,domain,gamma,const_phi):
    dz = Measure("dx", metadata={"quadrature_degree": 10})
    x = SpatialCoordinate(domain)
    T,S,k = get_fields(x[0])
    phi = Phi(N,const_phi=const_phi)
    N_v = N_visc(v,phi,gamma)

    pen = Constant(domain, PETSc.ScalarType(1))


    if gamma>0:
        F = N_v*Dx(v_t,0)*dz + ((1-phi*S)**2 / k)*v*v_t*dz - ((1+T)*Dx(phi*S,0) + 1-phi + Dx(N,0))*v_t*dz
    else:
        F = ((1-phi*S)**2 / k)*v*v_t*dz - ((1+T)*Dx(phi*S,0) + 1-phi + Dx(N,0))*v_t*dz
    
    F += pen*Dx(v,0)*Dx(v_t,0)*ds
    return F

def solve_pde(domain,sol_n,N_f,gamma,const_phi,v_i,dt,eps=1e-10,penalty = False,elliptic=False):
        # solves the momentum balance for a fixed porosity field
     
        # Define function space
        P1 = FiniteElement('P',domain.ufl_cell(),1)     
        element = P1*P1
        V = FunctionSpace(domain,element)       

        sol = Function(V)
        (N,v) = split(sol)
        (N_t,v_t) = TestFunctions(V)

        # define weak form
        F_N = weak_form_N(N,N_t,sol_n.sub(0),v,gamma,const_phi,v_i,domain,dt,eps,penalty=penalty,elliptic=elliptic)  
        F_v = weak_form_v(v,v_t,N,domain,gamma,const_phi)

        F = F_N + F_v

        # define effective stress boundary conditions
        x = SpatialCoordinate(domain)
        T,S,k = get_fields(x[0])
        z_,N_l = interp((1-Phi(sol_n.sub(0),const_phi=const_phi)*S)*(1+T),domain)
        N_l = N_l[-1]

        # effective stress beneath lens
        zl = domain.geometry.x[:,0].max()
        facets_l = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zl))
        dofs_l = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_l)
        bc_l = dirichletbc(PETSc.ScalarType(N_l), dofs_l,V.sub(0))  
        
        # effective stress at base 
        zb = domain.geometry.x[:,0].min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs_b = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(N_f), dofs_b,V.sub(0))      

        if penalty == False:
            # supply both boundary conditions if not
            # using penalty method
            bcs = [bc_b,bc_l]
     
        elif penalty == True:
            # only set boundary condition at base strongly 
            # BC at lens is set weakly with a penalty method
            bcs = [bc_b]

        # set initial guess for Newton solver
        sol.sub(0).interpolate(sol_n.sub(0))
        sol.sub(1).interpolate(sol_n.sub(1))

        # solve for effective stress N
        problem = NonlinearProblem(F, sol, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
        solver.error_on_nonconvergence = False

        solver.convergence_criterion = "residual"
        solver.atol = 1e-12
        solver.rtol = 1e-12


        set_log_level(LogLevel.ERROR)

        # return number of iterations n and 
        # convergence flag
        n,converged = solver.solve(sol)

        if gamma == 0:
            f = heave_plastic(sol.sub(0),domain,const_phi)
            sol.sub(1).interpolate(Expression(f, V.sub(1).element.interpolation_points()))
        
      
        return sol,converged

def time_stepping(domain,initial,N_f,N_c,gamma,const_phi,v_i,timesteps,eps=1e-10):
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
    # v: heave rate (solution)
    # z: domain coordinates (change over time due ice lens motion)
    # new_lens: array of time steps when a new lens forms (new lens == 1, otherwise 0)


    dt = timesteps[1]-timesteps[0]
    nt = np.size(timesteps)

    # create arrays for saving solution
    N = np.zeros((nt,nz+1))
    visc = np.zeros((nt,nz+1))
    z = np.zeros((nt,nz+1))
    new_lens = np.zeros(nt)
    heave = np.zeros((nt,nz+1))

    # initial condition
    P1 = FiniteElement('P',domain.ufl_cell(),1)     
    element = P1*P1
    V = FunctionSpace(domain,element)               
    sol_n = Function(V)
    sol_n.sub(0).interpolate(initial.sub(0))
    sol_n.sub(1).interpolate(initial.sub(1))


    # time-stepping loop
    for i in range(nt):

        if i>0:
            print('time step '+str(i+1)+' out of '+str(nt)+', N_min = '+'{:.2f}'.format(N_min)
                  +', ('+str(np.size(np.where(new_lens==1)))+' ice lenses)'+' \r',end='')
             
        
        # solve the compaction problem for sol = N
        sol,converged = solve_pde(domain,sol_n,N_f,gamma,const_phi,v_i,dt,eps,elliptic=False,penalty=False)

        if converged == False:
            # exit loop and remove nans
            N = np.nan_to_num(N)
            N[N<1e-7] = 1e-7
            print('\n failed to converge...')
            break

        # save the solution as numpy arrays
        z_i,N_i = interp(sol.sub(0),domain)
        z_i,heave_i = interp(sol.sub(1),domain)  
        z_i,visc_i = interp(N_visc(sol.sub(1),Phi(sol.sub(0),const_phi=const_phi),gamma),domain)
 
        N[i,:] = N_i
        heave[i,:] = heave_i 
        visc[i,:] = visc_i
        z[i,:] = z_i


        # compute minimum effective stress
        
        N_min = (N_i+visc_i).min()

        if N_min <= N_c: 
        #  initiate new lens position where effective stress reaches zero
           
            # set flag to 1 for lens initiation times 
            new_lens[i] = 1 

            # create new domain [z_b, z_n] where z_n corresponds to minimum
            # where the effective stress N reaches threshold          
            l = np.argmin(N_i+visc_i)
            z_n = z_i[l]
            z_b = z_i.min()
            domain = create_interval(MPI.COMM_WORLD,nz,[z_b,z_n])

            # initialize solution on new domain        
            sol_n = initialize(domain,N_f,gamma,const_phi,eps_min=eps)
    
            z_i,N_i = interp(sol_n.sub(0),domain)
            z_i,heave_i = interp(sol_n.sub(1),domain)  
            z_i,visc_i = interp(N_visc(sol_n.sub(1),Phi(sol_n.sub(0),const_phi=const_phi),gamma),domain)
    
            N[i,:] = N_i
            heave[i,:] = heave_i 
            visc[i,:] = visc_i
            z[i,:] = z_i
         
        else:
            # evolve domain geometry according to ice velocity
            # and heave rate
            v_s = v_i - heave_i[-1]

            Z = domain.geometry.x
            z_l = Z[:,0].max()
            z_b = Z[:,0].min()
            Z[:,0] += dt*v_s*(Z[:,0]-z_b)/(z_l-z_b)

            P1 = FiniteElement('P',domain.ufl_cell(),1)     
            element = P1*P1
            V = FunctionSpace(domain,element)               
            sol_n = Function(V)
            sol_n.sub(0).interpolate(sol.sub(0))
            sol_n.sub(1).interpolate(sol.sub(1))
    
    # remove nans at temporal jumps
    heave = np.nan_to_num(heave)
    N = np.nan_to_num(N)
    visc = np.nan_to_num(visc)        

    return N,heave,visc, z, new_lens, converged


def initialize(domain,N_f,gamma,const_phi,eps_min=1e-10):
    # initialize effective stress solution given a domain
    # and the boundary condition at the base
    # this assumes that solution is rigid (dN/dt + v_i dN/dz = 0)
    # so that N solves an elliptic pde 

    P1 = FiniteElement('P',domain.ufl_cell(),1)     
    element = P1*P1
    V = FunctionSpace(domain,element)    
    z_l = domain.geometry.x[:,0].max()
    z_b = domain.geometry.x[:,0].min()
    sol_n = Function(V)
    sol_n.sub(0).interpolate(lambda x: 0.001*(x[0]-z_b)/(z_l-z_b) + N_f)
    sol_n, converged = solve_pde(domain,sol_n,N_f,gamma,const_phi,0,1,eps=1e-2,penalty=True,elliptic=True)
   
    # decrease regularization parameter 
    eps_ = np.flipud(np.logspace(np.log10(eps_min),-2,20))
    for j in range(eps_.size):
        sol_n, converged = solve_pde(domain,sol_n,N_f,gamma,const_phi,0,1,eps=eps_[j],penalty=False,elliptic=True)

    return sol_n


def heave_plastic(N,domain,const_phi):
    phi = Phi(N,const_phi=const_phi)
    x = SpatialCoordinate(domain)
    T,S,k = get_fields(x[0])
    v_h =  (Dx(N,0)+ (1-phi) + (1+T)*Dx(phi*S,0))*k/((1-phi*S)**2)
    V = FunctionSpace(domain, ("CG", 1))
    v = Function(V)
    v.interpolate(Expression(v_h, V.element.interpolation_points()))
    return v


  