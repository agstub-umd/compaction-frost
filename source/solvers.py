# This file contains the functions needed for solving the compaction problem.
import numpy as np
from dolfinx.fem import (Expression, Function, FunctionSpace, dirichletbc,
                         locate_dofs_topological)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import locate_entities_boundary
from dolfinx.nls.petsc import NewtonSolver
from misc import C, D, K, interp, max, move_mesh, sat, temp, w_ice
from mpi4py import MPI
from petsc4py import PETSc
from params import Ts, dt, nt, nz, phi_min, theta
from ufl import (Dx, FiniteElement, SpatialCoordinate, TestFunction,
                 TestFunctions, ds, dx, split)


def weak_form(w,w_t,w_n,phi,phi_t,phi_n,Gamma,domain):
    # Weak form of the residual for the compaction problem
    w_theta = theta*w + (1-theta)*w_n
    phi_theta = theta*phi + (1-theta)*phi_n

    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T)
    wi = w_ice(w,phi,S,Gamma)

    # weak form of momentum balance:
    F_w =  D(phi)*Dx(w,0)*Dx(w_t,0)*dx  + (C(phi,S)-Gamma)*w_t*dx + (1/K(S))*((1-phi*S)*w + phi*S*wi)*w_t*dx
    F_w += phi*S*(1+T)*Dx(w_t,0)*dx - phi*S*(1+Ts)*w_t*ds 

    # weak form of porosity evolution:
    F_phi = (phi-phi_n)*phi_t*dx + dt*w_theta*Dx(phi_theta,0)*phi_t*dx - dt*(1-phi_theta)*Dx(w_theta,0)*phi_t*dx 

    # add constraint phi>phi_min:
    F_phi += (phi-max(phi_min,phi))*phi_t*dx
    return F_w + F_phi

def weak_form_vel(w,w_t,phi,Gamma,domain):
    # Non-coupled problem - solve momentum balance for a fixed porosity: 
    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T)
    wi = w_ice(w,phi,S,Gamma) 
    F_w =  D(phi)*Dx(w,0)*Dx(w_t,0)*dx  + (C(phi,S)-Gamma)*w_t*dx + (1/K(S))*((1-phi*S)*w + phi*S*wi)*w_t*dx
    F_w += phi*S*(1+T)*Dx(w_t,0)*dx - phi*S*(1+Ts)*w_t*ds 
    return F_w 


def solve_pde(domain,sol_n,Gamma):
        # solves the compaction PDE problem at a given time step

        # Define function space
        P1 = FiniteElement('P',domain.ufl_cell(),1)     
        element = P1*P1
        V = FunctionSpace(domain,element)       

        sol = Function(V)
        (w,phi) = split(sol)
        (w_n,phi_n) = split(sol_n)
        (w_t,phi_t) = TestFunctions(V)

        # # Define weak form:
        F = weak_form(w,w_t,w_n,phi,phi_t,phi_n,Gamma,domain)

        # Define Dirichlet boundary condition for base if you want
        L = domain.geometry.x.min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],L))
        dofs_b = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(0), dofs_b,V.sub(0))      # w = 0 at base  

        bcs = [] #[bc_b]    

        # set initial guess for Newton solver to be the solution 
        # from the previous time step:
        sol.sub(0).interpolate(sol_n.sub(0))
        sol.sub(1).interpolate(sol_n.sub(1))

        # Solve for sol = (w,phi):
        problem = NonlinearProblem(F, sol, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
        solver.solve(sol)

        # bound porosity below by ~phi_min: 
        V0 = FunctionSpace(domain, ("CG", 1))
        Phi = Function(V0)
        Phi.interpolate(sol.sub(1))
        Phi.x.array[Phi.x.array<phi_min] = 1.01*phi_min
        Phi.x.array[Phi.x.array>0.999] = 0.999
        sol.sub(1).interpolate(Phi)

        x = SpatialCoordinate(domain)
        T = temp(x[0])
        S = sat(T)
        wi_0 = w_ice(sol.sub(0),sol.sub(1),S,Gamma)

        V0 = FunctionSpace(domain, ("CG", 1)) 
        wi = Function(V0)
        wi.interpolate(Expression(wi_0, V0.element.interpolation_points()))

        return sol,wi

def full_solve(domain,initial,m,Gamma):
    # solve the compaction problem given:
    # domain: the computational domain
    # initial: initial conditions 
 
    # *see example.ipynb for an example of how to set these
    #
    # The solution sol = (u,phii) returns:
    # ws: vertical velocity (sediment)
    # phi: porosity
    # We also save:
    # z: domain coordinates (change over time due to compaction)
    # wi: ice velocity

    ws_arr = np.zeros((nt,nz+1))
    wi_arr = np.zeros((nt,nz+1))
    phi_arr = np.zeros((nt,nz+1))
    z_arr = np.zeros((nt,nz+1))

    sol_n = initial
    phi_i = phi_arr

    # time-stepping loop
    for i in range(nt):
        print('time step '+str(i+1)+' out of '+str(nt)+' \r',end='')

        # solve the compaction problem for sol = (w,phi)
        sol,wi = solve_pde(domain,sol_n,Gamma)
        
        # # displace the mesh according to the velocity solution
        domain = move_mesh(domain,wi,m)

        # set the solution at the previous time step
        sol_n.sub(0).interpolate(sol.sub(0))
        sol_n.sub(1).interpolate(sol.sub(1))
     
        # save the solution as numpy arrays
        z_i,ws_i = interp(sol_n.sub(0),domain)
        z_i,wi_i = interp(wi,domain)
        z_i,phi_i = interp(sol_n.sub(1),domain)

        ws_arr[i,:] = ws_i
        wi_arr[i,:] = wi_i
        phi_arr[i,:] = phi_i
        z_arr[i,:] = z_i
       

    return ws_arr,phi_arr, wi_arr, z_arr


def vel_solve(domain,phi,Gamma):
        # Solve the momentum balance for a fixed porosity

        # Define function space
        V = FunctionSpace(domain, ("CG", 1))
   
        w = Function(V)
        w_t = TestFunction(V)
 
   
        # Define weak form
        F = weak_form_vel(w,w_t,phi,Gamma,domain)

        # Define Dirichlet boundary condition for base if you want
        L = domain.geometry.x.min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],L))
        dofs_b = locate_dofs_topological(V, domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(0), dofs_b,V)    # w = 0 at base  

        bcs = [] #[bc_b] 

        # initial guess for Newton solver
        w.interpolate(lambda x: (x[0]-L)/L )
        # Solve for w
        problem = NonlinearProblem(F, w, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
      
        solver.solve(w)
      
        return w