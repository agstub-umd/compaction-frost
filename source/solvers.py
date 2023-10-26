# This file contains the functions needed for solving the compaction problem.
import numpy as np
from constitutive import Nc, perm, q, sat, temp, qr
from dolfinx.fem import (assemble_scalar, Expression, form,Function,FunctionSpace, dirichletbc,
                         locate_dofs_topological)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.mesh import locate_entities_boundary,create_interval
from dolfinx.nls.petsc import NewtonSolver
from fem_spaces import mixed_space, vel_space
from mpi4py import MPI
from params import G, nu, nz, theta, phi_min, phi_max
from petsc4py import PETSc
from post_process import interp
from scipy.interpolate import interp1d
from ufl import (Dx, SpatialCoordinate, TestFunction, TestFunctions,
                 conditional, ds, dx, ge, le, split)


def Nr(domain,phi,Nf):
    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T) 
    K = perm(S)
    V_heave = heave_rate(domain,phi,Nf)

    V = vel_space(domain)
    N = Function(V)
    N_t = TestFunction(V)
    weak_form = Dx(N,0)*N_t*dx  + G*(1-phi)*(nu-1)*N_t*dx - (1-phi*S)*(qr(V_heave,phi,S)/K)*N_t*dx + (1+T)*Dx(phi*S,0)*N_t*dx  

    zb = domain.geometry.x[:,0].min()
    facets = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
    dofs = locate_dofs_topological(V, domain.topology.dim-1, facets)
    bc = dirichletbc(PETSc.ScalarType(Nf), dofs,V)    
    bcs = [bc]   

    problem = NonlinearProblem(weak_form, N, bcs=bcs)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.solve(N) 

    return N, V_heave


def heave_rate(domain,phi,Nf):
    # rigid V_heave = v_i-v_s
    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T) 
    K = perm(S)
    int_0 = G*(1-phi)*(nu-1)*dx + (1+T)*Dx(phi*S,0)*dx 
    int_0 = assemble_scalar(form(int_0))
    int_1 = ((1-phi*S)**2/K)*dx
    int_1 = assemble_scalar(form(int_1))
    return  (1-Nf + int_0)/int_1

def weak_form(w,w_t,w_n,phi,phi_t,phi_n,V_s,domain,dt):
    # Weak form of the residual for the compaction problem
    w_theta = theta*w + (1-theta)*w_n
    phi_theta = theta*phi + (1-theta)*phi_n

    # weak form of momentum balance:
    F_w =  weak_form_vel(w,w_t,phi,domain)

    # weak form of porosity evolution:
    F_phi = dt*(w_theta+V_s)*Dx(phi_theta,0)*phi_t*dx - dt*(1-phi_theta)*Dx(w_theta,0)*phi_t*dx 
    F_phi += (phi-phi_n)*phi_t*dx 

    return F_w + F_phi

def weak_form_vel(w,w_t,phi,domain):
    # weak form for unfrozen (sediment + water) momentum balance
    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T) 
    K = perm(S)
    
    # body integral terms:
    F_w =  -Nc(phi,w)*Dx(w_t,0)*dx  - (1-phi*S)*(q(w,phi,S)/K)*w_t*dx

    # boundary integral terms:
    # F_w -= dNf*w_t*ds # effective stress perturbation at base of fringe

    return F_w 



def coupled_solve(domain,sol_n,q_in,V_s,dt):
        # solves the compaction PDE problem at a given time step

        # Define function space (in fem_spaces.py)
        V = mixed_space(domain)       

        sol = Function(V)
        (w,phi) = split(sol)
        (w_n,phi_n) = split(sol_n)
        (w_t,phi_t) = TestFunctions(V)

        # sediment velocity at base of ice lens relative to rigid solution (zero)
        zl = domain.geometry.x[:,0].max()
        facets_l = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zl))
        dofs_l = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_l)
        bc_l = dirichletbc(PETSc.ScalarType(0), dofs_l,V.sub(0))    

        zb = domain.geometry.x[:,0].min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs_b = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(-q_in), dofs_b,V.sub(0))  

        bcs = [bc_l,bc_b]

        # # Define weak form:
        F = weak_form(w,w_t,w_n,phi,phi_t,phi_n,V_s,domain,dt)

        # Solve for sol = (w,phi):
        problem = NonlinearProblem(F, sol, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)

        # set initial guess for Newton solver to be the solution 
        # from the previous time step:
        sol.sub(0).interpolate(sol_n.sub(0))
        sol.sub(1).interpolate(sol_n.sub(1))
   
        solver.solve(sol)

        V0 = vel_space(domain)
        N_c = Function(V0)

        N_expr = Nc(sol.sub(1),sol.sub(0))

        N_c.interpolate(Expression(N_expr, V0.element.interpolation_points()))
        
        # bound porosity below by phi_min and above by phi_max: 
        V0 = FunctionSpace(domain, ("CG", 1))
        Phi = Function(V0)
        Phi.interpolate(sol.sub(1))
        Phi.x.array[Phi.x.array<phi_min] = phi_min
        Phi.x.array[Phi.x.array>phi_max] = phi_max
        sol.sub(1).interpolate(Phi)

        return sol, N_c

def vel_solve(domain,phi,q_in):
        # solves the momentum balance for a fixed porosity field

        # Define function space
        V = vel_space(domain) 

        w = Function(V)
        w_t = TestFunction(V)
 
        # Define weak form
        F = weak_form_vel(w,w_t,phi,domain) 

        # sediment velocity at base of ice lens relative to rigid solution (zero)
        zl = domain.geometry.x[:,0].max()
        facets_l = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zl))
        dofs_l = locate_dofs_topological(V, domain.topology.dim-1, facets_l)
        bc_l = dirichletbc(PETSc.ScalarType(0), dofs_l,V)  

        zb = domain.geometry.x[:,0].min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs_b = locate_dofs_topological(V, domain.topology.dim-1, facets_b)
        bc_b = dirichletbc(PETSc.ScalarType(-q_in), dofs_b,V)      

        bcs = [bc_l,bc_b]

        # Solve for w
        problem = NonlinearProblem(F, w, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)

        solver.solve(w)

        N_c = Nc(phi,w)
      
        return w,N_c


def full_solve(domain,initial,Nf,V_pull,V_heave,q_f,timesteps):
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

    dt = timesteps[1]-timesteps[0]
    nt = np.size(timesteps)

    ws_arr = np.zeros((nt,nz+1))
    phi_arr = np.zeros((nt,nz+1))
    N_arr = np.zeros((nt,nz+1))
    Nr_arr = np.zeros((nt,nz+1))
    z_arr = np.zeros((nt,nz+1))
    Vh_arr = np.zeros(nt)

    sol_n = initial
    phi_i = phi_arr

    

    # time-stepping loop
    for i in range(nt):
        print('time step '+str(i+1)+' out of '+str(nt)+' \r',end='')

        # solve the compaction problem for sol = (w,phi)
        q_in = q_f - V_heave
        V_s = V_pull-V_heave
        sol,N_c = coupled_solve(domain,sol_n,q_in,V_s,dt)
        
        # set the solution at the previous time step
        sol_n.sub(0).interpolate(sol.sub(0))
        sol_n.sub(1).interpolate(sol.sub(1))
     
        N_r,V_heave = Nr(domain,sol.sub(1),Nf)
        N = N_r + N_c

        # save the solution as numpy arrays
        z_i,ws_i = interp(sol_n.sub(0),domain)
        z_i,phi_i = interp(sol_n.sub(1),domain)
        z_i,N_i = interp(N,domain)
        z_i,Nr_i = interp(N_r,domain)

        N_arr[i,:] = N_i
        Nr_arr[i,:] = Nr_i
        ws_arr[i,:] = ws_i
        phi_arr[i,:] = phi_i
        z_arr[i,:] = z_i
        Vh_arr[i] = V_heave 
      

        V_s = V_pull - V_heave + ws_i[-1]
       

        N_min = N_i.min()

        if N_min <= 0:
        # initiate new lens position where effective stress vanishes
            l = np.argmin(N_i)
            z_n = z_i[l]
            zb = z_i.min()
            z_new = np.linspace(zb,z_n,nz+1)
            phi_n = interp1d(z_i[0:l],phi_i[0:l],fill_value='extrapolate')
            ws_n = interp1d(z_i[0:l],ws_i[0:l],fill_value='extrapolate')
            domain = create_interval(MPI.COMM_WORLD,nz,[zb,z_n])
            V0 = FunctionSpace(domain, ("CG", 1))
            V = mixed_space(domain)
            phi_ = Function(V0)
            ws_ = Function(V0)
            phi_.x.array[:] = phi_n(z_new)
            ws_.x.array[:] = ws_n(z_new)
            sol_n = Function(V)
            sol_n.sub(0).interpolate(ws_)
            sol_n.sub(1).interpolate(phi_)
        else:
            z = domain.geometry.x
            zl = z[:,0].max()
            zb = z[:,0].min()
            z[:,0] += dt*V_s*(z[:,0]-zb)/(zl-zb)

    return ws_arr,phi_arr,N_arr,Nr_arr,Vh_arr, z_arr
