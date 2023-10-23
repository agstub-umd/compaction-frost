# This file contains the functions needed for solving the compaction problem.
import numpy as np
from constitutive import N, Q, Sigma, perm, q, sat, temp, Pi, Lamda
from dolfinx.fem import (Expression, Constant, Function, FunctionSpace, dirichletbc,
                         locate_dofs_topological)
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx import log
from dolfinx.mesh import create_interval,locate_entities_boundary
from dolfinx.nls.petsc import NewtonSolver
from mpi4py import MPI
from params import CFL, G, nu, nz, phi_max, phi_min, theta
from petsc4py import PETSc
from post_process import interp
from ufl import (Dx, SpatialCoordinate,
                 TestFunctions, conditional, ds, dx, ge, le, split)
from fem_spaces import mixed_space,vel_space
from scipy.interpolate import interp1d

def weak_form(w,w_t,w_n,phi,phi_t,phi_n,wi,wi_t,Nf,domain,dt):
    # Weak form of the residual for the compaction problem
    w_theta = theta*w + (1-theta)*w_n
    phi_theta = theta*phi + (1-theta)*phi_n

    # weak form of momentum balance:
    F_w =  weak_form_sed(w,w_t,phi,wi,Nf,domain) + weak_form_ice(w,phi,wi,wi_t,domain)

    # weak form of porosity evolution:
    F_phi = dt*w_theta*Dx(phi_theta,0)*phi_t*dx - dt*(1-phi_theta)*Dx(w_theta,0)*phi_t*dx 
    F_phi += (phi-phi_n)*phi_t*dx 

    return F_w + F_phi

def weak_form_sed(w,w_t,phi,wi,Nf,domain):
    # weak form for unfrozen (sediment + water) momentum balance
    x = SpatialCoordinate(domain)
    T = temp(x[0])
    S = sat(T) 
    K = perm(S)

    zl = domain.geometry.x[:,0].max()
    zb = domain.geometry.x[:,0].min()
    zm = 0.5*(zb+zl)
    
    # body integral terms:
    F_w =  -N(phi,w)*Dx(w_t,0)*dx  + G*(1-phi)*(nu-1)*w_t*dx - (1-phi*S)*(q(w,wi,phi,S)/K)*w_t*dx
    F_w += (1+T)*Dx(phi*S,0)*w_t*dx  

    # boundary integral terms:
    # Nf0 = Pi(phi) - Lamda(phi)*w*Dx(phi,0)
    F_w -= Nf*conditional(le(x[0],zm),1,0)*w_t*ds # effective stress at base of fringe

    Nl = (1-phi*S)*(1+T) + Sigma(wi,phi,S)        # effective stress below ice lens
    F_w += Nl*w_t*conditional(ge(x[0],zm),1,0)*ds 
    return F_w 

def weak_form_ice(w,phi,wi,wi_t,domain):
     # weak form for ice momentum balance
     x = SpatialCoordinate(domain)
     T = temp(x[0])
     S = sat(T)
     K = perm(S)

     zl = domain.geometry.x[:,0].max()
     zb = domain.geometry.x[:,0].min()
     zm = 0.5*(zb+zl)

     m = Q(phi,S)/(1-phi*S)
  
     # body integral terms:
     F_i = -Sigma(wi,phi,S)*Dx(wi_t,0)*dx - phi*S*Dx(T,0)*wi_t*dx  + phi*S*(q(w,wi,phi,S)/K)*wi_t*dx

     # boundary integral terms: 
     pen = Constant(domain, PETSc.ScalarType(1e6))
     F_i += pen*(wi-(w+m))*conditional(ge(x[0],zm),1,0)*wi_t*ds  
     return F_i 

def pde_solve(Nf,domain,sol_n,dt):
        # solves the compaction PDE problem at a given time step

        # Define function space (in fem_spaces.py)
        V = mixed_space(domain)       

        sol = Function(V)
        (w,wi,phi) = split(sol)
        (w_n,wi_n,phi_n) = split(sol_n)
        (w_t,wi_t,phi_t) = TestFunctions(V)

   
        w_b =  0.0 

        # set Dirichlet bc's
        zb = domain.geometry.x[:,0].min()
        facets = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets)
        bc_s = dirichletbc(PETSc.ScalarType(w_b), dofs,V.sub(0))    # ws at base of fringe 

        dofs_i = locate_dofs_topological(V.sub(1), domain.topology.dim-1, facets)
        bc_i = dirichletbc(PETSc.ScalarType(w_b), dofs_i,V.sub(1))    # wi at base of fringe 

        # bcs = [bc_s,bc_i]
        bcs = [bc_i]

        # # Define weak form:
        F = weak_form(w,w_t,w_n,phi,phi_t,phi_n,wi,wi_t,Nf,domain,dt)

        # Solve for sol = (w,phi):
        problem = NonlinearProblem(F, sol, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
        solver.max_it = 100
        solver.error_on_nonconvergence = False
   
        ksp = solver.krylov_solver
        opts = PETSc.Options()
        option_prefix = ksp.getOptionsPrefix()
        opts[f"{option_prefix}pc_type"] = "ksp"
        ksp.setFromOptions()

        # set initial guess for Newton solver to be the solution 
        # from the previous time step:
        sol.sub(0).interpolate(sol_n.sub(0))
        sol.sub(1).interpolate(sol_n.sub(1))
        sol.sub(2).interpolate(sol_n.sub(2))

        n, converged = solver.solve(sol)
      
        # bound porosity below by phi_min and above by phi_max: 
        V0 = FunctionSpace(domain, ("CG", 1))
        Phi = Function(V0)
        Phi.interpolate(sol.sub(2))
        Phi.x.array[Phi.x.array<phi_min] = phi_min
        Phi.x.array[Phi.x.array>phi_max] = phi_max
        sol.sub(2).interpolate(Phi)

        return sol,converged

def initial_solve(domain,phi,Nf):
        # solves the momentum balance for a fixed porosity field

        # Define function space
        V = vel_space(domain) 

        sol = Function(V)
        (w,wi) = split(sol)
        (w_t,wi_t) = TestFunctions(V)
 
        # Define weak form
        F = weak_form_sed(w,w_t,phi,wi,Nf,domain) + weak_form_ice(w,phi,wi,wi_t,domain)


        # set sediment velocity at base related to water flux bc 
        w_b = 0.0 
    
        # Define Dirichlet boundary conditions at base of fringe (wi=0)
        zb = domain.geometry.x[:,0].min()
        facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],zb))
        dofs_i = locate_dofs_topological(V.sub(1), domain.topology.dim-1, facets_b)
        bc_i = dirichletbc(PETSc.ScalarType(w_b), dofs_i,V.sub(1)) #wi at base of fringe
     
        dofs_s = locate_dofs_topological(V.sub(0), domain.topology.dim-1, facets_b)
        bc_s = dirichletbc(PETSc.ScalarType(w_b), dofs_s,V.sub(0))    # ws at base of fringe

        # bcs = [bc_i,bc_s]
        bcs = [bc_i]

        # Solve for w
        problem = NonlinearProblem(F, sol, bcs=bcs)
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
        solver.max_it = 100

        solver.solve(sol)

        N_ = N(phi,sol.sub(0))
      
        return sol,N_

def full_solve(domain,initial,Nf,t_0,t_max,dt0):
    # solve the compaction problem given:
    # domain: the computational domain
    # initial: initial conditions 
 
    # *see example.ipynb for an example of how to set these
    #
    # The solution sol = (w,wi,phi) returns:
    # w: vertical velocity (sediment)
    # wi: vertical velocity (ice)
    # phi: porosity
    # We also save:
    # z: domain coordinates (change over time due to compaction)
    log.set_log_level(log.LogLevel.ERROR)

    nt0 = int(t_max/dt0)
    nt = 10*nt0
    dt = dt0

    # arrays for saving results
    ws_arr = np.zeros((nt,nz+1))    # sediment velocity
    wi_arr = np.zeros((nt,nz+1))    # ice velocity 
    phi_arr = np.zeros((nt,nz+1))   # porosity
    z_arr = np.zeros((nt,nz+1))     # vertical coordinate
    N_arr = np.zeros((nt,nz+1))     # effective stress (unfrozen material)
    timesteps = np.zeros(nt)        # timesteps 

    sol_n = initial
    i = 0       # time step
    j = 0       # total number of solve attempts
    t = t_0
    converged = True  # convergence flag
    C = 0.5          # CFL constant

    C_fine0 = 0.1     # CFL constant for fine timestepping
    C_fine = C_fine0

    dt_prev = dt0

    N_min = Nf

    # time-stepping loop
    c_count = 16
    print('Time stepping...')
    while t<t_max and (i<nt and j<10*nt):

        # set the timestep size according to CFL condition
        z = domain.geometry.x[:,0]
        dz = np.abs(z[1]-z[0])
        vs_max = np.abs(ws_arr[i-1,:]).max()

        if c_count>15:
            C_fine = C_fine0
            dt = np.min([C*dz/(vs_max+1e-10),dt0])
        else:
            dt = np.min([C_fine*dz/(vs_max+1e-10),dt0]) 
                    
        if dt > 10*dt_prev and i>2:
            dt = 2*dt_prev   

        if i>0:
            print('iteration '+str(i)+': time = '+'{:.3f}'.format(t)+' out of '+'{:.2f}'.format(t_max)+',  dt = '+'{:.1f}'.format(100*dt/dt0)+'% of dt0, N_min = '+'{:.3f}'.format(N_min)+'  \r',end='')

        # solve the compaction problem for sol = (w_s,w_i,phi)
        sol,converged = pde_solve(Nf,domain,sol_n,dt)

        if converged == False:
             c_count = 0
             sol.sub(0).interpolate(sol_n.sub(0))
             sol.sub(1).interpolate(sol_n.sub(1))
             sol.sub(2).interpolate(sol_n.sub(2))
             if C_fine > 1e-9:
                C_fine *= 0.5
             else:
                print('\n failed to converge... :(')
                break
                     
        elif converged == True:
            c_count += 1    
            dt_prev = dt

            # set the solution at the previous time step
            sol_n.sub(0).interpolate(sol.sub(0))
            sol_n.sub(1).interpolate(sol.sub(1))
            sol_n.sub(2).interpolate(sol.sub(2))
        
            # save the solution as numpy arrays
            # x = SpatialCoordinate(domain)
            # T = temp(x[0])
            # S = sat(T,N(sol.sub(2),sol.sub(0)))
            z_i,ws_i = interp(sol_n.sub(0),domain)
            z_i,wi_i = interp(sol_n.sub(1),domain)
            z_i,phi_i = interp(sol_n.sub(2),domain)
            z_i,N_i = interp(N(sol.sub(2),sol.sub(0)),domain)

            ws_arr[i,:] = ws_i
            wi_arr[i,:] = wi_i
            N_arr[i,:] = N_i
            phi_arr[i,:] = phi_i
            z_arr[i,:] = z_i
            timesteps[i] = t

            # minimum effective stress
            N_min = N_i.min()
            
            if N_min <= 0:
                # initiate new lens position where effective stress vanishes
                l = np.argmin(N_i)
                z_n = z_i[l]
                zb = z_i.min()
                z_new = np.linspace(zb,z_n,nz+1)
                phi_n = interp1d(z_i[0:l],phi_i[0:l],fill_value='extrapolate')
                ws_n = interp1d(z_i[0:l],ws_i[0:l],fill_value='extrapolate')
                wi_n = interp1d(z_i[0:l],wi_i[0:l],fill_value='extrapolate')
                domain = create_interval(MPI.COMM_WORLD,nz,[zb,z_n])
                V0 = FunctionSpace(domain, ("CG", 1))
                V = mixed_space(domain)
                phi_ = Function(V0)
                ws_ = Function(V0)
                wi_ = Function(V0)
                phi_.x.array[:] = phi_n(z_new)
                ws_.x.array[:] = ws_n(z_new)
                wi_.x.array[:] = wi_n(z_new)
                sol_n = Function(V)
                sol_n.sub(0).interpolate(ws_)
                sol_n.sub(1).interpolate(wi_)
                sol_n.sub(2).interpolate(phi_)

            else:    
                # displace the mesh according to the velocity solution
                z = domain.geometry.x
                zl = z[:,0].max()
                zb = z[:,0].min()
                z[:,0] += dt*(ws_i[-1]+0*1)*(z[:,0]-zb)/(zl-zb)
       
            t += dt
            i += 1

        j += 1 # total number of solve attempts 
    ws_arr = ws_arr[0:i-1,:]
    wi_arr= wi_arr[0:i-1,:]
    N_arr = N_arr[0:i-1,:]
    phi_arr = phi_arr[0:i-1,:]
    z_arr = z_arr[0:i-1,:]
    timesteps = timesteps[0:i-1]    

    return ws_arr,phi_arr, wi_arr, N_arr, z_arr,sol,domain,timesteps