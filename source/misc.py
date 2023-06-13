# this file contains some helper functions for the compaction problem solver

import numpy as np
from dolfinx.fem import (Constant, Expression, Function, FunctionSpace,
                         dirichletbc, locate_dofs_topological)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.mesh import locate_entities_boundary
from params import (T0, G, Lh, Tf, Tm, a, b, delta, dt, g, nz, phi0, rho_i,
                    rho_s, rho_w)
from petsc4py import PETSc
from scipy.interpolate import griddata
from ufl import Dx, TestFunction, TrialFunction, dx


def K(S):
      # scaled permeability as a function of ice saturation
      return (1-S)**a

def max(f1,f2):
     # max function
     return 0.5*(f1+f2 + ((f1-f2)**2)**0.5)

def sign(f):
     # regularized sign function
     return f/(abs(f)+1e-20)

def sat(T):
     # ice saturation as a function of temperature
     f = 1 - abs(1/(1+T))**b
     return max(0*T,f)*(1+sign(T))/2.

def temp(z):
     # temperature as a gunction of z
     return (Tf-T0 + G*delta*z)/(Tm-Tf)

def C(phi,S):
     # coefficient 
     rho_b = (1-phi)*rho_s + phi*S*rho_i + phi*(1-S)*rho_w
     return (rho_b-rho_w)*g*delta/(rho_i*Lh*(1-Tf/Tm))

def D(phi):
     # coefficient on dw/dz in weak form
     return (1-phi)*(4./3. + 1/phi)/(4./3. + 1/phi0)

def w_ice(w,phi,S,Gamma):
     f = rho_w*Gamma + (1-phi)*(rho_s-rho_w)*w 
     f = f/(rho_i+phi*S*(rho_w-rho_i)) 
     return f 

def interp(f,domain):
    # returns a numpy array of a (dolfinx) function f that has been
    # evaluated at the mesh nodes
    V = FunctionSpace(domain, ("CG", 1))
    u = Function(V)
    u.interpolate(Expression(f, V.element.interpolation_points()))

    z = domain.geometry.x[:,0]
    vals = u.x.array

    Z = np.linspace(z.min(),z.max(),nz+1)

    points = (z)
    values = vals
    points_i = Z

    F = griddata(points, values, points_i, method='linear')    

    return points_i,F


def move_mesh(domain,wi,m):
    # this function computes the surface displacements and moves the mesh
    # by solving Laplace's equation for a smooth displacement function
    # defined for all mesh vertices. (*Solving Laplace is not really necessary 
    # for this 1d problem but would be helpful for higher dimensions)

    V = FunctionSpace(domain, ("CG", 1))
    L = domain.geometry.x.min()

    w_top = wi.x.array[-1] + m

    facets_t = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],0))
    facets_b = locate_entities_boundary(domain, domain.topology.dim-1, lambda x: np.isclose(x[0],L))
        
    dofs_t = locate_dofs_topological(V, domain.topology.dim-1, facets_t)
    dofs_b = locate_dofs_topological(V, domain.topology.dim-1, facets_b)

    bc_top = dirichletbc(PETSc.ScalarType(w_top), dofs_t,V)  # displacement = w at top
    bc_base = dirichletbc(PETSc.ScalarType(0), dofs_b,V)     # w = 0 at base  

    bcs = [bc_top,bc_base]

    # # solve Laplace's equation for a smooth displacement field on all vertices,
    # # given the boundary displacement bc's
    disp = TrialFunction(V)
    v = TestFunction(V)
    a = Dx(disp,0)*Dx(v,0)*dx 
    f = Constant(domain, PETSc.ScalarType(0.0))
    L = f*v*dx

    problem = LinearProblem(a,L, bcs=bcs)
    sol = problem.solve()

    disp_vv = sol.x.array

    z = domain.geometry.x

    z[:,0] += dt*disp_vv

    return domain

