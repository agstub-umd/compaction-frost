# this file contains some helper functions for the compaction problem solver

import numpy as np
from params import Tz_sc, a, alpha, b, d0, e0
from ufl import ln, Dx
from scipy.optimize import root_scalar

def Max(f1,f2):
     # max function
     return 0.5*(f1+f2 + ((f1-f2)**2)**0.5)

def sign(f):
     # regularized sign function
     return f/(abs(f)+1e-20)

def sat(T):
     # ice saturation as a function of temperature
     f = 1 - abs(1/(1+T))**b 
     S0 = Max(0*T,f)*(1+sign(T))/2.
     return S0

def sat_T(T):
     # derivative of saturation with respect to T
     return 0.5*sign(T + 1)/np.abs(T + 1)**1.5       

def perm(S):
      # scaled permeability as a function of ice saturation
      return (1-S)**a

def temp(z):
     # undercooling temperature as a function of z
     return Tz_sc*z

def D(phi,S,eps=1e-10):
     k = perm(S)
     d = ((1-phi*S)**2)/ k
     return (1-phi)/d +  eps

def Phi(N,log=ln,const_phi=False):
     # porosity as a function of scaled effective stress
     if const_phi == False:
          # empirical consolidation law
          e = e0 - d0*log(N/alpha)/np.log(10)
          P = e/(e+1)
     else:
          # constant porosity
          P = 0.325 + 1e-20*N 
     return P
   
def dPhi(N,log=ln,const_phi=False):
     # (absolute value of) derivative of Phi w.r.t. N
     if const_phi == False:
          # empirical consolidation law
          dP = d0*np.log(10)/(N*(-d0*log(N/alpha) + np.log(10)*(e0 + 1))**2)
     else:
          # constant porosity
          dP = 1e-20*N     
     return dP


def get_fields(z):
     # get saturation, temperature, and permeability fields
     # given the z coordinate
     T = temp(z)
     S = sat(T)
     k = perm(S)
     return T,S,k

def N_visc(v,phi,gamma):
     # viscous component of effective stress
     # due to ice creep in/out of pores
     if gamma>0:
          N_v = (1-phi)*gamma*Dx(v,0)
     else:
          N_v = 1e-20*v     
     return N_v

def Phi_lens(z):
    res = lambda phi,z: phi - Phi((1-phi*sat(temp(z)))*(1+temp(z)),log=np.log)
    return root_scalar(lambda phi: res(phi,z),x0 = 0.25,bracket = [0,1]).root   