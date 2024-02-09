# this file contains some helper functions for the compaction problem solver

import numpy as np
from params import Tz_sc, a, alpha, b, d0, e0
from ufl import ln

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

def Phi(N,log=ln):
     # porosity as a function of scaled effective stress
     e = e0 - d0*log(N/alpha)/np.log(10)
     return e/(e+1)
     # # constant porosity:
     # return 0.325 + 1e-20*N 

def dPhi(N,log=ln):
     # (absolute value of) derivative of Phi w.r.t. N
     return d0*np.log(10)/(N*(-d0*log(N/alpha) + np.log(10)*(e0 + 1))**2)
     # # constant porosity:
     # return 1e-20*N

def get_fields(z):
     # get saturation, temperature, and permeability fields
     # given the z coordinate
     T = temp(z)
     S = sat(T)
     k = perm(S)
     return T,S,k