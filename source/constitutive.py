# this file contains some helper functions for the compaction problem solver

import numpy as np
from params import (Ki, Ks, Kw, Tf, Tm, Ts, Tz0, Tz_sc, a, alpha, b, beta,
                    gamma, z_sc)
from ufl import Dx


def Max(f1,f2):
     # max function
     return 0.5*(f1+f2 + ((f1-f2)**2)**0.5)

def perm(S):
      # scaled permeability as a function of ice saturation
      return (1-S)**a

def sign(f):
     # regularized sign function
     return f/(abs(f)+1e-20)

def step(f):
     # smoothed step function
     return 1/(1+np.exp(1)**(-200*f))

def sat(T):
     # ice saturation as a function of temperature
     f = 1 - abs(1/(1+T))**b
     S0 = Max(0*T,f)*(1+sign(T))/2.
     return S0

def temp(z):
     # undercooling temperature as a function of z
     return (Tf-Ts + Tz0*z_sc*z)/(Tm-Tf)

def L(phi):
     # viscosity (coefficient on dw/dz)
     return gamma*(4./3. + 1/phi)  

def Ke(phi,S):
     # effective thermal conductivity
     return (Kw**(phi*(1-S)))*(Ks**(1-phi))/(Ki**(1-phi*S)) 

def Pi(phi):
     # sediment yield stress as a function of phi
     return alpha*((1-phi)**3)/phi**2

def Sigma(wi,phi,S):
     # scaled deviatoric ice stress
     return beta*phi*S*Dx(wi,0)

def Q(phi,S):
     # jump in heat flux across ice lenses
     return (1-Ke(phi,S))*Tz_sc

def N(phi,w):
     # effective pressure of unfrozen material
     return Pi(phi)  - (1-phi)*L(phi)*Dx(w,0) 

def q(w,wi,phi,S):
     # water flux
     return (1-phi*S)*(wi-w) 

