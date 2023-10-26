# this file contains some helper functions for the compaction problem solver

import numpy as np
from params import Tf, Tm, Ts, Tz0, a, b, gamma, z_sc
from ufl import Dx


def perm(S):
      # scaled permeability as a function of ice saturation
      return (1-S)**a

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

def temp(z):
     # undercooling temperature as a function of z
     return (Tf-Ts + Tz0*z_sc*z)/(Tm-Tf)

def Lamda(phi):
     # viscosity (coefficient on dw/dz)
     return gamma*(4./3. + 1/phi)  

def Nc(phi,w):
     # deviatoric effective stress (compaction)
     return - (1-phi)*Lamda(phi)*Dx(w,0)  

def q(w,phi,S):
     # water flux relative to rigid solution
     return -(1-phi*S)*w 

def qr(V_heave,phi,S):
     # water flux for rigid solution
     # V_heave = v_i - v_s
     return (1-phi*S)*V_heave

# def Q(phi,S):
#      # jump in heat flux across ice lenses
#      return 0#(1-Ke(phi,S))*Tz_sc
