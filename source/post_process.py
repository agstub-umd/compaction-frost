# interpolation function for converting dolfinx functions
# to numpy arrays for plotting
import numpy as np
from dolfinx.fem import Expression, Function, FunctionSpace
from params import nz
from scipy.interpolate import griddata


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
    F = griddata(points, values, points_i, method='linear',fill_value=0.0)    
    return points_i,F

