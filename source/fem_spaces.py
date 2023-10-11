from dolfinx.fem import FunctionSpace
from ufl import FiniteElement, MixedElement


def mixed_space(domain):
    P1 = FiniteElement('P',domain.ufl_cell(),1)  
    element = MixedElement([P1, P1, P1])
    V = FunctionSpace(domain,element)   
    return V    

def vel_space(domain):
    P1 = FiniteElement('P',domain.ufl_cell(),1)   
    element = P1*P1
    V = FunctionSpace(domain,element)    
    return V    