import numpy as np
import matplotlib.pyplot as plt

def torusF(x,y,z,rin,rout):
    """Return the Torus surface function over x,y,z
    with parameters rin and rout
    Equation taken from
    http://www.emeyex.com/site/projects/raytorus.pdf
    """
    return (x**2+y**z+z**2-rin**2-rout**2)**2 + \
           4*rout**2*(z**2-rin**2)

def torusF2(x,y,z,rin,rout):
    """Return the Torus surface function with the origin
    tangent to the outer surface.
    """
    return ((z+rin+rout)**2+y**2+x**2+rout**2-rin**2)**2 - \
           4*rout**2*((z+rin+rout)**2+y**2)

def torusGrad(x,y,z,rin,rout):
    """Return the derivatives of the Torus surface function
    with origin tangent to outer surface."""
    dfdz = 4*(rin+rout+z)*(2*rin*(rout+z)+2*rout*z+x**2+y**2+z**2)
    dfdy = 4*y*(2*rin*(rout+z)+2*rout*z+x**2+y**2+z**2)
    dfdx = 4*x*(-rin**2+(rin+rout+z)**2+rout**2+x**2+y**2)
    return dfdx,dfdy,dfdz
    
def constructQuartic(p,d,rin,rout):
    """Construct the quartic equation for intersection
    of a line with a torus. p is position, d is direction
    of line."""
    #Dot parameters
    alpha = np.dot(d,d)
    beta = 2*np.dot(p,d)
    gamma = np.dot(p,p)-rin**2-rout**2
    #Quartic parameters
    a4 = alpha**2
    a3 = 2*alpha*beta
    a2 = beta**2 + 2*alpha*gamma + 4*rout**2*d[2]**2
    a1 = 2*beta*gamma + 8*rout**2*p[2]*d[2]
    a0 = gamma**2 + 4*rout**2*p[2]**2 - 4*rout**2*rin**2
    return np.array([a4,a3,a2,a1,a0])

##def findRoots(a):
##    
