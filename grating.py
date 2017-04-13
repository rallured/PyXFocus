import numpy as np
from numpy import pi,sqrt,sin,cos,tan,exp
import matplotlib.pyplot as plt
import pdb

def blazeYaw(inc,wave,m,d):
    """Determine required yaw angle for blaze wavelength
    at order with groove period d after setting incidence
    angle inc."""
    return np.arcsin(m*wave/2/d/cos(inc))

def blazeAngle(inc,wave,m,d):
    """Determine the blaze angle for wavelength wave at order m
    with an incidence angle inc and groove period d
    """
    psi = blazeYaw(inc,wave,m,d)
    beta1 = cos(inc)*cos(psi)
    alpha1 = cos(inc)*sin(psi)-m*wave/d
    return np.arcsin(alpha1/cos(np.arcsin(beta1)))

def litBetaAlpha(inc,wave,m,d):
    """Determine the cone angle (gamma) at the Littrow condition"""
    psi = blazeYaw(inc,wave,m,d)
    beta1 = cos(inc)*cos(psi)
    alpha1 = cos(inc)*sin(psi)-m*wave/d
    return beta1,alpha1

def blazeAngle2(inc,wave,m,d):
    a = m*wave/2/d
    return np.arcsin(a/sqrt(1-cos(inc)**2+a**2))
    
def eta(phi0,theta0):
    return np.arcsin(np.cos(phi0)*np.cos(theta0))

def yaw(phi0,theta0):
    return np.arctan(np.sin(theta0)/np.tan(phi0))

def radialApproxEffect(hubdist1,hubdist2,width,length):
    """
    Determine effect of wrong radial approximation for X-ray test
    """
    #Grating coordinates
    x,y = np.meshgrid(np.linspace(-width,width,1000),\
                      np.linspace(-length,length,1000))
    y1 = y + hubdist1
    y2 = y + hubdist2

    #Convert to period and yaw angle
    period1 = np.sqrt(x**2+y1**2)/hubdist1*160. #nm
    period2 = np.sqrt(x**2+y2**2)/hubdist2*160. #nm
    yaw = blazeYaw(1.5*np.pi/180,2.4,3,160.)
    yaw1 = np.pi/2 - np.arctan(x/y1) + yaw
    yaw2 = np.pi/2 - np.arctan(x/y2) + yaw

    #Determine alpha and beta
    beta0,alpha0 = litBetaAlpha(1.5*np.pi/180,2.4,3,160.)
    alpha1 = alpha0 + 3*2.4/period1*np.sin(yaw1)
    alpha2 = alpha0 + 3*2.4/period2*np.sin(yaw2)
    beta1 = beta0 + (3*2.4/period1)*np.cos(yaw1)
    beta2 = beta0 + (3*2.4/period2)*np.cos(yaw2)

    #Determine spot shifts
    x1 = hubdist2*(alpha1/beta1)
    x2 = hubdist2*(alpha2/beta2)
    

    pdb.set_trace()
    
    return x1,x2
