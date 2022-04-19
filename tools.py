import numpy as np
from PyXFocus.surfaces import focusI
from PyXFocus.analyses import interpolateVec

def focalWavefront(rays,Nx,Ny,**kwargs):
    """
    Go to best focus of rays, and then interpolate a wavefront.
    This is to get rid of best fit spherical wavefront.
    Save ray x and y positions prior to going to best fit focus.
    """
    #Save a local copy of rays
    lrays = np.copy(rays)
    
    #Save x and y positions
    x = np.copy(lrays[1])
    y = np.copy(lrays[2])

    #Go to best focus
    focusI(lrays,nr=1.)

    #Reset x and y positions
    lrays[1] = x
    lrays[2] = y

    #Interpolate wavefront
    w = interpolateVec(lrays,0,Nx,Ny,**kwargs)
    
    return w
