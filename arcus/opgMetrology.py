import numpy as np
import matplotlib.pyplot as plt
import pdb
import reconstruct

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import traces.conicsolve as conic

import utilities.imaging.man as man

def littrow():
    """
    Trace rectangular beam into OPG in Littrow.
    400 nm groove period
    8.4 m groove convergence
    """
    #Set up beam
    rays = sources.rectArray(50.,50.,1e2)

    #Trace to Littrow and diffract
    tran.transform(rays,0,0,0,0,52.28*np.pi/180,0)
    surf.flat(rays)
    tran.transform(rays,0,8400.,0,0,0,0)
    tran.radgrat(rays,400./8400.,1,632.8)

    #Steer out beam tilt
    tran.steerX(rays)
    tran.steerY(rays)

    #Bring rays to common plane
    surf.flat(rays)

    #Interpolate slopes to regular grid
    y,dx,dy = anal.interpolateVec(rays,5,200,200)
    x,dx,dy = anal.interpolateVec(rays,4,200,200)

    #Prepare arrays for integration
    #Reconstruct requires nans be replaced by 100
    #Fortran functions require arrays to be packed in
    #fortran contiguous mode
    x = man.padRect(x)
    y = man.padRect(y)
    phase = np.zeros(np.shape(x),order='F')
    phase[np.isnan(x)] = 100.
    x[np.isnan(x)] = 100.
    y[np.isnan(y)] = 100.
    y = np.array(y,order='F')
    x = np.array(x,order='F')

    #Reconstruct and remove border    
    phase = reconstruct.reconstruct(x,y,1e-12,dx,phase)
    phase[phase==100] = np.nan
    x[x==100] = np.nan
    y[y==100] = np.nan

    return man.stripnans(phase),man.stripnans(x),man.stripnans(y)
