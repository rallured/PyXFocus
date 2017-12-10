import numpy as np
import matplotlib.pyplot as plt
import pdb
import PyXFocus.reconstruct as reconstruct

import PyXFocus.surfaces as surf
import PyXFocus.transformations as tran
import PyXFocus.analyses as anal
import PyXFocus.sources as sources
import PyXFocus.conicsolve as conic

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

    #Get wavefront
    r = anal.wavefront(rays,200,200)

    return r

def testwave():
    """
    Verify wavefront reconstruction is occurring properly
    Coefficients can be altered in zernsurf and then resulting
    wavefront can be examined for the appropriate structure.
    """
    #Set up Zernike wave
    rays = sources.circularbeam(1.,10000)
    surf.zernsurf(rays,[0,0,0,.00,.001],1.)
    tran.reflect(rays)

    r = anal.wavefront(rays,200,200)
    return r
