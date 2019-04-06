#This module will test the behavior of the Wolter Schwarzschild surface
#tracing. Rays are launched at a variety of radii and incidence angles to
#determine under what conditions the rays are actually traced to the mirror.

import numpy as np
import matplotlib.pyplot as plt
import PyXFocus.surfaces as surf
import PyXFocus.sources as sources
import PyXFocus.transformations as tran
import PyXFocus.analyses as anal

def testNormalPrimary(r0,z0,zeta=1.,offx=0.,dz=0.):
    """
    Launch a slit of rays toward the mirror from infinity.
    Return the 
    """
    #Set up source
    rays = sources.circularbeam(r0*2,1e4)
    tran.transform(rays,0,0,-z0+dz,0,0,0)

    #Apply offaxis angle
    rays[4] = rays[4] + offx
    rays[6] = -np.sqrt(1-rays[4]**2)

    #Trace to primary
    xi = np.copy(rays[1])
    yi = np.copy(rays[2])
    surf.wsPrimary(rays,r0,z0,zeta)

    return rays,xi,yi

def testTransversePrimary(r0,z0,zeta=1.,offz=0.,dx=0.):
    """
    Launch rays perpendicular to optical axis toward mirror.
    """
    #Set up source
    rays = sources.rectbeam(z0,r0,1e4)
    tran.transform(rays,0,0,0,0,np.pi/2,0)
    tran.transform(rays,-r0+dx,0,-z0,0,0,0)

    #Apply offaxis angle
    rays[6] = rays[6] + offz
    rays[4] = -np.sqrt(1-rays[6]**2)

    #Trace to primary
    zi = np.copy(rays[3])
    yi = np.copy(rays[2])
    surf.wsPrimary(rays,r0,z0,zeta)

    return rays,yi,zi

def testNormalSecondary(r0,z0,zeta=1.,offx=0.,dz=0.):
    """
    Launch a slit of rays toward the mirror from infinity.
    Return the 
    """
    #Set up source
    rays = sources.circularbeam(r0*2,1e4)
    tran.transform(rays,0,0,-z0+dz,0,0,0)

    #Apply offaxis angle
    rays[4] = rays[4] + offx
    rays[6] = -np.sqrt(1-rays[4]**2)

    #Trace to primary
    xi = np.copy(rays[1])
    yi = np.copy(rays[2])
    surf.wsSecondary(rays,r0,z0,zeta)

    return rays,xi,yi

def testTransverseSecondary(r0,z0,zeta=1.,offz=0.,dx=0.):
    """
    Launch rays perpendicular to optical axis toward mirror.
    """
    #Set up source
    rays = sources.rectbeam(z0,r0,1e4)
    tran.transform(rays,0,0,0,0,np.pi/2,0)
    tran.transform(rays,-r0+dx,0,-z0,0,0,0)

    #Apply offaxis angle
    rays[6] = rays[6] + offz
    rays[4] = -np.sqrt(1-rays[6]**2)

    #Trace to primary
    zi = np.copy(rays[3])
    yi = np.copy(rays[2])
    surf.wsSecondary(rays,r0,z0,zeta)

    return rays,yi,zi
