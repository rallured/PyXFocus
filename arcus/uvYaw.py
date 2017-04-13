import numpy as np
import matplotlib.pyplot as plt
import pdb

import traces.sources as sources
import traces.transformations as tran
import traces.surfaces as surf

#Set up incident beam trace and determine sensitivity to beam
#impact location.
#Trace nominal geometry (function of incidence angle) and
#record location of diffracted and reflected spots
#Perturb location and angle of beam impact and record
#change in spot locations

#What don't you know about geometry?
#Location and angle of source -> location and angle of beam impact
#Location of detectors (angle is very small contributor)

#Write raytrace with beam impact orientation and
#grating orientation as free parameters

def alignTrace(inc,impact,grating,detector,order=0):
    """Traces UV laser rays to grating. Beam impact misalignment
    is handled with a single coordinate transformations right after
    source definition. Grating orientation is handled with
    symmetric coordinate transformations.
    inc - nominal beam glancing angle, must be less than
        50.39 deg for 262 nm light
    impact - 6 element array giving beam impact transform
    grating - 6 element array giving grating misalignment
    """
    #Set up source with single ray, diffraction plane
    #is XZ, glancing angle from XY plane, ray starts out
    #pointing +x and -z
    rays = sources.pointsource(0.,1)
    tran.transform(rays,0,0,0,0,-np.pi/2-inc,0)
    #Perform beam impact misalignment transform, rotation first
    tran.transform(rays,*np.concatenate(((0,0,0),impact[3:])))
    tran.transform(rays,*np.concatenate((impact[:3],(0,0,0))))
    #Perform grating misalignment
    tran.transform(rays,*grating)
    #Linear grating
    surf.flat(rays)
    tran.reflect(rays)
    tran.grat(rays,160.,order,262.)
    #Reverse misalignment transformation
    tran.itransform(rays,*grating)
    #Go to detector depending on order
    if order is not 0:
        tran.transform(rays,-200.,0,0,0,0,0)
    else:
        tran.transform(rays,200.,0,0,0,0,0)
    #Trace to detector
    tran.transform(rays,0,0,0,0,-np.pi/2,0)
    tran.transform(rays,*detector)
    surf.flat(rays)
    #Return ray position
    return rays[1],rays[2]

def computeYaw(inc,impact,grating,detr,detd):
    """Traces both orders and computes yaw of grating
    Uses alignTrace
    Returns yaw angle assuming 
    """
    xr,yr = alignTrace(inc,impact,grating,detr,order=0)
    betar = 800./np.sqrt(800.**2+xr**2+yr**2)
    alphar = -yr/np.sqrt(800.**2+xr**2+yr**2)
    xd,yd = alignTrace(inc,impact,grating,detd,order=1)
    betad = -800./np.sqrt(800.**2+xd**2+yd**2)
    alphad = -yd/np.sqrt(800.**2+xd**2+yd**2)
    
    return np.arctan((alphad-alphar)/(betar-betad))*180./np.pi*60.**2


def dofSensitivity(inc,alignvector,obj='beam',dof=0):
    """Compute x,y positions of reflected and diffracted spots
    return as a function of alignvector
    """
    #Initialize misalignment vectors
    grating = np.zeros(6)
    impact = np.zeros(6)

    #Initialize output vectors
    xr = np.zeros(np.size(alignvector))
    yr = np.copy(xr)
    xd = np.copy(xr)
    yd = np.copy(xr)

    #Perform raytraces in loop
    for a in alignvector:
        #Adjust misalignments
        if obj is 'beam':
            impact[dof] = a
        else:
            grating[dof] = a

        #Perform trace and set appropriate output elements
        i = a==alignvector
        x,y = alignTrace(inc,impact,grating,order=0)
        xr[i] = x
        yr[i] = y
        x,y = alignTrace(inc,impact,grating,order=1)
        xd[i] = x
        yd[i] = y

    #Return
    return xr,yr,xd,yd

def diffractionAngle(inc):
    """Return the diffraction angle for the UV yaw system
    Input graze angle in degrees
    Output diffraction graze angle in degrees"""
    alpha0 = np.sin((90.-inc)*np.pi/180)
    alpha1 = alpha0 - 266e-9/160e-9
    dang = 90 - np.arcsin(np.abs(alpha1))*180/np.pi
    return dang
