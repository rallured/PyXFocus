import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pdb,os,pyfits
import utilities.imaging.fitting as fit
import traces.sources as sources
import traces.transformations as tran
import traces.surfaces as surf
import traces.analyses as anal

def distortTrace(cone1,sag1,roc1,cone2,sag2,roc2,despace=0.,secondaryTilt=0.,\
                 nominal=True):
    """
    Trace rays through a Wolter mirror pair with low
    order distortions. Distortions can be applied to
    either primary or secondary or both.
    axial sag, azimuthal sag, and cone angle are included.
    
    """
    #Define ray subannulus
    r1 = surf.con.primrad(8600.,1000.,8400.)
    ang = 260./1000. #arc length over radius is angular extent
    rays = sources.subannulus(1000.,r1,ang,10**3)
    tran.transform(rays,0,0,0,np.   pi,0,0) #Point in -z
    tran.transform(rays,0,0,-10000,0,0,0) #Converge from above

    #Trace to primary
    surf.primaryLL(rays,1000.,8400.,8600,8400,ang,[cone1,sag1,roc1],\
                   [1,2,0],[0,0,2])
    #Vignette rays missing
    ind = np.logical_and(rays[3]<8600.,rays[3]>8400.)
    rays = tran.vignette(rays,ind)
    numin = float(len(rays[1]))
    #Reflect
    tran.reflect(rays)
    
    #Apply secondary misalignment
    tran.transform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    tran.transform(rays,0,0,despace,0,secondaryTilt,0)
    tran.itransform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    #Trace to secondary
    surf.secondaryLL(rays,1000.,8400.,8400.,8200.,ang,[cone2,sag2,roc2],\
                     [1,2,0],[0,0,2])
    #Vignette rays missing
    ind = np.logical_and(rays[3]<8400.,rays[3]>8200.)
    rays = tran.vignette(rays,ind)
    numout = float(len(rays[1]))
    #Reflect
    tran.reflect(rays)
    #Reverse secondary misalignment
    tran.transform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    tran.itransform(rays,0,0,despace,0,secondaryTilt,0)
    tran.itransform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)

    #Go to focus
    if nominal is True:
        surf.flat(rays)
    else:
        surf.focusI(rays)

    #Return merit function
    return rays#anal.hpd(rays)/8400.,numout/numin

def sensitivityPlots(dof,mag):
    """
    Mark the RMS spread and centroid shift of focus
    as function of DoF
    """
    align = np.repeat(0.,6)
    misalign = np.linspace(0,mag,100)
    cx = []
    rms = []
    for m in misalign:
        align[dof] = m
        rays = distortTrace(*align)
        cx.append(anal.centroid(rays)[0])
        rms.append(anal.rmsCentroid(rays))
        
    return cx,rms

