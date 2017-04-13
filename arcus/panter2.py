import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan
import matplotlib.pyplot as plt

import traces.analyses as anal
import traces.surfaces as surf
import traces.transformations as tran
import traces.sources as sources
import traces.surfaces as surf
import traces.grating as grat

import pdb



def traceSource(N):
    """Trace rays from PANTER point source to position of
    SPO optics
    """
    #Set up subannulus
    r0 = 737.
    F = 12e3
    L = 4*F*.605/r0
    tg = .25*np.arctan(r0/F)
    r1 = r0 + L*tg
    dphi = 50./r0/2
    rays = sources.subannulus(r0,r1,dphi,N)

    #Shift annulus to zero x coordinate
    tran.transform(rays,r0,0,0,0,0,0)
    #Set ray cosines
    srcdist = 119471.
    raydist = sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    l = rays[1]/raydist
    m = rays[2]/raydist
    n = -sqrt(1. - l**2 - m**2)
    rays = [rays[0],rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]

    #Go to optical axis
    tran.transform(rays,-r0,0,0,0,0,0)
    
    return rays

def placeSPO(rays):
    """Place an SPO pair, presuming you are already at the intersection
    node
    """
    #Place primary
    surf.spoPrimary(rays,737.,12.e3)
    tran.reflect(rays)
    #Place secondary
    surf.spoSecondary(rays,737.,12.e3)
    tran.reflect(rays)

    return

def findGratingPosition(N,hubdist=11832.911,order=1,wave=4.4,disp=0.):
    """Place the SPO pair, find the focus, and then go back
    up to grating placement
    """
    #Set up SPO
    rays = traceSource(N)
    placeSPO(rays)
    #SPO intersection plane is global coordinate system

    #Go to focus
    gratc = [tran.tr.identity_matrix()]*4
    surf.focusI(rays,coords=gratc)
    #Go back up to grating point
    tran.transform(rays,0,0,np.mean(-11828.856*rays[6]),0,0,0,coords=gratc)
    surf.flat(rays)
    
    #Place grating
    #Get to XY centroid of beam, now at center of grating
    tran.transform(rays,np.mean(rays[1]),0,0,0,0,0,coords=gratc)
    gratc2 = np.copy(gratc)
    #Rotate to proper incidence angle
    tran.steerX(rays,coords=gratc2)
    tran.transform(rays,0,0,0,0,pi/2-1.5*pi/180,0,coords=gratc2)
    surf.flat(rays)
    tran.transform(rays,0,0,0,0,0,pi/2,coords=gratc2)
    #Go to hub and diffract
    #Add yaw
    yaw = grat.blazeYaw(1.5*pi/180,2.4,3,160.)
    tran.transform(rays,0,0,0,0,0,yaw,coords=gratc2)
    tran.transform(rays,0,-hubdist+disp,0,0,0,0,coords=gratc2)
    tran.reflect(rays)
    tran.radgrat(rays,160./hubdist,order,wave)
    pdb.set_trace()
    #Go back to reference frame of grating
    rays = tran.applyT(rays,gratc2,inverse=True) #Back to global
    #rays = tran.applyT(rays,gratc) #forward to grating

    #Go to focus
    focusc = [tran.tr.identity_matrix()]*4
    surf.focusY(rays,coords=focusc)
    #Get rid of mean X and Y
    tran.transform(rays,np.mean(rays[1]),np.mean(rays[2]),0,0,0,0,\
                   coords=focusc)

    return rays,[gratc[1][i][-1] for i in range(3)],\
           [focusc[1][i][-1] for i in range(3)]
