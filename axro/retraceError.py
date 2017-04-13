import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.optimize as opt

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import traces.conicsolve as conic

pcoeff,pax,paz = np.genfromtxt('/home/rallured/Dropbox/AXRO/'
                        'Alignment/CoarseAlignment/150615_OP1S09Coeffs.txt')
pcoeff = pcoeff/1000.
primc = [pcoeff,pax,pax]

def tracePrimary(primCoeffs=None,primalign=np.zeros(6)):
    """
    Trace rays from focus to primary, off retroreflector,
    then back to focus. Return spot centroids.
    """
    #Set up source
    primfoc = conic.primfocus(220.,8400.)
    r1 = conic.primrad(8500.,220.,8400.)
    rays = sources.subannulus(220.,r1,100./220,100000,zhat=1.)
    tran.pointTo(rays,0.,0.,-primfoc,reverse=1.)
    theta = np.arctan2(rays[2],rays[1])

    #Trace to primary
    tran.transform(rays,*primalign)
    tran.transform(rays,0.,0,-8400.,0,0,0)
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    tran.transform(rays,0,0,8400.,0,0,0)
    tran.itransform(rays,*primalign)
    tran.reflect(rays)

    #Reflect and come back
    tran.transform(rays,0,0,400.,0,0,0)
    surf.flat(rays)
    tran.reflect(rays)
    tran.transform(rays,0,0,-400.,0,0,0)

    #Trace to primary
    tran.transform(rays,*primalign)
    tran.transform(rays,0.,0,-8400.,0,0,0)
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    ind = np.logical_and(rays[3]>8400.,rays[3]<8500.)
    tran.vignette(rays,ind=ind)
    tran.transform(rays,0,0,8400.,0,0,0)
    tran.itransform(rays,*primalign)
    tran.reflect(rays)

    #Go to primary focus
    tran.transform(rays,0,0,-primfoc,0,0,0)
    surf.flat(rays)

    return rays,theta
