import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan
import matplotlib.pyplot as plt
import traces.PyTrace as PT
import pdb,sys

def traceSector(Rin,Rout,F,N,span=pi/12,d=.605,t=.775,gap=25.,\
                inc=1.5*pi/180,l=95.):
    """Trace Arcus sector
    """
    #Trace parameters
    R = np.arange(Rin,Rout,t) #Vector of shell radii
    tg = .25*np.arctan((R+d/2)/F) #Primary graze angles
    L = d/tan(tg) #Primary length
    #L = 4*F*d/R #Vector of mirror lengths
    M = np.size(R) #Number of shells

    #Trace rays through SPO modules
    traceSPO(R,L,F,N,M,span=span,d=d,t=t)

    #Determine outermost radius of grating array
    PT.transform(0,0,-(L.max()+gap+95.),0,0,0)
    PT.flat()
    outerrad = np.max(sqrt(PT.x**2+PT.y**2))
    hubdist = sqrt(outerrad**2 + (1e4-(L.max()+gap+95.))**2)
    angle = np.arctan(outerrad/(1e4-(L.max()+gap+95.)))

    #Trace grating array
    gratArray(outerrad,hubdist,angle,inc,l=l)
    
    return

def traceSPO(R,L,F,N,M,span=pi/12,d=.605,t=.775):
    """Trace SPO surfaces sequentially. Collect rays from
    each SPO shell and set them to the PT rays at the end.
    Start at the inner radius, use the wafer and pore thicknesses
    to vignette and compute the next radius, loop while
    radius is less than Rout.
    """
    #Ray bookkeeping arrays
    tx = np.zeros(M*N)
    ty = np.zeros(M*N)
    tz = np.zeros(M*N)
    tl = np.zeros(M*N)
    tm = np.zeros(M*N)
    tn = np.zeros(M*N)
    tux = np.zeros(M*N)
    tuy = np.zeros(M*N)
    tuz = np.zeros(M*N)

    #Loop through shell radii and collect rays
    for i in range(M):
        #Set up source annulus
        PT.subannulus(R[i],R[i]+d,span,N)
        #Transform rays to be above xy plane
        PT.n = -PT.n
        PT.transform(0,0,-100.,0,0,0)
        #Trace to primary
        PT.spoPrimary(R[i],F)
        PT.reflect()
        #Vignette
        ind  = np.logical_and(PT.z<=L[i],PT.z>=0.)
        if np.sum(ind) < N:
            pdb.set_trace()
        PT.vignette(np.logical_and(PT.z<=L[i],PT.z>=0.))
        #Trace to secondary
        PT.spoSecondary(R[i],F)
        PT.reflect()
        #Vignette
        ind  = np.logical_and(PT.z<=0.,PT.z>=-L[i])
        if np.sum(ind) < N:
            pdb.set_trace()
        PT.vignette(ind)
        #Collect rays
        try:
            tx[i*N:(i+1)*N] = PT.x
            ty[i*N:(i+1)*N] = PT.y
            tz[i*N:(i+1)*N] = PT.z
            tl[i*N:(i+1)*N] = PT.l
            tm[i*N:(i+1)*N] = PT.m
            tn[i*N:(i+1)*N] = PT.n
            tux[i*N:(i+1)*N] = PT.ux
            tuy[i*N:(i+1)*N] = PT.uy
            tuz[i*N:(i+1)*N] = PT.uz
        except:
            pdb.set_trace()

    #Set to PT rays
    try:
        PT.x = tx
        PT.y = ty
        PT.z = tz
        PT.l = tl
        PT.m = tm
        PT.n = tn
        PT.ux = tux
        PT.uy = tuy
        PT.uz = tuz
    except:
        pdb.set_trace()

    return

def gratArray(outerrad,hubdist,angle,inc,l=95.):
    """Trace rays leaving SPO petal to the fanned grating array.
    Start with outermost radius and rotate grating array about
    the hub. Define outermost grating position by max ray radius
    at desired axial height.
    Rays have been traced to bottom of outermost grating.
    """
    #Put origin at bottom of outermost grating
    PT.transform(outerrad,0,0,0,0,0)
    #Go to proper incidence angle of grating
    PT.transform(0,0,0,0,0,-pi/2)
    PT.transform(0,0,0,-pi/2-angle+inc,0,0)
    #Go to hub
    PT.transform(0,hubdist,0,0,0,0)
    #Trace out gratings until no rays hit a grating
    #Flat
    #Indices
    #Reflect
    #Apply Grating
    #Next
    PT.flat()
    ind = np.logical_and(PT.y<-hubdist,PT.y>-l-hubdist)
    ang = l*sin(inc)/hubdist
    i = 0
    while np.sum(ind)>0:
        i = i+1
        PT.reflect(ind=ind)
        PT.radgrat(0.,160./hubdist,0,1.,ind=ind)
        PT.transform(0,0,0,ang,0,0)
        PT.flat()
        ind = np.logical_and(PT.y<-hubdist,PT.y>-l-hubdist)
        sys.stdout.write('%i \r' % i)
        sys.stdout.flush()
    pdb.set_trace()
    


    
