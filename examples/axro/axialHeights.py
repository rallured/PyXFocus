import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.interpolate as inte
import astropy.io.fits as pyfits

import PyXFocus.surfaces as surf
import PyXFocus.transformations as tran
import PyXFocus.analyses as anal
import PyXFocus.sources as sources
import PyXFocus.conicsolve as conic

import inducedPolarization as pol

cons = pol.readNK('/home/rallured/Dropbox/inducedPol/Ir_llnl_cxro.nk')

def showZetaEffect(zeta=1.,col='b'):
    """
    Trace zeta > 1 and zeta < 1 for a single shell.
    Demonstrate that zeta > 1 causes rays to spread out
    thus needing a longer secondary
    """
    #Set up source
    rays = sources.annulus(220.,221.,10000)
    rays[6] = rays[6]*-1
    tran.transform(rays,0,0,-8400.,0,0,0)

    #Set up figure
    plt.figure('ZetaEffect')

    #Trace
    surf.wolterprimary(rays,220.,8400.,psi=zeta)
    plt.plot(np.sqrt(rays[1]**2+rays[2]**2),rays[3],col+'.')
    tran.reflect(rays)
    surf.woltersecondary(rays,220.,8400.,psi=zeta)
    plt.plot(np.sqrt(rays[1]**2+rays[2]**2),rays[3],col+'.')
    

def traceZeta(pmin,R0=220.,Z0=1e4,psi=1.,offaxis=0.,L=100.,az=100.):
    """
    Set initial aperture based on height of bottom of mirror
    above node (pmin). Assume 100 mm long mirrors.
    Then, trace to secondary and mark smin and smax for
    where the rays strike.
    Then trace out off axis field positions and determine
    RMS and HPD vs angle.
    """
    #Set up aperture
    a0 = conic.primrad(pmin,R0,Z0)
    a1 = conic.primrad(pmin+L,R0,Z0) 
    rays = sources.subannulus(a0,a1,az/220,1e4)
    tran.transform(rays,0,0,-Z0,0,0,0)

    #Trace to primary and add off-axis angle
    surf.wsPrimary(rays,R0,Z0,psi)
    rays[4] = rays[4]+np.sin(offaxis)
    rays[6] = -np.sqrt(1.-rays[4]**2)
    tran.reflect(rays)

    #Trace to secondary
    surf.wsSecondary(rays,R0,Z0,psi)
    tran.reflect(rays)
    smax = np.nanmax(rays[3])
    smin = np.nanmin(rays[3])

    #Go to focus
    f = surf.focusI(rays)

    #Compute merit functions
    hpd = anal.hpd(rays)
    rms = anal.rmsCentroid(rays)

    return smin,smax,f,hpd,rms

def onaxisMerit(pmin,smax,R0=220.,Z0=1e4,L=200.,psi=1.):
    """
    Trace rays for a given shell and determine amount of rays vignetted
    """
    #Set up aperture
    a0 = conic.primrad(pmin,R0,Z0)
    a1 = conic.primrad(pmin+L,R0,Z0) 
    rays = sources.annulus(a0,a1,1e3)
    tran.transform(rays,0,0,-Z0,0,0,0)

    #Trace to primary
    surf.wsPrimary(rays,R0,Z0,psi)
    tran.reflect(rays)

    #Trace to secondary
    surf.wsSecondary(rays,R0,Z0,psi)

    #Determine fraction of rays vignetted
    return np.sum(np.logical_or(rays[3]<smax-L,rays[3]>smax))
    
def determineZeta(pmin,smax,R0=220.,Z0=1e4,L=200.):
    """
    Determine best value of zeta to put rays at the needed
    axial position.
    This will be maximum zeta such that the vignetting is zero
    """
    prange = np.linspace(.1,1.,100)
    v = np.array([onaxisMerit(pmin,smax,R0=R0,Z0=Z0,L=L,psi=p) for p in prange])
    try:
        return max(prange[v==0.]),0.
    except:
        return prange[np.argmin(v)],v.min()/1e3

def defineRx(N=3,L=200.,nodegap=25.,rnodes=False):
    """
    Determine node radii and axial positions
    Divide nodes into N sections.
    For each section, define Pmin and Smax.
    For each node, determine optimal zeta.
    Finally, trace each node, multiply by R**2,
    and determine weighted resolution and effective area.
    Determine off-axis resolution and effective area including vignetting.
    """
    #Determine node positions
    rad = np.arange(200.,1501.,2.)
    z = np.sqrt(1e4**2-rad**2)

    #Split into N sections and determine pmin and smax
    smax = np.zeros(len(rad))
    pmin = np.zeros(len(rad))
    smaxu = []
    pminu = []
    for i in range(N):
        ind = np.logical_and(rad>=200.+(1300./N)*i,rad<=200.+(1300./N)*(i+1))
        smax[ind] = np.min(z[ind])-nodegap
        pmin[ind] = 2*np.max(z[ind])-smax[ind]
        smaxu.append(np.min(z[ind])-nodegap)
        pminu.append(2*np.max(z[ind])-smaxu[i])

    #Determine proper zeta for each node
    zeta = np.array([determineZeta(pmin[i],smax[i],R0=rad[i],Z0=z[i],L=L) \
            for i in range(len(rad))])

    #Create an interpolation function for zeta for each section
    fun = []
    for i in range(N):
        ind = np.logical_and(rad>=200.+(1300./N)*i,rad<=200.+(1300./N)*(i+1))
        fun.append(inte.interp1d(rad[ind],np.transpose(zeta)[0][ind],\
                                 kind='linear',bounds_error=False,\
                                 fill_value="extrapolate"))

    if rnodes is True:
        return smax,pmin,np.transpose(zeta)[0],np.transpose(zeta)[1]

    return smaxu,pminu,fun



def traceXRS(smax,pmin,fun,L=200.,nodegap=25.,Nshell=1e3,energy=1000.,\
             rough=1.,offaxis=0.,rrays=False):
    """
    Using output from defineRx, trace the nodes in a Lynx design.
    Provide number of sections, and the zeta as a function of radius
    interpolation function.
    Calculate node radii iteratively, providing appropriate gaps between
    nodes in order to limit off-axis vignetting.
    """
    #Loop through sections and construct node positions
    N = len(smax)
    rsec = []
    for sec in range(N):
        #Establish radius vector
        rad = np.array([])
        rout = 0.
        #Compute shell gap
        gap = (pmin[sec]+L-1e4)*3e-3+0.4 #4 mm glass thickness plus vignetting
        #First node position
        rad = np.append(rad,200.+(1300./N)*sec)
        rout = conic.primrad(pmin[sec]+L,rad[-1],\
                             np.sqrt(1e4**2-rad[-1]**2))
        while rout+gap < 200.+(1300./N)*(sec+1):
            rad = np.append(rad,rout+gap)
            rout = conic.primrad(pmin[sec]+L,rad[-1],\
                                 np.sqrt(1e4**2-rad[-1]**2))
        #Add to list of node positions
        rsec.append(rad)

    #Trace through all shells, computing reflectivity and geometric area
    #for each shell
    for i in range(N):
        for r in rsec[i]:
            #Set up aperture
            z = np.sqrt(1e4**2-r**2)
            a0 = conic.primrad(pmin[i],r,z)
            a1 = conic.primrad(pmin[i]+L,r,z) 
            rays = sources.annulus(a0,a1,Nshell)
            tran.transform(rays,0,0,-z,0,0,0)

            #Set up weights
            weights = np.repeat((a1**2-a0**2) * np.pi / 100. / Nshell,Nshell)

            #Trace to primary
            psi = fun[i](r)
            surf.wsPrimary(rays,r,z,psi)
            rays[4] = rays[4]+np.sin(offaxis)
            rays[6] = -np.sqrt(1.-rays[4]**2)
            tran.reflect(rays)
            ang = anal.grazeAngle(rays)
            weights = weights*\
                      pol.computeReflectivities(ang,energy,rough,1,cons)[0]

            #Trace to secondary
            surf.wsSecondary(rays,r,z,psi)
            tran.reflect(rays)
            ang = anal.grazeAngle(rays)
            weights = weights*\
                      pol.computeReflectivities(ang,energy,rough,1,cons)[0]

            #Handle vignetting
            ind = np.logical_and(rays[3]>smax[i]-L,rays[3]<smax[i])
            if sum(ind) == 0:
                pdb.set_trace()
            rays = tran.vignette(rays,ind=ind)
            weights = weights[ind]

            #Go to focus
            try:
                surf.focusI(rays,weights=weights)
            except:
                pdb.set_trace()

            #Accumulate master rays
            try:
                mrays = [np.concatenate([mrays[ti],rays[ti]]) for ti in range(10)]
                mweights = np.concatenate([mweights,weights])
            except:
                mrays = rays
                mweights = weights

    if rrays is True:
        return mrays,mweights

    return anal.hpd(mrays,weights=mweights)/1e4*180/np.pi*60**2,\
           anal.rmsCentroid(mrays,weights=mweights)/1e4*180/np.pi*60**2

"""
Factors to consider for this project:
- Impact on off-axis resolution
- Impact on inherent area due to increased graze angle, does this
impact aperture diffraction significantly?
- Mirror thickness needs to be added
- Impact on reflectivity due to increased graze angle on one of the mirrors

What to report?
Weighted Diffraction limit - just use analytic equations
Effective area, HPD, RMSD vs. off-axis angle
Need to compare with nominal Rx with varying mirror axial heights
"""
