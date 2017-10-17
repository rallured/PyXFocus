import numpy as np
import matplotlib.pyplot as plt
import pdb,time
import scipy.interpolate as inte
import astropy.io.fits as pyfits
import scipy.optimize as opt

import PyXFocus.surfaces as surf
import PyXFocus.transformations as tran
import PyXFocus.analyses as anal
import PyXFocus.sources as sources
import PyXFocus.conicsolve as conic

import inducedPolarization as pol

cons = pol.readNK('/home/rallured/Dropbox/inducedPol/Ir_llnl_cxro.nk')

def showZetaEffect(zeta=1.,col='b',label=''):
    """
    Trace zeta > 1 and zeta < 1 for a single shell.
    Demonstrate that zeta > 1 causes rays to spread out
    thus needing a longer secondary
    """
    #Set up source
    rp1 = wsPrimrad(1e4+250.,220.,1e4,psi=zeta)
    rp2 = wsPrimrad(1e4+50.,220.,1e4,psi=zeta)

    rs1,zs1 = traceSingleRay(rp1,8000.,220.,1e4,zeta,rpos=True)
    rs2,zs2 = traceSingleRay(rp2,8000.,220.,1e4,zeta,rpos=True)

    plt.figure('Zeta')
    plt.plot([rp1,rp1],[1e4+500.,1e4+250.],col+'--',linewidth=.5)
    plt.plot([rp1,rs1],[1e4+250.,zs1],col+'--',linewidth=.5)
    plt.plot([rp2,rp2],[1e4+500.,1e4+50.],col+'--',linewidth=.5)
    plt.plot([rp2,rs2],[1e4+50.,zs2],col+'--',linewidth=.5)
    plt.plot([rs1,0],[zs1,0],col+'--',linewidth=.5)
    plt.plot([rs2,0],[zs2,0],col+'--',linewidth=.5)

    plt.plot([rp1,rp2],[1e4+250.,1e4+50.],col,label=label,linewidth=3)
    plt.plot([rs1,rs2],[zs1,zs2],col,linewidth=3)

    plt.xlim([215.,225])
    plt.ylim([1e4+-300,1e4+500])

    return None 

def wsPrimrad(z,r0,z0,psi=1.):
    """
    Determine radius of WS primary
    """
    #Set up ray
    ray = sources.pointsource(0.,1)
    tran.transform(ray,0,0,0,0,-np.pi/2,0)
    tran.transform(ray,-r0-2.,0,-z,0,0,0)

    surf.wsPrimary(ray,r0,z0,psi)

    return ray[1][0]

def wsSecrad(z,r0,z0,psi=1.):
    """
    Determine radius of WS primary
    """
    #Set up ray
    ray = sources.pointsource(0.,1)
    tran.transform(ray,0,0,0,0,-np.pi/2,0)
    tran.transform(ray,-r0-2.,0,-z,0,0,0)

    surf.wsSecondary(ray,r0,z0,psi)

    return ray[1][0]

def traceZeta(pmin,R0=220.,Z0=1e4,psi=1.,offaxis=0.,L=200.,az=100.,pause=False):
    """
    Set initial aperture based on height of bottom of mirror
    above node (pmin). Assume 100 mm long mirrors.
    Then, trace to secondary and mark smin and smax for
    where the rays strike.
    Then trace out off axis field positions and determine
    RMS and HPD vs angle.
    """
    #Set up aperture
    a0 = wsPrimrad(pmin,R0,Z0,psi)
    a1 = wsPrimrad(pmin+L,R0,Z0,psi) 
    rays = sources.subannulus(a0,a1,az/R0,1e4)
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
    if pause is True:
        pdb.set_trace()

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
    a0 = wsPrimrad(pmin,R0,Z0,psi)
    a1 = wsPrimrad(pmin+L,R0,Z0,psi)
    rays = sources.annulus(a0,a1,1e3)
    tran.transform(rays,0,0,-Z0,0,0,0)

    #Trace to primary
    surf.wsPrimary(rays,R0,Z0,psi)
    tran.reflect(rays)

    #Trace to secondary
    surf.wsSecondary(rays,R0,Z0,psi)

    #Determine fraction of rays vignetted
    try:
        return np.sum(np.logical_or(rays[3]<smax-L,rays[3]>smax)),np.mean(rays[3])
    except:
        pdb.set_trace()
    
def determineZeta(pmin,smax,R0=220.,Z0=1e4,L=200.):
    """
    Determine best value of zeta to put rays at the needed
    axial position.
    This will be maximum zeta such that the vignetting is zero
    """
    prange = np.linspace(.01,1.1,100)
    v = np.array([onaxisMerit(pmin,smax,R0=R0,Z0=Z0,L=L,psi=p) for p in prange])
    com = np.polyval(np.polyfit(prange,v[:,1],2),prange)
    ind = np.argmin(np.abs(com-(smax-L/2)))
    print R0
    return prange[ind]#,v[ind,0]

def defineRx(N=3,L=200.,nodegap=25.,rnodes=False,rzeta=False,\
             zeta=None,smaxu=None,pminu=None,snodegap=None,rad=None):
    """
    Determine node radii and axial positions
    Divide nodes into N sections.
    For each section, define Pmin and Smax.
    For each node, determine optimal zeta.
    Finally, trace each node, multiply by R**2,
    and determine weighted resolution and effective area.
    Determine off-axis resolution and effective area including vignetting.
    """
    if zeta is None:
        #Determine node positions
        rad = np.arange(200.,1501.,5.)
        z = np.sqrt(1e4**2-rad**2)

        #Split into N sections and determine pmin and smax
        smax = np.zeros(len(rad))
        pmin = np.zeros(len(rad))
        smaxu = []
        pminu = []
        snodegap = []
        for i in range(N):
            ind = np.logical_and(rad>=200.+(1300./N)*i,\
                                 rad<=200.+(1300./N)*(i+1))
            zm = np.max(z[ind])
            zmin = np.min(z[ind])
            if 2*(zm-zmin) > nodegap:
                print 'Invalid nodegap in section ' + str(i+1)
                nodegap = 2*(zm-zmin)
                print 'Increased to ' + str(nodegap)
            smax[ind] = zm-nodegap/2.
            pmin[ind] = zm+nodegap/2.
            smaxu.append(zm-nodegap/2.)
            pminu.append(zm+nodegap/2.)
            snodegap.append(nodegap)

        #Determine proper zeta for each node
        try:
            zeta = np.array([determineZeta(pmin[i],smax[i],R0=rad[i],Z0=z[i],L=L) \
                for i in range(len(rad))])
        except:
            pdb.set_trace()

        if rzeta is True:
            return rad,zeta,smaxu,pminu,snodegap

    #Switch this to a quadratic function fit
    #Create an interpolation function for zeta for each section
    fun = []
    for i in range(N):
        ind = np.logical_and(rad>=200.+(1300./N)*i,rad<200.+(1300./N)*(i+1))
##        fun.append(inte.interp1d(rad[ind],zeta[ind],\
##                                 kind='linear',bounds_error=False,\
##                                 fill_value="extrapolate"))
        fit = np.polyfit(rad[ind],zeta[ind],2)
        fun.append(fit)

    if rnodes is True:
        return smax,pmin,zeta

    return smaxu,pminu,fun,snodegap


def tracePerfectXRS(L=200.,nodegap=50.,Nshell=1e3,energy=1000.,\
                    rough=1.,offaxis=0.,rrays=False,rnodes=False):
    """
    Trace rays through a perfect Lynx design where all the shells
    are on the spherical principle surface, with zeta equal to unity.
    """
    #Construct node positions
    rad = np.array([200.])
    z = np.sqrt(1e4**2-rad[-1]**2)
    #Gap is projected radial shell plus thickness plus vignetting gap
    rout = wsPrimrad(z+L+nodegap/2.,rad[-1],z)
    gap = L*3e-3 + 0.4
    #Establish radius vector
    rad = np.array([])
    for sec in range(3):
        rout = 0.
        #First node position
        rad = np.append(rad,200.+(1300./3)*sec)
        #rguess = np.linspace(
        while rout+gap < 200.+(1300./3)*(sec+1)-10.: #Need to adjust this condition
            #Compute parameters for current node
            z = np.sqrt(1e4**2-rad[-1]**2)
            rout = wsPrimrad(z+L+nodegap/2.,rad[-1],z)
            rad = np.append(rad,rout+gap)
        rad = rad[:-1]

    if rnodes is True:
        return rad

    #Use radial nodes and trace Lynx, keeping track of effective area
    previousrho = 0.
    for r in rad:
        #Set up aperture
        z = np.sqrt(1e4**2-r**2)
        a0 = wsPrimrad(z+nodegap/2.,r,z)
        a1 = wsPrimrad(z+nodegap/2.+L,r,z) 
        rays = sources.annulus(a0,a1,Nshell)
        tran.transform(rays,0,0,-z,0,0,0)

        #Set up weights (cm^2)
        weights = np.repeat((a1**2-a0**2) * np.pi / 100. / Nshell,Nshell)

        #Trace to primary
        surf.wsPrimary(rays,r,z,1.)
        rays[4] = rays[4]+np.sin(offaxis)
        rays[6] = -np.sqrt(1.-rays[4]**2)
        tran.reflect(rays)
        ang = anal.grazeAngle(rays)
        weights = weights*\
                  pol.computeReflectivities(ang,energy,rough,1,cons)[0]

        #Trace to secondary
        surf.wsSecondary(rays,r,z,1.)
        tran.reflect(rays)
        ang = anal.grazeAngle(rays)
        weights = weights*\
                  pol.computeReflectivities(ang,energy,rough,1,cons)[0]

        #Handle vignetting
        ind = np.logical_and(rays[3]>z-nodegap/2.-L,rays[3]<z-nodegap/2.)
        if sum(ind) == 0:
            pdb.set_trace()
        rays = tran.vignette(rays,ind=ind)
        weights = weights[ind]

        #Go to exit aperture and confirm rays don't - EDIT HERE!
        #hit back of previous shell
        tran.transform(rays,0,0,z-nodegap/2-L,0,0,0)
        surf.flat(rays)
        rho = np.sqrt(rays[1]**2+rays[2]**2)
        ind = rho > previousrho
        #if np.sum(ind)==0:
            #pdb.set_trace()
        if np.sum(~ind) > 100:
            print '%i rays hit back of shell' % np.sum(~ind)
            print r
            #pdb.set_trace()
        rays = tran.vignette(rays,ind=ind)
        weights = weights[ind]
        previousrho = wsSecrad(z-nodegap/2-L,r,z)+.4

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
           anal.rmsCentroid(mrays,weights=mweights)/1e4*180/np.pi*60**2,\
           np.sum(mweights)

def traceSingleRay(rtrace,zexit,R0,Z0,zeta,rpos=False):
    """
    Trace single ray through WS shell and return
    ray at a given exit aperture z position
    """
    #Set up ray
    ray = sources.pointsource(0.,1)
    ray[6] = -ray[6]
    ray[1][0] = rtrace
    tran.transform(ray,0,0,-Z0,0,0,0)

    #Trace to shell
    surf.wsPrimary(ray,R0,Z0,zeta)
    tran.reflect(ray)
    surf.wsSecondary(ray,R0,Z0,zeta)
    tran.reflect(ray)
    if rpos is True:
        return ray[1][0],ray[3][0]

    #Go to exit aperture
    tran.transform(ray,0,0,zexit,0,0,0)
    surf.flat(ray)
    
    return ray

def findNextNode(rsmindesired,rguess,smax,pmin,fun,L=200.):
    """
    Use numerical approach to find position of next node
    Given pmin and smax, algorithm will trace extreme ray
    and adjust node radius until it converges to desired
    radius rsmindesired.
    """
    #Set up optimization function - this might not be behaving
    rsminout = lambda r: np.abs(rsmindesired - \
                    traceSingleRay(wsPrimrad(pmin,r,np.sqrt(1e4**2-r**2),\
                                        psi=np.max([np.polyval(fun,r),.01])),\
                                        smax-L,r,np.sqrt(1e4**2-r**2),\
                                        np.max([np.polyval(fun,r),.01]))[1][0])
    #Optimize node radius
    res = opt.minimize(rsminout,rguess,method='Nelder-Mead')

    return res['x'][0],res['nfev']


def defineNodePositions(smax,pmin,fun,nodegap,L=200.):
    """
    Loop through sections and determine optimal placement
    of nodes in order to prevent vignetting/collisions
    """
    #Loop through sections and construct node positions
    N = len(smax)
    rsec = []
    rext = []
    gap = L*3e-3+0.4 #.4 mm glass thickness plus vignetting
    for sec in range(N):
        #Establish radius vector
        rad = np.array([])
        rout = 0.
        #First node position
        rad = np.append(rad,200.+(1300./N)*sec)
        #rguess = np.linspace(
        while rout+gap < 200.+(1300./N)*(sec+1)-10.: #Need to adjust this condition
            #Compute parameters for current node
            z = np.sqrt(1e4**2-rad[-1]**2)
            psi = np.max([np.polyval(fun[sec],rad[-1]),.01])
            #Compute next node radius
            rsmin = wsSecrad(smax[sec]-L,rad[-1],z,\
                             psi=psi)+gap
            nrad,fev = findNextNode(rsmin,rad[-1],\
                                    smax[sec],pmin[sec],fun[sec],L=L)
            rad = np.append(rad,nrad)#rout+gap)
            rout = wsPrimrad(pmin[sec]+L,rad[-1],\
                                 np.sqrt(1e4**2-rad[-1]**2),\
                             psi=np.max([np.polyval(fun[sec],rad[-1]),.01]))
        rad = rad[:-1]

        #Add to list of node positions
        rsec.append(rad[:-1])
        
    return rsec#,rext

def traceXRS(rsec,smax,pmin,fun,nodegap,L=200.,Nshell=1e3,energy=1000.,\
             rough=1.,offaxis=0.,rrays=False):
    """
    Using output from defineRx, trace the nodes in a Lynx design.
    Provide number of sections, and the zeta as a function of radius
    interpolation function.
    Calculate node radii iteratively, providing appropriate gaps between
    nodes in order to limit off-axis vignetting.
    """
    gap = L*3e-3+0.4 #.4 mm glass thickness plus vignetting
    #Trace through all shells, computing reflectivity and geometric area
    #for each shell
    previousrho = []
    for i in range(len(smax)):
        #Variable to store radial position of bottom of previous
        #secondary mirror
        previousrho.append(0.)
        for r in rsec[i]:
            #Determine zeta for this shell...must be at least .01
            psi = np.polyval(fun[i],r)
            psi = np.max([.01,psi])
            #Set up aperture
            z = np.sqrt(1e4**2-r**2)
            a0 = wsPrimrad(pmin[i],r,z,psi=psi)
            a1 = wsPrimrad(pmin[i]+L,r,z,psi=psi) 
            rays = sources.annulus(a0,a1,Nshell)
            tran.transform(rays,0,0,-z,0,0,0)

            #Set up weights (cm^2)
            weights = np.repeat((a1**2-a0**2) * np.pi / 100. / Nshell,Nshell)

            #Trace to primary
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

            #Go to exit aperture and confirm rays don't
            #hit back of previous shell
            tran.transform(rays,0,0,smax[i]-L,0,0,0)
            surf.flat(rays)
            rho = np.sqrt(rays[1]**2+rays[2]**2)
            ind = rho > previousrho[-1]
            #if np.sum(ind)==0:
                #pdb.set_trace()
            if np.sum(~ind) > 100:
                print '%i rays hit back of shell' % np.sum(~ind)
                print r,psi
                #pdb.set_trace()
            rays = tran.vignette(rays,ind=ind)
            weights = weights[ind]
            previousrho.append(wsSecrad(smax[i]-L,r,z,psi=psi)+.4)

            #Go to focus
            if rrays is False:
                tran.transform(rays,0,0,-smax[i]+L,0,0,0)
                surf.flat(rays)

            #Accumulate master rays
            try:
                mrays = [np.concatenate([mrays[ti],rays[ti]]) for ti in range(10)]
                mweights = np.concatenate([mweights,weights])
            except:
                mrays = rays
                mweights = weights


    if rrays is True:
        return mrays,mweights,previousrho

    #Go to focus
    try:
        surf.focusI(rays,weights=weights)
    except:
        pdb.set_trace()

    return anal.hpd(mrays,weights=mweights)/1e4*180/np.pi*60**2,\
           anal.rmsCentroid(mrays,weights=mweights)/1e4*180/np.pi*60**2,\
           np.sum(mweights)

def plotNodeSpacing(rsec,smax,pmin,fun,nodegap,L=200.,pbeam=True):
    """
    Plot mirror positions and beam extent for each mirror node
    """
    plt.figure('Nodes')
    plt.clf()
    for i in range(len(smax)):
        for r in rsec[i]:
            #Mirror parameters
            z = np.sqrt(1e4**2-r**2)
            psi = np.polyval(fun[i],r)
            psi = np.max([.01,psi])
            #Plot mirrors
            rp1 = wsPrimrad(pmin[i]+L,r,z,psi)
            rp2 = wsPrimrad(pmin[i],r,z,psi)
            rs1 = wsSecrad(smax[i],r,z,psi)
            rs2 = wsSecrad(smax[i]-L,r,z,psi)
            plt.plot([rp1,rp2],[pmin[i]+L,pmin[i]],'k')
            plt.plot([rs1,rs2],[smax[i],smax[i]-L],'k')

            if np.where(r==rsec[i])[0] == 136:
                pdb.set_trace()

            if pbeam is True:
                #Plot beam
                #Set up ray
                ray = sources.pointsource(0.,1)
                ray[6] = -ray[6]
                ray[1][0] = rp1
                tran.transform(ray,0,0,-1e4,0,0,0)

                #Trace to shell
                surf.wsPrimary(ray,r,z,psi)
                tran.reflect(ray)
                surf.wsSecondary(ray,r,z,psi)
                tran.reflect(ray)
                rb1 = ray[1][0]
                z1 = ray[3][0]

                #Set up ray
                ray = sources.pointsource(0.,1)
                ray[6] = -ray[6]
                ray[1][0] = rp2
                tran.transform(ray,0,0,-1e4,0,0,0)

                #Trace to shell
                surf.wsPrimary(ray,r,z,psi)
                tran.reflect(ray)
                surf.wsSecondary(ray,r,z,psi)
                tran.reflect(ray)
                rb2 = ray[1][0]
                z2 = ray[3][0]

                plt.plot([rp1,rp1,rb1,0],[1e4+500.,pmin[i]+L,z1,0],'b--')
                plt.plot([rp2,rp2,rb2,0],[1e4+500.,pmin[i],z2,0],'b--')
    return None

def plotZetaRx(rsec,rx,fmt=''):
    """
    Plot zeta as a function of r for the nodes defined
    by the prescription scripts.
    """
    plt.figure('ZetaRx')
    for i in range(len(rsec)):
        plt.plot(rsec[i],np.polyval(rx[2][i],rsec[i]),fmt)
    return None
    
def traceNode(sec,nodeindex,rsec,smax,pmin,fun,nodegap,L=200.):
    """
    Trace individual node from prescription.
    """
    r = rsec[sec][nodeindex]
    z = np.sqrt(1e4**2-r**2)
    psi = np.max([np.polyval(fun[sec],r),.01])

    traceZeta(pmin[sec],R0=r,Z0=z,psi=psi,offaxis=0.,L=200.,az=100.,\
              pause=True)
    

def fullAnalysis():
    angle = np.linspace(0.,3e-3,11)
    energy = np.array([.5,1.,2.,4.,6.,10.])*1000.
    rx = [defineRx(N=ni,nodegap=50.) for ni in range(3,6)]
    rnodes = [defineNodePositions(*r) for r in rx]
    res = [[[traceXRS(rsec,*r,energy=e,rough=.3,offaxis=ang) \
          for rsec,r in zip(rnodes,rx)] for e in energy] \
           for ang in angle]
    return rx,rnodes,res

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
