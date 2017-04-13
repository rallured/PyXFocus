import traces.PyTrace as PT
import numpy as np
import traces.conicsolve as con
import sys,pdb

def traceChaseParam(num,psi,theta,alpha,L1,z0,bestFocus=False,chaseFocus=False):
    """Trace a WS mirror pair using the parameters in Eq. 13
    of Chase & VanSpeybroeck
    Return the RMS radius of the focus
    Can specify whether to find the best focus or not
    """
    #Define annulus of source rays
    r0 = np.tan(4*alpha)*z0
    alphap = 2*alpha/(1+1/psi)
    r1 = PT.wsPrimRad(z0+L1,1.,r0,z0)#np.tan(alphap/2.)*L1 + r0
    PT.annulus(r0,r1,num)
    PT.transform(0,0,0,np.pi,0,0)
    PT.transform(0,0,-10000.,0,0,0)

    pdb.set_trace()

    #Trace to primary
    PT.wsPrimary(r0,z0,psi)
    #Handle vignetting
    PT.vignette()
    ind = np.logical_and(PT.z<z0+L1,PT.z>z0)
    PT.vignette(ind=ind)
##    pdb.set_trace()
    #Apply pointing error
    PT.l = PT.l + np.sin(theta)
    PT.n = -np.sqrt(1 - PT.l**2)
    #Anything hitting backside goes away
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    #Reflect
    PT.reflect()

    N1 = np.size(PT.x)

    #Trace to secondary
    PT.wsSecondary(r0,z0,psi)
    #Vignette anything that did not converge
    PT.vignette()
    #Vignette anything hitting the backside
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
##    #Vignette anything outside the physical range of the mirror
##    ind = np.logical_and(PT.z>z0-L1,PT.z<z0)
##    PT.vignette(ind=ind)
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    PT.reflect()

    #Trace to focal plane
    PT.flat()
    cx,cy = PT.centroid()
    r = np.sqrt(cx**2+cy**2)

    #Find best focus
    delta = 0.
    if chaseFocus or bestFocus:
        delta = .0625*(psi+1)*(r**2*L1/z0**2)*(1/np.tan(alpha))**2
        PT.transform(0,0,delta,0,0,0)
        PT.flat()

    delta2 = 0.
    delta3 = 0.
    if bestFocus:
        try:
            delta2 = PT.findimageplane(20.,100)
            PT.transform(0,0,delta2,0,0,0)
            delta3 = PT.findimageplane(1.,100)
            PT.transform(0,0,delta3,0,0,0)
        except:
            pdb.set_trace()
        PT.flat()
    
    return PT.rmsCentroid()/z0, delta+delta2+delta3, r

def reproduceFig1(L1,z0,alpha=None,bestFocus=False,chaseFocus=False):
    """Attempts to reproduce Fig. 1 from Chase and Van Speybroeck
    """
    if alpha is None:
        alpha = np.array([.5,3/2.,5/2.,7/2.]) * np.pi/180.
    t = np.linspace(0.,30.,100.)/60.*np.pi/180.
    blur = np.zeros((np.size(alpha),100))
    delta = np.zeros((np.size(alpha),100))
    r = np.zeros((np.size(alpha),100))
    for a in alpha:
        for t0 in t:
            blur[a==alpha,t0==t],delta[a==alpha,t0==t],r[a==alpha,t0==t] = \
                        traceChaseParam(10000,1.,t0,a,L1,z0,\
                        bestFocus=bestFocus,chaseFocus=chaseFocus)
            sys.stdout.write(str(a)+'\t'+str(t0)+'\n')
            sys.stdout.flush()
    return blur,delta,r
