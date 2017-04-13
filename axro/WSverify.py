import traces.PyTrace as PT
import numpy as np
import traces.conicsolve as con
import sys,pdb
import matplotlib.pyplot as plt

def traceChaseParam(num,psi,theta,alpha,L1,z0,bestFocus=False,chaseFocus=False):
    """Trace a WS mirror pair using the parameters in Eq. 13
    of Chase & VanSpeybroeck
    Return the RMS radius of the focus
    Can specify whether to find the best focus or not
    """
    #Define annulus of source rays
    r0 = z0*np.tan(alpha)
    r1 = PT.wsPrimRad(z0+L1,1.,r0,z0)
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
    ind = dot > 0.
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

    pdb.set_trace()
    
    return PT.rmsCentroid()/z0, delta+delta2+delta3, r

def traceWSShell(num,theta,r0,z0,phigh,plow,shigh,slow,\
                 chaseFocus=False,bestFocus=False):
    """Trace a WS mirror pair with 10 m focal length and
    mirror axial cutoffs defined by phigh,plow
    """
    #Define annulus of source rays
    a,p,d,e = con.woltparam(r0,z0)
    r1 = PT.wsPrimRad(plow,1.,r0,z0)#np.tan(a/2.)*(plow-10000.) + r0
    r2 = PT.wsPrimRad(phigh,1.,r0,z0)#np.tan(a/2.)*(phigh-10000.) + r0
##    r2 = np.mean([r1,r2])
    PT.annulus(r1,r2,num)
    PT.transform(0,0,0,np.pi,0,0)
    PT.transform(0,0,z0,0,0,0)
##    pdb.set_trace()

    #Trace to primary
    PT.wsPrimary(r0,z0,1.)
    #Handle vignetting
    PT.vignette()
    ind = np.logical_and(PT.z<phigh,PT.z>plow)
    PT.vignette(ind=ind)
    #Vignette rays hitting backside of mirror
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
    #If all rays are vignetted, return
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    #Apply pointing error
    PT.l = PT.l + np.sin(theta)
    PT.n = -np.sqrt(1 - PT.l**2)
    #Reflect
    PT.reflect()
    #Compute mean incidence angle for reflectivity
##    ang = np.abs(np.mean(np.arcsin(dot))) #radians
##    refl1 = CXCreflIr(ang,energy,rough)
    #Total rays entering primary aperture
    N1 = np.size(PT.x)

    #Trace to secondary
    PT.wsSecondary(r0,z0,1.)
    #Vignette anything that did not converge
    PT.vignette()
    #Vignette anything outside the physical range of the mirror
    ind = np.logical_and(PT.z>slow,PT.z<shigh)
    PT.vignette(ind=ind)
    #Vignette anything hitting the backside
    dot = PT.l*PT.ux+PT.m*PT.uy+PT.n*PT.uz
    ind = dot < 0.
    PT.vignette(ind=ind)
    if np.size(PT.x) < 1:
        return 0.,0.,0.
    PT.reflect()
    #Compute mean incidence angle for reflectivity
##    ang = np.abs(np.mean(np.arcsin(dot))) #radians
##    refl2 = CXCreflIr(ang,energy,rough)

    #Trace to focal plane
    PT.flat()

    #Find Chase focus
    delta = 0.
    if chaseFocus or bestFocus:
        cx,cy = PT.centroid()
        r = np.sqrt(cx**2+cy**2)
        delta = .0625*(1.+1)*(r**2*(phigh-plow)/10000.**2)\
                *(1/np.tan(a))**2
        PT.transform(0,0,delta,0,0,0)
        PT.flat()

    #Find best focus
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

    #return refl1*refl2
    return PT.hpd(), PT.rmsCentroid(), delta

def evaluateShell(theta,alpha,bestFocus=False,chaseFocus=False):
    """Compute vignetting factor and HPD as a function of
    pointing error. Supply pointing errors theta and intersection
    radius of shell. Assumption is 200 mm long segments.
    """
    r0 = 10000.*np.tan(4.*alpha)
    hpd = np.zeros(np.size(theta))
    rms = np.copy(hpd)
    delta = np.copy(hpd)
    a,p,d,e = con.woltparam(r0,10000.)
    for t in theta:
        hpd[t==theta],rms[t==theta],delta[t==theta] = \
            traceWSShell(10000,t,r0,10000.,11000.,\
                         10000.,10000.,9000.,\
                         chaseFocus=chaseFocus,bestFocus=bestFocus)
##        pdb.set_trace()

    return hpd,rms,delta

def reproduceFig1():
    """Attempts to reproduce Fig. 1 from Chase and Van Speybroeck
    """
    theta = np.linspace(0.,30.,31)*np.pi/180/60.

    hpd,rms,delta = evaluateShell(theta,.5*np.pi/180,chaseFocus=True)
    pdb.set_trace()
    plt.plot(theta*180/np.pi*60,rms/1e4,'b--',label='Chase Raytrace')
    hpd,rms,delta = evaluateShell(theta,.5*np.pi/180,bestFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'b:',label='Best Raytrace')
    chase = con.wsRMS(1.,theta,.5*np.pi/180,1000.,10000.)
    plt.plot(theta*180/np.pi*60,chase,'b',label='Eq. 13')

    pdb.set_trace()

    hpd,rms,delta = evaluateShell(theta,1.5*np.pi/180,chaseFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'r--')
    hpd,rms,delta = evaluateShell(theta,1.5*np.pi/180,bestFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'r:')
    chase = con.wsRMS(1.,theta,1.5*np.pi/180,1000.,10000.)
    plt.plot(theta*180/np.pi*60,chase,'r')

    hpd,rms,delta = evaluateShell(theta,2.5*np.pi/180,chaseFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'g--')
    hpd,rms,delta = evaluateShell(theta,2.5*np.pi/180,bestFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'g:')
    chase = con.wsRMS(1.,theta,2.5*np.pi/180,1000.,10000.)
    plt.plot(theta*180/np.pi*60,chase,'g')

    hpd,rms,delta = evaluateShell(theta,3.5*np.pi/180,chaseFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'m--')
    hpd,rms,delta = evaluateShell(theta,3.5*np.pi/180,bestFocus=True)
    plt.plot(theta*180/np.pi*60,rms/1e4,'m:')
    chase = con.wsRMS(1.,theta,3.5*np.pi/180,1000.,10000.)
    plt.plot(theta*180/np.pi*60,chase,'m')

    #Plot formatting
    plt.ylim([0,9e-5])
    plt.xlabel('Field Position (arcmin)')
    plt.ylabel('RMS Radius (radians)')
    plt.title('Reproduction of CVS73 Fig 1')
    plt.grid()
    plt.text(9.,4e-5,r'$\alpha$=0.5')
    plt.text(18.5,4e-5,r'$\alpha$=1.5')
    plt.text(24.,4e-5,r'$\alpha$=2.5')
    plt.text(25.,2e-5,r'$\alpha$=3.5')
