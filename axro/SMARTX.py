import traces.PyTrace as PT
import numpy as np
import traces.conicsolve as con
import pdb,sys
from traces.axro.WSverify import traceChaseParam
import matplotlib.pyplot as plt
import time
import utilities.plotting as uplt

#Need to determine resolution and effective area
#as a function of shell radius
#Loop through all ~260 shells and determine vignetting
#and resolution functions vs. pointing error
#Will also need reflectivity vs. theta

#IMD Ir reflectivity at 1 keV
irang,irref = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO'
                            '/WSTracing/'
                            '150504IrReflectivity.txt',comments=';'))
#CXC Ir optical constants
ener,delta,beta,alpha,gamma = np.transpose(np.genfromtxt('/home/rallured'
                        '/Dropbox/AXRO/WSTracing/chandraConstants.txt')[2:])
#Paul's thermal filter transmission


def reflectivityIr(ang):
    """Return reflectivity with 0.5 nm RMS roughness
    calculated with IMD (1 keV)"""
    return irref[np.argmin(np.abs(irang-ang))]

def CXCreflIr(ang,energy,rough):
    """Return reflectivity with 0.5 nm RMS roughness
    calculated using Fresnel equations, Chandra
    optical constants, and "Strehl" factor (Nevot-Croce)
    Energy supplied in eV
    Roughness in RMS nm
    """
    #Get proper optical constants
    if np.size(energy) == 1:
        ind = np.argmin(abs(ener-energy/1000.))
        b = beta[ind]*.95
        d = delta[ind]*.95
    else:
        b,d = np.zeros(len(energy)),np.zeros(len(energy))
        for i in range(len(energy)):
            ind = np.argmin(abs(ener-energy[i]/1000.))
            b[i] = beta[ind]*.95
            d[i] = delta[ind]*.95
    n = 1 - d + 1j*b
    n2 = abs(n)**2
    #Compute reflectivity in each polarization plane
    #Return mean value
    Rp = abs(n*n*np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(n*n*np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    Rs = abs(np.sin(ang)-np.sqrt(n*n-np.cos(ang)**2))**2/\
         abs(np.sin(ang)+np.sqrt(n*n-np.cos(ang)**2))**2
    R = np.mean([Rp,Rs],axis=0)
    wave = 1240./energy #wavelength in nm
    k = 2*np.pi/wave
    strehl = np.exp(-4*k**2*np.sin(ang)**2*rough**2)

    return R*strehl

def traceWSShell(num,theta,r0,z0,phigh,plow,shigh,slow,\
                 energy,rough,chaseFocus=False,bestFocus=False):
    """Trace a WS mirror pair with 10 m focal length and
    mirror axial cutoffs defined by phigh,plow
    """
    #Define annulus of source rays
    a,p,d,e = con.woltparam(r0,z0)
    r1 = PT.wsPrimRad(plow,1.,r0,z0)#np.tan(a/2.)*(plow-10000.) + r0
    r2 = PT.wsPrimRad(phigh,1.,r0,z0)#np.tan(a/2.)*(phigh-10000.) + r0
    rays = PT.annulus(r1,r2,num)
    PT.transform(rays,0,0,0,np.pi,0,0)
    PT.transform(rays,0,0,z0,0,0,0)

    #Trace to primary
    PT.wsPrimary(rays,r0,z0,1.)
    #Handle vignetting
    ind = np.logical_and(rays[3]<phigh,rays[3]>plow)
    rays = PT.vignette(rays,ind=ind)
    #Vignette rays hitting backside of mirror
    dot = rays[4]*rays[7]+rays[5]*rays[8]+rays[6]*rays[9]
    ind = dot < 0.
    rays = PT.vignette(rays,ind=ind)
    #If all rays are vignetted, return
    if np.size(rays[1]) < 1:
        return 0.,0.,0.
    #Apply pointing error
    rays = [rays[0],rays[1],rays[2],rays[3],\
            rays[4]+np.sin(theta),rays[5],-np.sqrt(1-np.sin(theta)**2),\
            rays[7],rays[8],rays[9]]
##    PT.l = PT.l + np.sin(theta)
##    PT.n = -np.sqrt(1 - PT.l**2)
    #Reflect
    PT.reflect()
    #Compute mean incidence angle for reflectivity
    ang = np.abs(np.mean(np.arcsin(dot))) #radians
    refl1 = CXCreflIr(ang,energy,rough)
    #Total rays entering primary aperture
    N1 = np.size(rays[1])

    #Trace to secondary
    PT.wsSecondary(r0,z0,1.)
    #Vignette anything outside the physical range of the mirror
    ind = np.logical_and(rays[3]>slow,rays[3]<shigh)
    rays = PT.vignette(rays,ind=ind)
    #Vignette anything hitting the backside
    dot = rays[4]*rays[7]+rays[5]*rays[8]+rays[6]*rays[9]
    ind = dot < 0.
    rays = PT.vignette(rays,ind=ind)
    if np.size(rays[1]) < 1:
        return 0.,0.,0.
    PT.reflect()
    #Compute mean incidence angle for reflectivity
    ang = np.abs(np.mean(np.arcsin(dot))) #radians
    refl2 = CXCreflIr(ang,energy,rough)

    #Trace to focal plane
    rays = PT.flat(rays)

##    #Find Chase focus
##    delta = 0.
##    if chaseFocus or bestFocus:
##        cx,cy = PT.centroid()
##        r = np.sqrt(cx**2+cy**2)
##        delta = .0625*(1.+1)*(r**2*(phigh-plow)/10000.**2)\
##                *(1/np.tan(a))**2
##        PT.transform(0,0,delta,0,0,0)
##        PT.flat()
##
##    #Find best focus
##    delta2 = 0.
##    delta3 = 0.
##    if bestFocus:
##        try:
##            tran.focusI(rays,weights=
##        except:
##            pdb.set_trace()
##        PT.flat()

    return refl1*refl2,rays
    #return PT.hpd(), PT.rmsCentroid(), delta

def evaluateShell(theta,alpha):
    """Compute vignetting factor and HPD as a function of
    pointing error. Supply pointing errors theta and intersection
    radius of shell. Assumption is 200 mm long segments.
    """
    r0 = 10000.*np.tan(alpha)
    hpd = np.zeros(np.size(theta))
    rms = np.copy(hpd)
    delta = np.copy(hpd)
    a,p,d,e = con.woltparam(r0,10000.)
    for t in theta:
        hpd[t==theta],rms[t==theta],delta[t==theta] = \
            traceWSShell(10000,t,r0,10000.,11000.,\
                         10000.,10000.,8000.,1000.,.5,\
                         chaseFocus=True)

    return hpd,rms,delta

def SXperformance(theta,energy,rough,bestsurface=False,optsurface=False):
    """Go through a SMART-X prescription file and compute
    area weighted performance for a flat focal plane
    """
    #Load in rx data
##    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
##        'mirror-design-260sh-200mmlong-040mmthick'
##        '-3mdiam-10mfl-10arcmin-fov-planarIntersept032713.csv',\
##                       delimiter=','))
    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                    '150528_Pauls_Rx.csv',delimiter=','))
    geo = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                     'geometric_transmission_102711.txt'))
    therm = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/'
                                       'WSTracing/thermal_shield_transmission_102711.txt'))
    
    f = np.sqrt(rx[1][-1]**2+10000.**2)
    z = np.sqrt(f**2-rx[1]**2) #spherical
##    z = np.repeat(10000.,np.size(rx[1]))
##    ind = rx[0] > 210.
##    rx = rx[:,ind]
    Ns = np.shape(rx)[1]

    #Loop through and compute a resolution and a weight for each shell
    hpdTelescope = np.zeros(np.size(theta))
    rmsTelescope = np.zeros(np.size(theta))
    delta = np.zeros(np.size(theta))
    cent = np.zeros(np.size(theta))
    platefrac =  np.zeros(np.size(theta))
    #fig = plt.figure()
    for t in theta[:]:
        xi = np.array([])
        yi = np.array([])
        l = np.array([])
        m = np.array([])
        n = np.array([])
        weights = np.array([])
        #plt.clf()
        tstart = time.time()
        plate = np.zeros(Ns)
        for s in np.arange(0,Ns):
            if geo[1][s] > 0.:
                sys.stdout.write('Shell: %03i \r' % s)
                sys.stdout.flush()
                r,rays = traceWSShell(1000,t,rx[1][s],z[s],z[s]+225.,z[s]+25.,\
                                     z[s]-25.,z[s]-225.,energy,rough)
                r = r*geo[1][s]*rx[9][s] #Reflectivity*area*alignmentbars*vign
                #Account for thermal shield in shells 220-321
                if s > 219:
                    r = r * therm[1][np.abs(energy/1000.-therm[0]).argmin()]
                r = np.repeat(r,np.size(rays[1]))
                weights = np.append(weights,r)
                PT.conic(1107.799202,-1.)
                xi = np.append(xi,rays[1])
                yi = np.append(yi,rays[2])
                l = np.append(l,rays[4])
                m = np.append(m,rays[5])
                n = np.append(n,rays[6])
                if s%10==0:
                    plt.plot(rays[1][:100],rays[2][:100],'.')
                plate[s] = anal.centroid(rays,weights=r)[0]
        print time.time()-tstart

        #Have list of photon positions and weights
        #Need to compute centroid and then FoM
        #Normalize weights
        weights = weights/np.sum(weights)
        xi = np.array(xi,order='F')
        yi = np.array(yi,order='F')
        zi = np.zeros(np.size(xi)).astype('float')
        li = np.array(l,order='F')
        mi = np.array(m,order='F')
        ni = np.array(n,order='F')
        uxi = np.zeros(np.size(xi)).astype('float')
        uyi = np.zeros(np.size(xi)).astype('float')
        uzi = np.zeros(np.size(xi)).astype('float')
        rays = [np.zeros(np.size(xi)).astype('float'),\
                xi,yi,zi,\
                li,mi,ni,\
                uxi,uyi,uzi]
        if bestsurface:
            rays = tran.transform(rays,0,0,.25,0,0,0)
            surf.focusI(rays,weights=weights)
        if optsurface:
            PT.conic(1107.799202,-1.) #Emprically found best surface 1128.058314
        
        #Compute FoM
        rmsTelescope[t==theta] = PT.rmsCentroid(weights=weights)/10000.
        hpdTelescope[t==theta] = PT.hpd(weights=weights)/10000.
        cx,cy = PT.centroid()
        cent[t==theta] = cx
        ind = geo[1] > 0.
        platefrac[t==theta] = np.std(plate[ind]/1e4)/rmsTelescope[t==theta]
        
        print hpdTelescope[t==theta],rmsTelescope[t==theta]

    return hpdTelescope,rmsTelescope,delta,cent,plate

def sphericalNodes(rin,z0,fov,Nshells,N):
    """This function will iteratively scan node positions
    about a sphere around the focus. Node will start in obvious
    vignetting position. Extreme rays will be traced including
    FoV. Node will be nudged outward until vignetting no longer
    occurs. Node will then be moved by the designated mechanical
    gap. Then the next node is traced in the same fashion.
    Assumptions: 50 mm symmetric gap
    """
    #Bookkeeping parameters
    f = np.sqrt(rin**2+z0**2)
    fov = fov/60.*np.pi/180. #fov to radians
    zlist = []
    rlist = []
    for i in range(Nshells):
        #Starting radius for next shell node
        rstart = PT.wsPrimRad(z0+225.,1.,rin,z0)
        #Reduce rstart until vignetting is reached
        flag = 0
        while flag==0:
            zstart = np.sqrt(f**2-rstart**2)
            #Set up rays
            r1 = PT.wsPrimRad(zstart+25.,1.,rstart,zstart)
            r2 = PT.wsPrimRad(zstart+225.,1.,rstart,zstart)
            PT.pointsource(0.,N)
            PT.z = np.repeat(10500.,N)
            PT.x = np.linspace(r1,r2,N)
            PT.n = np.repeat(-1.,N)
            #Perform trace and add FoV deflections to rays
            PT.wsPrimary(rstart,zstart,1.)
            PT.l = np.repeat(np.sin(fov),N)
            PT.n = -np.sqrt(1 - PT.l**2)
            #Verify that rays do not hit prior primary
            PT.wsPrimary(rin,z0,1.)
            if np.sum(PT.z<z0+225.) != 0:
                #Ray has hit
                print 'Ray hits prior primary!'
                flag = 1
            #Verify that rays do not hit prior secondary
            PT.wsPrimary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rin,z0,1.)
            if np.sum(PT.z > z0-225.) != 0:
                print 'Ray hits prior secondary!'
                flag = 1
            #Look at other deflection
            PT.pointsource(0.,N)
            PT.z = np.repeat(10500.,N)
            PT.x = np.linspace(r1,r2,N)
            PT.n = np.repeat(-1.,N)
            #Perform trace and add FoV deflections to rays
            PT.wsPrimary(rstart,zstart,1.)
            PT.l = np.repeat(-np.sin(fov),N)
            PT.n = -np.sqrt(1 - PT.l**2)
            #Verify that rays do not hit prior primary
            PT.wsPrimary(rin,z0,1.)
            if np.sum(PT.z<z0+225.) != 0:
                #Ray has hit
                print 'Ray hits prior primary!'
                flag = 1
            #Verify that rays do not hit prior secondary
            PT.wsPrimary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rstart,zstart,1.)
            PT.reflect()
            PT.wsSecondary(rin,z0,1.)
            if np.sum(PT.z > z0-225.) != 0:
                print 'Ray hits prior secondary!'
                flag = 1
            if flag==0:
                rstart = rstart - .01 #Take off 10 microns
##                sys.stdout.write(str(rstart)+'\n')
##                sys.stdout.flush()
        #Vignetting has been reached, append rstart and zstart
        #to list of node positions
        rlist.append(rstart)
        zlist.append(zstart)
        rin = rstart
        z0 = zstart
    return rlist,zlist

def rxPlot():
    #Get Rx
    rx = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                    '150528_Pauls_Rx.csv',delimiter=','))
    geo = np.transpose(np.genfromtxt('/home/rallured/Dropbox/AXRO/WSTracing/'
                                     'geometric_transmission_102711.txt'))
    rx = rx[:,geo[1]>0]
    geo = geo[1][geo[1]>0]
    f = np.sqrt(rx[1][-1]**2+10000.**2)
    z = np.sqrt(f**2-rx[1]**2) #spherical

    #Make plot
    plt.figure('SX')
    plt.clf()
    for i in np.arange(0,len(geo),3):
        rp1 = con.primrad(z[i]+50.,rx[1][i],1e4)
        rp2 = con.primrad(z[i]+250.,rx[1][i],1e4)
        rh1 = con.secrad(z[i]-50.,rx[1][i],1e4)
        rh2 = con.secrad(z[i]-250.,rx[1][i],1e4)
        uplt.isoplot([rp1,rp2],[z[i]+50.-1e4,z[i]+250.-1e4],'b')
        uplt.isoplot([rh1,rh2],[z[i]-50.-1e4,z[i]-250.-1e4],'b')

    return rx,geo
