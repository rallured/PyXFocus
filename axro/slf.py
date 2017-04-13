import numpy as np
from numpy import pi,cos,sin,tan,exp,sqrt
import matplotlib.pyplot as plt
import traces.sources as sources
import traces.surfaces as surf
import traces.transformations as tran
import traces.conicsolve as conic
import traces.analyses as anal
import pdb
import scipy.signal as sig

psiE = 1.1844342872222389

def singleEllipse(n,misalign=np.zeros(6),srcdist=89.61e3+1.5e3,az=100.,\
                 returnRays=False,f=None,\
                 plist=[[0],[0],[0]],\
                 ax=100.,psi=psiE):
    """Alternative SLF finite source trace"""
    #Establish subannulus of rays
    r0 = conic.primrad(8426.,220.,8400.)
    r1 = conic.primrad(8426.+ax,220.,8400.)
    rays = sources.subannulus(r0,r1,az/220.,n,zhat=-1.)
    tran.pointTo(rays,0.,0.,srcdist,reverse=1.)
    #Transform to node position
    tran.transform(rays,220,0,0,0,0,0)
    #Apply misalignment
    tran.transform(rays,*misalign)
    #Place mirror
    tran.transform(rays,-220.,0,-8400.,0,0,0)
##    surf.wolterprimarynode(rays,220,8400.)
    surf.ellipsoidPrimaryLL(rays,220.,8400.,srcdist,psi,8426.+ax,8426.,\
                            az/220.,*plist)
    tran.itransform(rays,-220.,0.,-8400.,0,0,0)
    #Vignette rays not landing in active mirror area
    ind = np.logical_and(rays[3]>26.,rays[3]<(26.+ax))
##    ind = np.logical_and(np.abs(rays[2])<az/2.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Reverse misalignment
    tran.itransform(rays,*misalign)
    #Reflect and go to surface
    tran.reflect(rays)
    if f is None:
        f = surf.focusI(rays)
    else:
        tran.transform(rays,0,0,f,0,0,0)
        surf.flat(rays)
    #Get centroid
    cx,cy = anal.centroid(rays)

    if returnRays is True:
        return rays
    
    return anal.hpd(rays)/abs(f)*180/pi*60**2#,f,cx

def ellipsoidPair(N,srcdist=89.61e3+1.5e3,primalign=np.zeros(6),\
               secalign=np.zeros(6),rrays=False,f=None,\
                  plist=[[0],[0],[0]],hlist=[[0],[0],[0]]):
    """
    Trace an ellipsoid-hyperboloid telescope in SLF geometry.
    plist is [pcoeff,pax,paz]
    """
    #Establish subannulus of rays
    r1 = conic.ellipsoidRad(srcdist,1.,220.,8400.,8500.)
    rays = sources.subannulus(220.,r1,100./220.,N,zhat=-1.)
    tran.pointTo(rays,0,0,srcdist,reverse=1.)
##    #Transform to node position
##    tran.transform(rays,220,0,0,0,0,0)
##    #Set up finite source distance
##    raydist = sqrt(srcdist**2+rays[1]**2+rays[2]**2)
##    rays[4] = rays[1]/raydist
##    rays[5] = rays[2]/raydist
##    rays[6] = -sqrt(1.-rays[4]**2-rays[5]**2)

    #Place mirror pair
    coords = [tran.tr.identity_matrix()]*4
    prad = conic.ellipsoidRad(srcdist,1.,220.,8400.,8450.)
    tran.transform(rays,prad,0,50.,0,0,0,\
                   coords=coords)
    tran.transform(rays,*primalign,coords=coords)
    tran.transform(rays,-prad,0,-8450.,0,0,0,\
                   coords=coords)
    surf.ellipsoidPrimaryLL(rays,220.,8400.,srcdist,1.,8500.,8400.,100./220,\
                            *plist)
    #Vignette any rays outside of active area
    rays = tran.vignette(rays,ind=np.logical_and(rays[3]<8500.,\
                                                 rays[3]>8400.))
##    surf.ellipsoidPrimary(rays,220.,8400.,srcdist,1.)
    tran.reflect(rays)
    #Place secondary in primary's reference frame
    srad = conic.ehSecRad(srcdist,1.,220.,8400.,8350.)
    tran.transform(rays,srad,0,8350.,0,0,0,\
                   coords=coords)
    tran.transform(rays,*secalign,coords=coords)
    tran.itransform(rays,srad,0,8350.,0,0,0,\
                   coords=coords)
##    surf.ellipsoidSecondary(rays,220.,8400.,srcdist,1.)
    surf.ellipsoidSecondaryLL(rays,220.,8400.,srcdist,1.,8400.,8300.,100./220,\
                              *hlist)
    rays = tran.vignette(rays,ind=np.logical_and(rays[3]<8400.,\
                                                 rays[3]>8300.))
    ang = anal.grazeAngle(rays)
    tran.reflect(rays)

    #Go back to nominal node reference frame and down to focus
    rays = tran.applyT(rays,coords,inverse=True)

    if f is None:
        f = -surf.focusI(rays)
        print f
    else:
        tran.transform(rays,0,0,-f,0,0,0)
        surf.flat(rays)

    if rrays is True:
        return rays
    
    return anal.hpd(rays)/f * 180/np.pi * 60.**2

def singleOptic(N,misalign=np.zeros(6)):
    """Trace single primary mirror from SLF finite
    source distance.
    """
    #Define some Wolter parameters
    r1 = conic.primrad(8600.,220.,8400.)
    dphi = 100./220./2
    #Set up subannulus
    rays = sources.subannulus(220.,r1,dphi*1.25,N)
##    #Set direction cosines
##    srcdist = 89.61e3+(1.5e3-misalign[2])
##    raydist = sqrt(srcdist**2+\
##                   (rays[1]-misalign[0])**2+\
##                   (rays[2]-misalign[1])**2)
##    l = (rays[1]-misalign[0])/raydist
##    m = (rays[2]-misalign[1])/raydist
##    n = -sqrt(1. - l**2 - m**2)
##    rays = [rays[0],rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Go to mirror node and apply rotational misalignment
    tran.transform(rays,220.,0,0,misalign[3],misalign[4],misalign[5])
    tran.transform(rays,-220.,0,0,0,0,0)
    #Place Wolter surfaces
    tran.transform(rays,0,0,-8400.,0,0,0)
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8426.,rays[3]<8526.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Go to focus
    surf.flat(rays)
    f = surf.focusI(rays)-8400.

    return rays#anal.hpd(rays)/abs(f)*180/pi*60**2,abs(f)

def singleOptic2(n,misalign=np.zeros(6),srcdist=89.61e3+1.5e3,az=100.,\
                 returnRays=False,f=None,\
                 plist=[[0],[0],[0]],\
                 ax=100.):
    """Alternative SLF finite source trace"""
    #Establish subannulus of rays
    r0 = conic.primrad(8426.,220.,8400.)
    r1 = conic.primrad(8426.+ax,220.,8400.)
    rays = sources.subannulus(r0,r1,az/220.,n,zhat=-1.)
    #Transform to node position
    tran.transform(rays,220,0,0,0,0,0)
    #Set up finite source distance
    raydist = sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    l = rays[1]/raydist
    m = rays[2]/raydist
    n = -sqrt(1.-l**2-m**2)
    rays = [raydist,rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Align perfectly to beam
    tran.steerX(rays)
    #Apply misalignment
    tran.transform(rays,*misalign)
    #Place mirror
    tran.transform(rays,-220.,0,-8400.,0,0,0)
##    surf.wolterprimarynode(rays,220,8400.)
    surf.primaryLL(rays,220.,8400.,8426.+ax,8426.,az/220.,*plist)
    rays = tran.vignette(rays,ind=np.logical_and(rays[3]<8400.+ax,\
                                                 rays[3]>8400.))
    tran.itransform(rays,-220.,0.,-8400.,0,0,0)
    #Vignette rays not landing in active mirror area
    ind = np.logical_and(rays[3]>26.,rays[3]<(26.+ax))
##    ind = np.logical_and(np.abs(rays[2])<az/2.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Reverse misalignment
    tran.itransform(rays,*misalign)
    #Reflect and go to surface
    tran.reflect(rays)
    if f is None:
        f = surf.focusI(rays)
    else:
        tran.transform(rays,0,0,f,0,0,0)
        surf.flat(rays)
    #Get centroid
    cx,cy = anal.centroid(rays)

    if returnRays is True:
        return rays
    
    return anal.hpd(rays)/abs(f)*180/pi*60**2,f,cx

def examineAzimuthalStripSize(strip=np.linspace(2.5,15.,5)):
    pitch = np.linspace(0,10*.3e-3,200)

    foc = []
    perf = []
    lat = []
    angle = []
    yawbound = []
    pbound = []
    dbound = []
    for s in strip:
        #Compute performance as function of pitch
        res = np.transpose(\
            np.array(\
                [singleOptic2(10000,misalign=[0,0,0,0,t,0],\
                              az=s) for t in pitch]))
        #Find optimal performance
        per = sig.savgol_filter(res[0],11,3)
        ind = np.argmin(per)
        angle.append(pitch[ind])
        perf.append(res[0][ind])
        foc.append(res[1][ind])
        lat.append(res[2][ind])
        #Find alignment sensitivities
        drange = np.linspace(-500.,500.,101)
        angrange = np.linspace(-.3e-3,.3e-3,101)
        yawres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,y,angle[-1],0],\
                                                  f=foc[-1],az=s) \
                                     for y in angrange]))
        plt.figure('yaw')
        plt.plot(angrange*180/pi*60**2,sig.savgol_filter(yawres[0],11,3),\
                 label=str(s))
        yawbound.append(abs(angrange[np.argmin(abs(yawres[0]-1.))]))
        pres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,0,p+angle[-1],0],\
                                                  f=foc[-1],az=s) \
                                     for p in angrange]))
        plt.figure('pitch')
        plt.plot(angrange*180/pi*60**2,sig.savgol_filter(pres[0],11,3),\
                 label=str(s))
        pbound.append(abs(angrange[np.argmin(abs(pres[0]-1))]))
        dres = np.transpose(np.array([singleOptic2(1000,misalign=\
                                                  [0,0,0,0,angle[-1],0],\
                                                  f=foc[-1]+d,az=s) \
                                     for d in drange]))
        plt.figure('d')
        plt.plot(drange,sig.savgol_filter(dres[0],11,3),label=str(s))
        dbound.append(abs(drange[np.argmin(abs(dres[0]-1))]))

    return np.array([perf,foc,lat,angle,yawbound,pbound,dbound])
        
def mirrorPair(N,srcdist=89.61e3+1.5e3,primalign=np.zeros(6),\
               secalign=np.zeros(6),rrays=False,f=None,\
               plist=[[0],[0],[0]],hlist=[[0],[0],[0]]):
    """
    SLF finite source trace
    """
    #Establish subannulus of rays
    rays = sources.subannulus(220.,221.,100./220.,N,zhat=-1.)
    #Transform to node position
    tran.transform(rays,220,0,0,0,0,0)
    #Set up finite source distance
    raydist = sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    rays[4] = rays[1]/raydist
    rays[5] = rays[2]/raydist
    rays[6] = -sqrt(1.-rays[4]**2-rays[5]**2)

    #Place mirror pair
    coords = [tran.tr.identity_matrix()]*4
    tran.transform(rays,-220+conic.primrad(8450.,220.,8400.),0,50.,0,0,0,\
                   coords=coords)
    tran.transform(rays,*primalign,coords=coords)
    tran.transform(rays,-conic.primrad(8450.,220.,8400.),0,-8450.,0,0,0,\
                   coords=coords)
##    surf.wolterprimary(rays,220.,8400.)
    surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220,*plist)
    rays = tran.vignette(rays,ind=np.logical_and(rays[3]<8500.,\
                                                 rays[3]>8400.))
    tran.reflect(rays)
    #Place secondary in primary's reference frame
    tran.transform(rays,conic.secrad(8350.,220.,8400.),0,8350.,0,0,0,\
                   coords=coords)
    tran.transform(rays,*secalign,coords=coords)
    tran.itransform(rays,conic.secrad(8350.,220.,8400.),0,8350.,0,0,0,\
                   coords=coords)
##    surf.woltersecondary(rays,220.,8400.)
    surf.secondaryLL(rays,220.,8400.,1.,8400.,8300.,100./220,*hlist)
    rays = tran.vignette(rays,ind=np.logical_and(rays[3]<8400.,\
                                                 rays[3]>8300.))
    tran.reflect(rays)

    #Go back to nominal node reference frame and down to focus
    rays = tran.applyT(rays,coords,inverse=True)

    if f is None:
        f = -surf.focusI(rays)
        print f
    else:
        tran.transform(rays,0,0,-f,0,0,0)
        surf.flat(rays)

    if rrays is True:
        return rays
    
    return anal.hpd(rays)/f * 180/np.pi * 60.**2, \
           airnp.mean(rays[1]), np.mean(rays[2])

def sourceToChamber(N,misalign=np.zeros(6)):
    """
    Trace randomly sampled rays from the TruFocus X-ray source
    to the 1.22 m diameter entrance to the test chamber.
    A-B from Jeff K.'s memo is 89.61
    Use oversized sub-apertured annulus, applying translations
    """
    #Define some Wolter parameters
    r1 = conic.primrad(8600.,220.,8400.)
    dphi = 100./220./2
    #Set up subannulus
    rays = sources.subannulus(220.,r1,dphi*1.25,N)
    #Set direction cosines
    srcdist = 89.61e3+(1.5e3-misalign[2])
    raydist = sqrt(srcdist**2+\
                   (rays[1]-misalign[0])**2+\
                   (rays[2]-misalign[1])**2)
    l = (rays[1]-misalign[0])/raydist
    m = (rays[2]-misalign[1])/raydist
    n = -sqrt(1. - l**2 - m**2)
    rays = [rays[0],rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Go to mirror node and apply rotational misalignment
    tran.transform(rays,220.,0,0,misalign[3],misalign[4],misalign[5])
    tran.transform(rays,-220.,0,0,0,0,0)
    #Place Wolter surfaces
    tran.transform(rays,0,0,-8400.,0,0,0)
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8426.,rays[3]<8526.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Place secondary
    surf.woltersecondary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8276.,rays[3]<8376.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Go back up to intersection plane
    tran.transform(rays,0,0,8400,0,0,0)
    #Reverse misalignments
    tran.itransform(rays,-220.,0,0,0,0,0)
    tran.itransform(rays,0,0,0,misalign[3],misalign[4],misalign[5])
    tran.itransform(rays,220,0,0,0,0,0)
    #Now back in nominal intersection coordinate system
    #Go to focus
    f = -9253.3858232
    tran.transform(rays,0,0,f,0,0,0)
    surf.flat(rays)
    
    return rays#anal.hpd(rays)/abs(f)*60**2*180/pi

def placeWolterPair(rays,misalign=np.zeros(6)):
    """
    Place the X-ray test mirror pair in the beam path.
    Assume rays are at XY plane with -z mean direction
    Nominal position of intersection plane is 1.5 m
    past chamber entrance with mirror optical axis
    coincident with chamber optical axis.
    Can supply misalignment about X=0,Y=0 in intersection plane.
    """
    #Go to nominal intersection plane
    tran.transform(rays,0,0,-1500.,0,0,0)
    #Apply misalignments
    tran.transform(rays,*misalign)
    #Go to focus and place primary
    tran.transform(rays,0,0,-8400,0,0,0)
    pdb.set_trace()
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    pdb.set_trace()
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8426.,rays[3]<8526.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Place secondary
    pdb.set_trace()
    surf.woltersecondary(rays,220.,8400.)
    tran.reflect(rays)
    #Vignette rays not landing in active mirror area
    indz = np.logical_and(rays[3]>8276.,rays[3]<8376.)
    ind = np.logical_and(np.abs(rays[2])<50.,indz)
    rays = tran.vignette(rays,ind=ind)
    #Go back up to nominal intersection plane
    tran.transform(rays,0,0,8400,0,0,0)
    tran.itransform(rays,*misalign)
    return rays

def nominalSensitivities():
    """
    Compute the alignment sensitivities for the parabola/hyperbola
    mirror pair in the SLF.
    """
    #Scan ranges
    ang = np.linspace(-5*.3e-3,5*.3e-3,100)
    tx = np.linspace(-.3,.3,100)

    #Mirror Pair Sensitivities
    pitch2 = [mirrorPair(1000,primalign=[0,0,0,0,a,0]) for a in ang]
    yaw2 = [mirrorPair(1000,primalign=[0,0,0,a,0,0]) for a in ang]
    plt.figure('Pair')
    plt.plot(ang*180/pi*60,pitch2,label='Pitch')
    plt.plot(ang*180/pi*60,yaw2,label='Yaw')
    plt.title('SLF Mirror Pair Alignment Sensitivities')
    plt.grid()
    plt.legend(loc='upper center')
    plt.xlabel('Angular Error (arcmin)')
    plt.ylabel('HPD (arcsec)')

    #Secondary Sensitivities
    pitch = [mirrorPair(1000,secalign=[0,0,0,0,a,0]) for a in ang/20]
    yaw = [mirrorPair(1000,secalign=[0,0,0,a,0,0]) for a in ang/20]
    roll = [mirrorPair(1000,secalign=[0,0,0,0,0,a]) for a in ang/20]
    plt.figure('SecondaryAng')
    plt.semilogy(ang/20.*180/pi*60**2,pitch,label='Pitch')
    plt.plot(ang/20.*180/pi*60**2,yaw,label='Yaw')
    plt.plot(ang/20.*180/pi*60**2,roll,label='Roll')
    plt.grid()
    plt.legend(loc='upper center')
    plt.xlabel('Angular Error (arcsec)')
    plt.ylabel('HPD (arcsec)')
    plt.xlim([-5,5])
    plt.ylim([0,3])
    plt.title('SLF Secondary Alignment Sensitivities')
    decenter = [mirrorPair(1000,secalign=[t,0,0,0,0,0]) for t in tx]
    lateral = [mirrorPair(1000,secalign=[0,t,0,0,0,0]) for t in tx]
    despace = [mirrorPair(1000,secalign=[0,0,t,0,0,0]) for t in tx]
    plt.figure('SecondaryTx')
    plt.semilogy(tx,decenter,label='Decenter')
    plt.plot(tx,despace,label='Despace')
    plt.plot(tx,lateral,label='Lateral')
    plt.grid()
    plt.legend(loc='upper center')
    plt.xlabel('Translation Error (mm)')
    plt.ylabel('HPD (arcsec)')
    plt.title('SLF Secondary Translation Sensitivities')
    
    
    #Compensating behavior...
    

    return [pitch2,yaw2,pitch,yaw,decenter,lateral,despace]
