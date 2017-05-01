import numpy as np
from numpy import cos,sin,sqrt,pi
import matplotlib.pyplot as plt
import pdb

import traces.sources as source
import traces.surfaces as surf
import traces.analyses as anal
import traces.transformations as tran
import traces.grating as grat
import traces.conicsolve as conic

def traceSPO(N,rin=700.,rout=737.,azwidth=66.,srcdist=89.61e3+1.5e3,\
             scatter=False):
    """
    Trace a set of rays through an SPO module using a
    finite source distance. Ignore vignetting, we are
    only interested in aberrations.
    Set up a subannulus, apply finite source effect,
    and then simply translate to inner SPO radius and
    trace outward.
    Let x be the radial direction, y the azimuthal
    """
    #Establish subannulus of rays
    rays = source.subannulus(rin,rout,azwidth/rin,N,zhat=-1.)
    #Transform to node position
    mx = np.mean([rin,rout])
    tran.transform(rays,mx,0,0,0,0,0)
    #Set up finite source distance
    raydist = np.sqrt(srcdist**2+rays[1]**2+rays[2]**2)
    l = rays[1]/raydist
    m = rays[2]/raydist
    n = -np.sqrt(1.-l**2-m**2)
    rays = [raydist,rays[1],rays[2],rays[3],l,m,n,rays[7],rays[8],rays[9]]
    #Align perfectly to beam
    #tran.steerX(rays)
    tran.transform(rays,0,0,0,0,-np.mean(rays[4]),0)
    
    #Move to SPO optical axis and trace through shells
    tran.transform(rays,-mx,0,0,0,0,0)
    R = np.arange(rin,rout+.605,.605)
    rad = np.sqrt(rays[1]**2+rays[2]**2)
    for r in R:
        #Collect relevant rays
        ind = np.logical_and(rad>r,rad<r+.605)
        if np.sum(ind)==0:
            continue
        #Trace them through system
        surf.spoPrimary(rays,r,12e3,ind=ind)
        tran.reflect(rays,ind=ind)
        surf.spoSecondary(rays,r,12e3,ind=ind)
        tran.reflect(rays,ind=ind)
    #Rays are now at secondary surfaces,
    #Add scatter
    if scatter is True:
        rays[4] = rays[4] + np.random.normal(scale=15./2.35*5e-6,size=N)
        rays[5] = rays[5] + np.random.normal(scale=1.5/2.35*5e-6,size=N)
        rays[6] = -np.sqrt(1.-rays[5]**2-rays[4]**2)
    return rays

def traceOPG(rays,hubdist=11832.911,yaw=0.,order=1,wave=1.,ang=2.5/11832.911,\
             gpitch=0.,gyaw=0.,groll=0.,\
             radapprox=False):
    """
    Trace the OPG module. Probably ignore vignetting again.
    Place perfect OPG surfaces at the correct angular distance
    to make this a reasonable approximation.
    Assume reference frame is in center of module with -z
    pointing toward hub - achieved with steerX/steerY and
    rotate inc before this function call
    Create vector to keep track of which grating each ray
    diffracts from. Separate LSFs can be identified using
    this vector.
    """
    #Establish starting coordinate system
    coords = [tran.tr.identity_matrix()]*4
    #Get -x pointing to hub
    #Question whether to rotate about z to swap x and y
    tran.transform(rays,0,0,0,0,0,-pi/2,coords=coords)
    tran.transform(rays,0,0,0,pi/2,0,0,coords=coords)
    #Go to hub, then rotate to extreme grating surface
    tran.transform(rays,0,0,0,0,0,yaw,coords=coords) #possible blaze
    tran.transform(rays,0,-11832.911,0,0,0,0,coords=coords)
    tran.transform(rays,0,0,0,-ang*7,0,0,coords=coords) #minus sign ambiguity
    #Loop through gratings, tracing rays
    left = np.repeat(True,len(rays[1]))
    record = np.zeros(len(rays[1]))
    for i in range(15):
        #If no rays left, we are done
        if np.sum(left) == 0:
            continue
        #Rays with small incidence angle removed
        indg = np.abs(np.arcsin(rays[6])) > .001
        ind = np.logical_and(left,indg)
        if np.sum(ind)==0:
            tran.transform(rays,0,0,0,ang,0,0,coords=coords)
            continue
        #Trace rays to surface
        tyaw = np.random.uniform(low=-gyaw,high=gyaw)
        tpitch = np.random.uniform(low=-gpitch,high=gpitch)
        troll = np.random.uniform(low=-groll,high=groll)
        tran.transform(rays,0,11832.911,0,0,0,0,ind=ind)
        tran.transform(rays,0,0,0,tpitch,troll,tyaw,ind=ind)
        surf.flat(rays,ind=ind)
        tran.itransform(rays,0,0,0,tpitch,troll,tyaw,ind=ind)
        tran.itransform(rays,0,11832.911,0,0,0,0,ind=ind)
        #Identify relevant rays
        ind = np.logical_and(rays[2]>11832.911-96./2,rays[2]<11832.911+96./2)
        ind = np.logical_and(ind,left)
        #Remove these rays from the set that remain
        left = np.logical_and(left,np.invert(ind))
        if np.sum(ind)==0:
            tran.transform(rays,0,0,0,ang,0,0,coords=coords)
            continue
        #Record which grating these rays diffracted from
        record[ind] = i+1
        #Diffract this set of rays
        tran.reflect(rays,ind=ind)
        tran.transform(rays,0,11832.911-hubdist,0,0,0,0,coords=coords)
        
        if radapprox is False:
            tran.radgrat(rays,160./hubdist,order,wave,ind=ind)
        else:
            ind3 = np.logical_and(rays[2]<11832.911+48.,\
                                 rays[2]>11832.911+48-9.282)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,160.,order,wave,ind=ind4)

            ind3 = np.logical_and(rays[2]<11832.911+48.-9.282,\
                                 rays[2]>11832.911+48-9.282-18.564)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,159.75,order,wave,ind=ind4)

            ind3 = np.logical_and(rays[2]<11832.911+48.-9.282-18.564,\
                                 rays[2]>11832.911+48-9.282-18.564*2)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,159.5,order,wave,ind=ind4)

            ind3 = np.logical_and(rays[2]<11832.911+48.-9.282-18.564*2,\
                                 rays[2]>11832.911+48-9.282-18.564*3)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,159.25,order,wave,ind=ind4)

            ind3 = np.logical_and(rays[2]<11832.911+48.-9.282-18.564*3,\
                                 rays[2]>11832.911+48-9.282-18.564*4)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,159.,order,wave,ind=ind4)

            ind3 = np.logical_and(rays[2]<11832.911+48.-9.282-18.564*4,\
                                 rays[2]>11832.911+48-9.282-18.564*4-12.462)
            ind4 = np.logical_and(ind3,ind)
            if np.sum(ind4)>0:
                tran.grat(rays,158.75,order,wave,ind=ind4)
            #pdb.set_trace()
            
        tran.transform(rays,0,hubdist-11832.911,0,0,0,0,coords=coords)
        #Rotate to next grating
        tran.transform(rays,0,0,0,ang,0,0,coords=coords)
    #Go back to original coordinate system
    rays = tran.applyT(rays,coords,inverse=True)

    return rays,record

def test(N,rin=700.,rout=737.,azwidth=66.,srcdist=89.61e3+1.5e3,\
         hubdist=11832.911,yaw=0.,wave=6.,order=1,\
         opgalign=[0,0,0,0,0,0],f=None,\
         rrays=False,glob=False,rcen=False,\
         groll=0.,gyaw=0.,gpitch=0.,\
         scatter=False,coordin=None,\
         radapprox=False):
    """
    Trace through the SPO module, then place the OPG module
    at its nominal position, allowing for misalignments about the
    center of the OPG module. The module tolerances can be
    investigated by a coordinate transformation around the
    OPG module placement.
    """
    #Trace through SPO module
    rays = traceSPO(N,rin=rin,rout=rout,azwidth=azwidth,srcdist=srcdist,\
                    scatter=scatter)

    #Find the nominal OPG module location using formalism
    #from Flanagan's SPIE paper
    #Go to focus, steer out X and Y, then go up a distance
    #defined using Flangan formula, this should leave you
    #at the center of the beam, therefore the center of the
    #OPG module
    if coordin is None:
        coords = [tran.tr.identity_matrix()]*4
        tran.transform(rays,0,0,0,0,-np.mean(rays[4]),0,coords=coords)
        #tran.steerX(rays,coords=coords)
        #tran.steerY(rays,coords=coords)
        tran.transform(rays,0,0,0,pi-np.mean(rays[5]),0,0,coords=coords)
        f0 = surf.focusI(rays,coords=coords)
        tran.transform(rays,np.mean(rays[1]),np.mean(rays[2]),0,0,0,0,\
                       coords=coords)
        tran.transform(rays,0,0,0,0,pi,0,coords=coords)
        tran.transform(rays,0,0,11832.911*np.cos(1.5*np.pi/180),0,0,0,coords=coords)
        tran.transform(rays,0,0,0,0,1.5*np.pi/180,0,coords=coords)
    else:
        rays = tran.applyT(rays,coordin)
        coords = np.copy(coordin)
    surf.flat(rays)
    #Now at center of central grating, with -z pointing toward hub
    tran.transform(rays,*opgalign,coords=coords)
    rays,record = traceOPG(rays,hubdist=hubdist,yaw=yaw,wave=wave,order=order,\
                           gyaw=gyaw,groll=groll,gpitch=gpitch,\
                           radapprox=radapprox)
    tran.itransform(rays,*opgalign,coords=coords)
    #Should be at same reference frame, with rays now diffracted
    if np.sum(record)==0:
        pdb.set_trace()
    rays = tran.vignette(rays,ind=record>0)
    record = record[record>0]

    #Trace to detector and determine LSF
    rays = tran.applyT(rays,coords,inverse=True)
    #surf.focusI(rays)
    if f is not None:
        try:
            tran.transform(rays,0,0,-f,0,0,0)
            surf.flat(rays)
        except:
            pdb.set_trace()

    if rcen is True:
        return anal.centroid(rays)

    if rrays is True:
        if glob is True:
            tran.transform(rays,0,0,f,0,0,0)
        return rays,record

    #Return LSF in arcseconds
    return anal.hpdY(rays)/12e3*180/pi*60**2

def misalign(dof,off):
    """
    Return array of zeros where arr[dof]=off
    Used to apply an offset to a specific dof
    """
    arr = np.zeros(6)
    arr[dof] = off
    return arr

def alignmentSensitivities():
    """
    Trace Al K Littrow configuration. For each DoF,
    apply +- offsets and measure centroids. Fit resulting
    centroids vs. offset and determine slope and non-linearity
    """
    #Get alignment parameters
    yaw = grat.blazeYaw(1.5*np.pi/180,1.575757,3,160.)
    rays,rec = test(1000,yaw=yaw,order=6,wave=1.575757/2.,rrays=True)
    f = -surf.focusY(rays)

    #Perform raytrace scans
    tr = np.linspace(-5.,5.,100)
    rot = np.linspace(-.6e-3,.6e-3,100)
    tx = np.array([[test(1000,yaw=yaw,wave=1240./1490,order=6,\
               opgalign=misalign(dof,off),\
               f=f,rcen=True) for off in tr] \
           for dof in range(3)])
    rx = np.array([[test(1000,yaw=yaw,wave=1240./1490,order=6,\
               opgalign=misalign(dof,off),\
               f=f,rcen=True) for off in rot] \
          for dof in range(3,6)])
    
    #Fit results to 2nd or 3rd order polynomial
    #Get derivative at off=0.
    #Compute peak-to-valley departure from this line
    #Compute peak-to-valley departure from best fit line
    #np.polyfit(tr,t[:,0]
    tx = tx.transpose(0,2,1)
    rx = rx.transpose(0,2,1)
    tfit = [[np.polyfit(tr,m,3) for m in t] for t in tx]
    rfit = [[np.polyfit(rot,m,3) for m in t] for t in rx]
    tslope = [[m[-2]/12e3*180/pi*60**2 for m in t] for t in tfit] #arcsec/mm
    rslope = [[m[-2]/12e3 for m in t] for t in rfit] #rad/rad=arcsec/arcsec
    tdiff = [[(m[1]*tr[-1]**2+m[0]*tr[-1]**3)/\
              (m[2]*tr[-1]) for m in t] for t in tfit]
    rdiff = [[(m[1]*rot[-1]**2+m[0]*rot[-1]**3)/\
              (m[2]*rot[-1]) for m in t] for t in rfit]
    pdb.set_trace()
    
    return [tslope,rslope,tdiff,rdiff]

def linearity(polyOrder=1,wavelength=np.linspace(0.,1.57*4.,100),\
              opgalign=[0,0,0,0,0,0]):
    """
    Trace a range of wavelengths, determine spot centroids,
    and then fit x centroid vs wavelength to a polynomial.
    Return polynomial coefficients and both RMS and PV
    deviation.
    """
    #Get alignment parameters
    yaw = grat.blazeYaw(1.5*np.pi/180,1.575757,3,160.)
    rays,rec = test(1000,yaw=yaw,order=6,wave=1.575757/2.,rrays=True)
    f = -surf.focusY(rays)

    #Perform linearity scan
    cen = [test(10000,yaw=yaw,order=1,rcen=True,opgalign=opgalign,wave=w,f=f) \
           for w in wavelength]
    cen = np.transpose(np.array(cen))

    #Perform polynomial fit
    fit = np.polyfit(wavelength,cen[1],polyOrder)
    recon = np.polyval(fit,wavelength)

    rms = np.sqrt(np.mean((recon-cen[1])**2))
    pv = np.max(np.abs(recon-cen[1]))

    #Plot
    plt.plot(wavelength,cen[1]-recon,label=str(opgalign[-1]*180/np.pi))
    
    return cen,(rms,pv)

def testRadApprox(num,order=1,wave=1.,radapprox=False,N=3,f=None,yaw=0.,\
                  azwidth=66.*.68,autofocus=False,returnMet=False,axwidth=2.5):
    """

    """
    #Set up converging source
    rays = source.convergingbeam2(12e3,-azwidth/2,azwidth/2,\
                                  -axwidth/2,axwidth/2,num,0.)
    tran.transform(rays,0,0,12e3,0,0,0)
    tran.transform(rays,0,0,0,88.5*np.pi/180.,0,0)
    tran.transform(rays,0,0,0,0,0,yaw)
    surf.flat(rays)

    #Place grating
    if radapprox is False:
        tran.reflect(rays)
        tran.transform(rays,0,-12e3,0,0,0,0)
        tran.radgrat(rays,160./12e3,order,wave)
        tran.transform(rays,0,12e3,0,0,0,0)
        tran.transform(rays,0,0,0,0,0,-yaw)
    else:
        tran.reflect(rays)
        gratedges = np.linspace(-50.,50.,N+1)
        for i in range(N):
            ind = np.logical_and(rays[2]>gratedges[i],\
                             rays[2]<gratedges[i+1])
            d = (12e3+np.mean(gratedges[i:i+2]))/12e3*160.
            if np.sum(ind)>0:
                tran.grat(rays,d,-order,wave,ind=ind)
        tran.transform(rays,0,0,0,0,0,-yaw)


    #Go to focal plane
    tran.transform(rays,0,0,0,-88.5*np.pi/180.,0,0)
    tran.transform(rays,0,0,0,0,0,np.pi/2)

    if f is not None:
        try:
            tran.transform(rays,0,0,-f,0,0,0)
            surf.flat(rays)
        except:
            pdb.set_trace()

    if autofocus is True:
        surf.focusY(rays)

    if returnMet is True:
        return anal.hpdY(rays)/12e3*180/np.pi*60**2
    
    return rays

def reproduceChevron(num,rin=220.,axlength=100.,azwidth=50.,F=8.4e3,\
                     hubdist=8e3,radapprox=False,order=1,wave=.83,f=None,\
                     autofocus=False,returnMet=False,yaw=0.,N=1,\
                     gratalign=np.zeros(6)):
    #Create Wolter beam
    rout = conic.primrad(F+axlength,rin,F)
    rays = source.subannulus(rin,rout,azwidth/rin,num,zhat=-1.)
    surf.wolterprimary(rays,rin,F)
    tran.reflect(rays)
    surf.woltersecondary(rays,rin,F)
    tran.reflect(rays)
    tran.transform(rays,0,0,0,0,0,np.pi/2)

    #Go to focus
    surf.focusI(rays)
    #Steer beam
    coords = [tran.tr.identity_matrix()]*4
    tran.transform(rays,0,0,0,np.mean(rays[5]),-np.mean(rays[4]),0)
    pdb.set_trace()
    #Go up to grating
    tran.transform(rays,0,0,hubdist/np.cos(1.5*np.pi/180),0,0,0)
    #Go to incidence angle
    tran.transform(rays,0,0,0,91.5*np.pi/180,0,0)
    tran.transform(rays,0,0,0,0,0,yaw)
    #Apply grating misalignment
    tran.transform(rays,*gratalign)
    surf.flat(rays)
    #Get rid of rays outside grating
    ind = np.logical_and(np.abs(rays[2])<16,np.abs(rays[1])<25/2.)
    rays = tran.vignette(rays,ind=ind)

    plt.figure('grat')
    plt.clf()
    plt.plot(rays[1],rays[2],'.')
    plt.title('Beam Footprint')

    #Place grating
    if radapprox is False:
        tran.reflect(rays)
        tran.transform(rays,0,-hubdist,0,0,0,0)
        tran.radgrat(rays,160./hubdist,order,wave)
        tran.transform(rays,0,hubdist,0,0,0,0)
    else:
        tran.reflect(rays)
        gratedges = np.linspace(-16.,16.,N+1)
        for i in range(N):
            ind = np.logical_and(rays[2]>gratedges[i],\
                             rays[2]<gratedges[i+1])
            d = (hubdist+np.mean(gratedges[i:i+2]))/hubdist*160.
            if np.sum(ind)>0:
                tran.grat(rays,d,-order,wave,ind=ind)


    #Go to focal plane
    tran.transform(rays,*gratalign)
    tran.transform(rays,0,0,0,0,0,-yaw)
    tran.transform(rays,0,0,0,-91.5*np.pi/180.,0,0)
    tran.transform(rays,0,0,0,0,0,np.pi/2)

    if f is not None:
        try:
            tran.transform(rays,0,0,-f,0,0,0)
            surf.flat(rays)
        except:
            pdb.set_trace()

    if autofocus is True:
        surf.focusY(rays)

    if returnMet is True:
        return anal.hpdY(rays)/F*180/np.pi*60**2

    plt.figure('LSF')
    plt.clf()
    plt.plot(rays[1],rays[2],'.')
    plt.title('LSF')
    
    return rays
