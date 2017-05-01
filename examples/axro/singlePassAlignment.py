import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.optimize as opt

import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import traces.conicsolve as conic
import traces.axro.slf as slf

import utilities.imaging.analysis as ana
import utilities.imaging.man as man
import utilities.imaging.fitting as fit

pcoeff,pax,paz = np.genfromtxt('/home/rallured/Dropbox/AXRO/'
                        'Alignment/CoarseAlignment/150615_OP1S09Coeffs.txt')
pcoeff = pcoeff/1000.
primc = [pcoeff,pax,pax]

def createWavefront(rad,num,coeff,rorder=None,aorder=None,\
                    slitwidth=3.,masknum=15,trans=np.zeros(2)):
    """Bounce rays off of Zernike surface. Use flat to
    bring rays to a common plane, leaving the OPD as twice
    the figure error of the Zernike surface.
    Use subannulus so as not to waste rays in beginning
    of simulation
    Group rays for a given mask slit together using a
    Hartmann vector. Vignette everything else.
    Assume masknum slits 3 mm wide distributed evenly
    over the mirror aperture
    """
    #Create set of rays
    r1 = conic.primrad(8500.,220.,8400.)
    #Loop through Hartmann mask
    maskcenters = np.linspace(-48.5/220.,48.5/220.,masknum)
    for i in range(masknum):
        trays = sources.subannulus(220.,r1,slitwidth/220.,round(num/masknum))
        tran.transform(trays,0,0,0,0,0,maskcenters[i])
        try:
            rays = [np.concatenate([rays[ti],trays[ti]]) for ti in range(10)]
            mask = np.concatenate([mask,np.repeat(i,round(num/masknum))])
        except:
            rays = trays
            mask = np.repeat(i,round(num/masknum))

    tran.transform(rays,220.3+trans[0],trans[1],0,0,0,0)
    #Reflect to Zernike surface
    surf.zernsurf(rays,coeff,rad,nr=1.,rorder=rorder,aorder=aorder)
    tran.reflect(rays)
    surf.flat(rays,nr=1.)
    tran.transform(rays,-220.3,0,0,0,0,0)
    #Wavefront now has the proper Zernike form, rays pointing in
    #-z direction
    return rays,mask

def traceThroughPrimary(rays,mask,primalign=np.zeros(6),\
                        detalign=np.zeros(6),primCoeffs=None,cenSig=0.):
    """
    Trace rays through the primary mirror and then down to a focus.
    Need to simulate an initial misalignment and then applying
    an optimization algorithm to align primary to beam.
    Merit function should include the random error in spot centroiding
    primCoeffs is a list of coefficients, axial orders, and azimuthal orders
    Use global coordinate systems to determine sign conventions
    """
    #Move to primary reference frame - rays 200 mm above node
    tran.transform(rays,0,0,-200.,0,0,0)
    glo = [tran.tr.identity_matrix()]*4
    #Move to mirror tangent point and apply misalignment
    tran.transform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,0,*primalign[3:],coords=glo)
    tran.itransform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    #Trace to Wolter surface
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    rays = tran.applyT(rays,glo,inverse=True)
    #Rays are now at primary in global coordinate system
    #(origin on optical axis and at nominal node height)
    #Now reflect and trace down to the detector
    tran.reflect(rays)
    tran.transform(rays,0,0,-conic.primfocus(220.,8400.),0,0,0)
    #Apply detector misalignment
    tran.transform(rays,*detalign)
    surf.flat(rays)
    #Pick out spot centroids
    cen = [anal.centroid(rays,weights=mask==i) for i in range(mask[-1]+1)]
    cen = np.transpose(np.array(cen))
    #Add centroiding error
    if cenSig > 0:
        cen = cen + np.random.normal(scale=cenSig,size=np.shape(cen))

    return cen
    
def primaryTrace(rad,num,coeff,primalign=np.zeros(6),detalign=np.zeros(6),\
                 pcoeffs=None):
    """
    Function to create source rays and trace through the primary
    to the focus detector.
    """
    #Trace
    rays,mask = createWavefront(rad,num,coeff,trans=primalign[:2])
    cen = traceThroughPrimary(rays,mask,primalign=primalign,detalign=detalign,\
                              primCoeffs=pcoeffs)
    #Return deviation from centroid of centroid
    cx = np.mean(cen[0])
    cy = np.mean(cen[1])
    return cen#np.sqrt(np.mean((cx-cen[0])**2+(cy-cen[1])**2))

def alignPrimary(rad,num,coeff,initial=np.zeros(6),detalign=np.zeros(6)):
    """
    Simulate the process of aligning the primary mirror to the
    plane wave. Assume a starting misalignment of the primary as
    input to the optimizer, and then allow the angular degrees
    of freedom to vary.
    Most of this can be done using scipy.optimize.minimize
    Question of which algorithm to use, likely Nelder-Mead
    """
    fun = lambda p: primaryTrace(rad,num,coeff,\
                                 primalign=[initial[0],initial[1],initial[2],\
                                            p[0],p[1],initial[-1]],\
                                 detalign=detalign)
    
    res = opt.minimize(fun,initial[3:5],method='Nelder-Mead')

    return res

def traceThroughPair(rays,mask,primalign=np.zeros(6),\
                        detalign=np.zeros(6),cenSig=0.,\
                     primCoeffs=None,secCoeffs=None,\
                     secalign=np.zeros(6)):
    """
    Trace rays through the primary mirror and then down to a focus.
    Need to simulate an initial misalignment and then applying
    an optimization algorithm to align primary to beam.
    Merit function should include the random error in spot centroiding
    primCoeffs is a list of coefficients, axial orders, and azimuthal orders
    Use global coordinate systems to determine sign conventions
    """
    #Move to primary reference frame - rays 200 mm above node
    tran.transform(rays,0,0,-200.,0,0,0)
    glo = [tran.tr.identity_matrix()]*4
    #Move to mirror tangent point and apply misalignment
    tran.transform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,0,*primalign[3:],coords=glo)
    tran.itransform(rays,conic.primrad(8450.,220.,8400.),0,50,0,0,0,coords=glo)
    tran.transform(rays,0,0,-8400.,0,0,0,coords=glo)
    #Trace to Wolter surface
    if primCoeffs is None:
        surf.wolterprimary(rays,220.,8400.)
    else:
        surf.primaryLL(rays,220.,8400.,8500.,8400.,100./220.,\
                       *primCoeffs)
    tran.reflect(rays)
    #Place secondary in primary reference frame
    tran.transform(rays,conic.secrad(8350.,220.,8400.),0,8350.,0,0,0,\
                   coords=glo)
    tran.transform(rays,*secalign,coords=glo)
    tran.itransform(rays,conic.secrad(8350.,220.,8400.),0,8350.,0,0,0,\
                    coords=glo)
    if secCoeffs is None:
        surf.woltersecondary(rays,220.,8400.)
    else:
        surf.secondaryLL(rays,220.,8400.,8400.,8300.,100./220.,\
                         *secCoeffs)
    rays = tran.applyT(rays,glo,inverse=True)
    #Rays are now at secondary in global coordinate system
    #(origin on optical axis and at nominal node height)
    #Now reflect and go to detector
    tran.reflect(rays)
    tran.transform(rays,0,0,-8400.,0,0,0)
    tran.transform(rays,*detalign)
    surf.flat(rays)

    #Pick out spot centroids
    cen = [anal.centroid(rays,weights=mask==i) for i in range(mask[-1]+1)]
    cen = np.transpose(np.array(cen))
    #Add centroiding error
    if cenSig > 0:
        cen = cen + np.random.normal(scale=cenSig,size=np.shape(cen))

    return cen

def pairTrace(rad,num,coeff,primalign=np.zeros(6),secalign=np.zeros(6),\
              detalign=np.zeros(6)):
    """
    Function to create source rays and trace through the mirror pair
    to the focus detector.
    """
    #Trace
    rays,mask = createWavefront(rad,num,coeff)
    cen = traceThroughPair(rays,mask,primalign=primalign,detalign=detalign,\
                           secalign=secalign)
    #Return deviation from centroid of centroid
    cx = np.mean(cen[0])
    cy = np.mean(cen[1])
    return cen#np.sqrt(np.mean((cx-cen[0])**2+(cy-cen[1])**2))

def alignSecondary(rad,num,coeff,initial=np.zeros(6),detalign=np.zeros(6),\
                   primalign=np.zeros(6)):
    """
    Simulate the process of aligning the secondary mirror to the
    primary. Assume a starting misalignment of the secondary as
    input to the optimizer, and then allow the angular degrees
    of freedom to vary.
    Primary alignment taken as fixed, output of alignPrimary.
    Most of this can be done using scipy.optimize.minimize
    Question of which algorithm to use, likely Nelder-Mead
    """
    #Optimize Hartmann test either with minimization or
    #fitting of misalignments
    fun = lambda p: pairTrace(rad,num,coeff,\
                                 secalign=[initial[0],initial[1],initial[2],\
                                            p[0],p[1],p[2]],\
                                 detalign=detalign,\
                              primalign=primalign)
    
    res = opt.minimize(fun,initial[3:6],method='Nelder-Mead')
    final = np.concatenate((initial[:3],res['x']))
    perf = slf.mirrorPair(1000,secalign=final)

    return final,perf

def zernSensitivity():
    """
    Loop through Zernike terms and determine sensitivity and
    maximum peak-to-valley allowable in the term
    """
    x,y = np.meshgrid(np.linspace(-1,1,100),np.linspace(-1,1,100))
    sens = np.zeros(64)
    for i in range(3,64):
        coeff = np.zeros(64)
        coeff[i] = .0001
        res = alignSecondary(125./2,1000,coeff,\
                             initial=[0,0,0,.3e-3,.3e-3,.3e-3])
        z = surf.zernikemod.zernsurf(x,y,0,0,1,coeff)
        sens[i] = res[1]/anal.ptov(z)

    return sens

def fullShellPair(ang=2*np.pi,secalign=[0,0,0,0,0,0]):
    """
    Trace a full shell optic with misalignments of the secondary
    with respect to primary. Allow azimuthal extent to be a variable.
    """
    #Set up ray bundle
    rays = sources.subannulus(220.,220.6,ang,100000,zhat=-1.)
    tran.transform(rays,0,0,-8400.,0,0,0)
    theta = np.arctan2(rays[2],rays[1])
    
    #Trace through optics
    surf.wolterprimary(rays,220.,8400.)
    tran.reflect(rays)
    tran.transform(rays,0,0,8400.,0,0,0)
    tran.transform(rays,*secalign)
    tran.transform(rays,0,0,-8400.,0,0,0)
    surf.woltersecondary(rays,220.,8400.)
    tran.reflect(rays)
    tran.transform(rays,0,0,8400.,0,0,0)
    tran.itransform(rays,*secalign)

    #Go to focus
    tran.transform(rays,0,0,-8400.,0,0,0)
    surf.flat(rays)

    #Plot misalignment curves
##    plt.plot(theta,rays[1],'.')
##    plt.plot(theta,rays[2],'.')
    plt.plot(rays[1],rays[2],'.')

    return rays,theta

def sourceAlignment(dx,dy,dz):
    """
    Set up a trace of rays from the fiber source to
    the OAP. Determine wavefront error of collimated
    beam.
    """
    #Source
    rays = sources.circularbeam(125./4,10000)
    tran.pointTo(rays,dx,dy,-775./2+dz,reverse=1)
    rays[0] = np.sqrt((rays[1]+dx)**2+(rays[2]+dy)**2+(775./2+dz)**2)
    pdb.set_trace()

    #Go to focus
    tran.transform(rays,0,0,0,np.pi/2,0,0)
    tran.transform(rays,0,-775./2,-775./2,0,0,0)

    #Trace to parabola
    surf.conic(rays,775.,-1.,nr=1.)
    tran.reflect(rays)

    #Restrict to 5" diameter
##    ind = np.logical_and(rays[3]>387.5-62.5,rays[3]<387.5+62.5)
##    tran.vignette(rays,ind=ind)

    #Reflect
##    for i in range(7,10):
##        rays[i] = -rays[i]
    tran.transform(rays,0,0,775./2,0,0,0)
    surf.flat(rays,nr=1.)

    pdb.set_trace()


    #Get OPD
    opd,dx0,dy0 = anal.interpolateVec(rays,0,200,200)
    opd = man.remove2DLeg(opd,xo=1,yo=0)
    opd = man.remove2DLeg(opd,xo=0,yo=1)
    pv = ana.ptov(opd)
    
    plt.figure('OPD')
    plt.imshow(opd)
    plt.colorbar()

    wavesl = np.gradient(opd,dx0)
    resy = fit.legendre2d(wavesl[0],xo=2,yo=2)
    resx = fit.legendre2d(wavesl[1],xo=2,yo=2)
    resy[0][np.isnan(opd)] = np.nan
    resx[0][np.isnan(opd)] = np.nan
    
    plt.figure('Y')
    plt.imshow(resy[0]*180/np.pi*60**2)
    plt.colorbar()

    plt.figure('X')
    plt.imshow(resx[0]*180/np.pi*60**2)
    plt.colorbar()
    
    

    return pv*1e6
