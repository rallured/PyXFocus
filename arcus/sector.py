import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pdb,sys,pickle,os
import astropy.io.fits as pyfits
import traces.grating as grat
import utilities.plotting as plotting
from traces.axro.SMARTX import CXCreflIr
from scipy import interpolate
import utilities.imaging.man as man
import utilities.transformations as tr
import utilities.imaging.fitting as fit
import astropy.convolution as conv

import traces.PyTrace as PT

import traces.analyses as anal
import traces.surfaces as surf
import traces.transformations as tran
import traces.sources as sources

#Load CCD QE data and define interpolation function
##ccd = np.genfromtxt('/home/rallured/Dropbox/'
##                    'Arcus/Raytrace/151002_CCDQE.csv',delimiter=',')
##ccd = np.transpose(ccd)
ccd = np.genfromtxt('/home/rallured/Dropbox/Arcus/Raytrace/160513_CCDInfo.csv',\
                    delimiter=',')
ccd = np.transpose(ccd)
detector = np.prod(ccd[1:],axis=0)
ccdQE = interpolate.interp1d(1239.8/ccd[0],detector,kind='linear',\
                             bounds_error=False,fill_value=0.)

#Nominal OPG design
bestfoc = -11778.498893609752 + 1.3342160017350473 - 5.
#Standard deviation in microns for allowable
#contribution to LSF of everything aside from
#design aberrations (i.e. SPO quality, alignments, everything else)
marg = 0.048027499983957431
outerradNom = 798.895210199+2.
yaw = 0.022509613654884453

###Load efficiency data and return interpolation function
##w = np.genfromtxt('/home/rallured/Dropbox/'
##                'Arcus/Raytrace/Wave.txt')
##o = np.genfromtxt('/home/rallured/Dropbox/'
##                  'Arcus/Raytrace/Order.txt')
##eff = np.genfromtxt('/home/rallured/Dropbox/'
##                    'Arcus/Raytrace/Eff.txt')
##def gratEff(order):
##    return interpolate.interp1d(w,eff[:,o==order],kind='linear',axis=0)

##d = np.genfromtxt('/home/rallured/Dropbox/Arcus/'
##                  'Raytrace/160613_RandyEfficiencies.csv',delimiter=',')
d = np.genfromtxt('/home/rallured/Dropbox/Arcus/'
                  'Raytrace/160721_3nmEfficiencies.csv',delimiter=',')
d = np.transpose(d)
w,eff,o = d
w = w/10.
def gratEff(order):
    return interpolate.interp1d(w[o==order],eff[o==order],kind='linear',\
                                fill_value=0.,bounds_error=False)

#Define SPO reflectivity function
ref = pyfits.getdata('/home/rallured/Dropbox/Arcus/Raytrace/160512_B4C_Ir.fits')
theta = np.linspace(.2,1.3,401)
energy = np.linspace(200,2000,401)
sporef = interpolate.interp2d(theta,energy,ref,fill_value=0.)

##How to investigate background?
##Turn wavelength into a vector. Draw from distribution Andy sent you.
##Also apply random wavevector kick. But how big? This should be limited by
##geometric stray light approximations. Then trace things through and
##return list of background photons at the focal plane. This shouldn't
##be too bad.

def meritFunction():
    """Go through the science cases, trace the relevant lines,
    sum the observation times
    """
    #Set up system
    yaw = grat.blazeYaw(1.5*pi/180,2.4,3,160.)
    bestfoc = sec.traceSector(300.,800.,12e3,100,span=375.,\
                              blazeYaw=yaw,bestFocus=None,order=-3,wave=2.4)
    marg = sec.traceSector(300.,800.,12e3,100,span=375.,blazeYaw=yaw,\
                           bestFocus=bestfoc,order=-3,wave=2.4,\
                           findMargin=True)
    
    #Case G1-1, column densities as function of radial distance
    #from cluster center
    wave = 2.16
    o = np.arange(1,8)
    o = o[(4.4<wave*o) & (wave*o<8.8)]

    return None
    


def investigateSector(Rin,Rout,F,N,wave,span=20.,d=.605,t=.775,gap=50.,\
                inc=1.5*pi/180,l=95.,\
                blazeYaw=0.,order=-1,marg=0.):
    #Get axial offset at design wavelength
    bestFoc = traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
                inc=inc,l=l,order=-3,\
                blazeYaw=blazeYaw,wave=2.4)
    
    #Get axial offset at nominal focus - this sets the
    #focal plane tilt angle
##    bestFoc0 = traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
##                inc=inc,l=l,order=-1,\
##                blazeYaw=blazeYaw,wave=4.9)[0]
    
    #Get grating efficiency function
    geff = gratEff(order)
    
    #Set up loops to investigate resolution and
    #effective area as a function of wavelength
    res = np.zeros(np.size(wave))
    area = np.copy(res)
    for i in range(np.size(wave)):
        b,r,a = 0.,0.,0.
        if geff(wave[i]) > 0.:
            r,a = traceSector(Rin,Rout,F,N,span=span,d=d,\
                            t=t,gap=gap,inc=inc,l=l,\
                            wave=wave[i],blazeYaw=blazeYaw,\
                            bestFocus=bestFoc,order=order,marg=marg)[:2]
        #Add in grating efficiency and CCD QE
##        a = a * geff(wave[i]) * ccdQE(wave[i])
        res[i] = r
        area[i] = a
        sys.stdout.write('Wave: %.2f, Area: %.2f, Res: %.2f\r'\
                         % (wave[i],area[i],res[i]))
        sys.stdout.flush()
        
    return res,area

def traceSector(Rin,Rout,F,N,span=20.,d=.605,t=.775,gap=50.,\
                inc=1.5*pi/180,l=95.,bestFocus=None,order=0,\
                blazeYaw=0.,wave=1.,marg=0.,findMargin=False,\
                analyzeLine=True,offX=0.,offY=0.):
    """Trace Arcus sector
    """
    #Trace parameters
    R = np.arange(Rin,Rout,t) #Vector of shell radii
    tg = .25*np.arctan((R+d/2)/F) #Primary graze angles
    L = d/tan(tg) #Primary length
    #L = 4*F*d/R #Vector of mirror lengths
    M = np.size(R) #Number of shells

    #Create focal length vector to implement spherical principle surface
    #This should be changed based on nominal focal lengths of
    #modules from cosine.
    focConst = F**2#+Rin**2
    focVec = sqrt(focConst-R**2)

    #Weight vector for shell radii
    weights = np.zeros(M*N)
    spanv = np.zeros(M)
    for i in range(M):
        #Full angular span of each shell
        spanv[i] = 2*np.arcsin(span/2/R[i])
        #Geometric area in each shell - 10^4 is unit conversion
        weights[i*N:(i+1)*N] = ((R[i]+d)**2-R[i]**2) * spanv[i]/2 / 10000.

    #Assign random wavelength to each ray
    if wave=='uniform':
        wave = np.random.uniform(3.6,3.6*2,size=M*N)

    #Trace rays through SPO modules
    rays,refl = traceSPO(R,L,focVec,N,M,spanv,wave,d=d,t=t,\
                         offX=offX,offY=offY)
    #Can call traceSPO repeatedly for each individual module
    #Then send rays through grating array
    
    #Refl needs to be array if wave is array
    weights = weights*refl

    #Determine outermost radius of grating array
    PT.transform(rays,0,0,focVec[-1]-(L.max()+gap+95.),0,0,0)
    PT.flat(rays)
    #outerrad should be fixed
    outerrad = np.max(sqrt(rays[1]**2+rays[2]**2))
    hubdist = sqrt(outerrad**2 + (focVec[-1]-(L.max()+gap+95.))**2)
    angle = np.arctan(outerrad/(focVec[-1]-(L.max()+gap+95.)))
    thetag = angle - 1.5*pi/180.
##    print 'Outerrad: %f\nHubdist: %f\nLmax: %f\nOuter Focus: %f\n' % \
##          (outerrad,hubdist,L.max(),focVec[-1])
    pdb.set_trace()
    #Trace grating array - need to add in grating efficiency and
    #CCD quantum efficiency
    if bestFocus is None:
        return gratArray(rays,outerrad,hubdist,angle,inc,l=l,\
                         weights=weights,order=-3,blazeYaw=blazeYaw,\
                         wave=2.4)
##        return traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
##                inc=inc,l=l,bestFocus=bestFocus,order=order,\
##                blazeYaw=blazeYaw,wave=wave,marg=marg)
        
    gratArray(rays,outerrad,hubdist,angle,inc,l=l,bestFocus=bestFocus,\
                  weights=weights,order=order,blazeYaw=blazeYaw,wave=wave)

    #Account for grating efficiency
    geff = gratEff(order)
    #Add in grating efficiency and CCD QE
    #geff and ccdQE need to be arrays if wave is array
    g = geff(wave)
    g = g.reshape(np.size(g))
    weights = weights * g * ccdQE(wave)

    #Go to focal plane
    tran.transform(rays,0,0,bestFocus,0,0,0)
    surf.flat(rays)

    if analyzeLine is False:
        return rays,weights,wave

    #Get rid of rays that made it through
##    ind = rays[1] > 0
##    rays = PT.vignette(rays,ind=ind)
##    weights = weights[ind]
    #Get rid of outliers
    ind = np.abs(rays[2]-np.average(rays[2]))<10.
    rays = PT.vignette(rays,ind=ind)
    weights = weights[ind]

    #If no rays made it through
    if np.size(rays[1])==0:
        return 0,0,0,0

    #Look at resolution of this line
    try:
        cy = np.average(rays[2],weights=weights)
    except:
        pdb.set_trace()
    if findMargin is True:
        #Compute FWHM after Gaussian convolution
        fwhms = np.linspace(1,3,201)/60.**2*pi/180*12e3
        tot = np.array([convolveLSF(rays,.001,m,weights=weights)\
                        for m in fwhms/2.35])
        #Convert to arcsec, return std of margin convolution
        marg = fwhms[np.argmin(np.abs(tot/12e3*180/pi*60**2-3.))]/2.35
        return marg,rays
    fwhm = convolveLSF(rays,.001,marg,weights=weights)
    resolution = cy/fwhm
    area = np.sum(weights)*.799*.83*.8*.9#Grat plates, azimuthal ribs, usable area
    #print resolution
    #print area
    #print sqrt((cy/3000)**2 - lsf**2)/F * 180/pi*60**2
    print 'Done'
    return resolution, area, np.nanmean(rays[1]), np.nanmean(rays[2]),\
           rays,weights

def traceArcus(N,span=20.,d=.605,t=.775,gap=50.,\
                inc=1.5*pi/180,l=95.,bestFocus=bestfoc,order=0,\
                blazeYaw=yaw,wave=1.,marg=marg,findMargin=False,\
                analyzeLine=True,offX=0.,offY=0.,calcCentroid=False,\
               vis=False):
    """Trace Arcus sector
    """

    sys.stdout.write(str(wave)+'\t'+str(order)+'\r')
    sys.stdout.flush()

    if wave is not 'uniform':
        #Compute efficiency and determine whether to proceed
        geff = gratEff(order)
        #Add in grating efficiency and CCD QE
        #geff and ccdQE need to be arrays if wave is array
        g = geff(wave)
        det = ccdQE(wave)
        if g*det==0.:
            return 0,0
        g = g.reshape(np.size(g))

    #Assign random wavelength to each ray
    if wave=='uniform':
        rays,weights,minfoc,lmax,wave,coords = defineSPOaperture(N,wave,\
                                offX=offX,offY=offY,\
                                vis=vis)
        
    else:
        rays,weights,minfoc,lmax,coords = defineSPOaperture(N,wave,offX=offX,\
                                offY=offY,\
                                vis=vis)

    #Trace rays through SPO modules

##    rays = plotting.pload('/home/rallured/Dropbox/Arcus/'
##                 'Raytrace/Performance/160516_SPORays.pkl')
##    weights = pyfits.getdata('/home/rallured/Dropbox/Arcus/'
##                             'Raytrace/Performance/160516_Weights.fits')
##    minfoc = 11972.53195
##    lmax = 115.97497
##    pdb.set_trace()

    #Determine outermost radius of grating array
    #outerrad should be fixed
    outerrad = outerradNom#np.max(sqrt(rays[1]**2+rays[2]**2))
    foc = tran.applyTPos(0,0,0,coords,inverse=True)
    hubdist = sqrt(outerrad**2 + foc[2]**2)
    angle = np.arctan(outerrad/(minfoc-(lmax+gap+95.)))
    thetag = angle - 1.5*pi/180.
    
##    print 'Outerrad: %f\nHubdist: %f\nLmax: %f\nOuter Focus: %f\n' % \
##          (outerrad,hubdist,L.max(),focVec[-1])

    #Trace grating array - need to add in grating efficiency and
    #CCD quantum efficiency
    if bestFocus is None:
        return gratArray(rays,outerrad,hubdist,angle,inc,l=l,\
                         weights=weights,order=-3,blazeYaw=blazeYaw,\
                         wave=2.4,coords=coords)
##        return traceSector(Rin,Rout,F,N,span=span,d=d,t=t,gap=gap,\
##                inc=inc,l=l,bestFocus=bestFocus,order=order,\
##                blazeYaw=blazeYaw,wave=wave,marg=marg)
        
    gratArray(rays,outerrad,hubdist,angle,inc,l=l,bestFocus=bestFocus,\
                  weights=weights,order=order,blazeYaw=blazeYaw,wave=wave,\
              offX=offX,vis=vis,coords=coords)

    #Grating vignetting
    weights = weights * (1 - abs(offX)/inc)

    if np.size(wave)==1:
        #Account for grating efficiency
        weights = weights * g * det

    #Vignette rays with no weight (evanescence)
##    rays = tran.vignette(rays,weights>0.)
##    if len(rays[1])==0:
##        return 0.,0.

    #Go to focal plane
    tran.transform(rays,0,0,bestFocus,0,0,0,coords=coords)
    surf.flat(rays)

    if vis is True:
        pyfits.writeto('FocalPlanePos.fits',\
                       np.array([rays[1],rays[2],rays[3]]),\
                       clobber=True)
        pyfits.writeto('Wavelengths.fits',\
                       wave,clobber=True)

    if analyzeLine is False:
        return rays,weights

    #Get rid of outliers
    ind = np.abs(rays[2]-np.average(rays[2]))<10.
    rays = PT.vignette(rays,ind=ind)
    weights = weights[ind]

    if calcCentroid is True:
        cent = anal.centroid(rays,weights=weights)
        if cent[0] < 100 or cent[1] < 100:
            pdb.set_trace()
        return anal.centroid(rays,weights=weights)

    #Get rid of rays that made it through
##    ind = rays[1] > 0
##    rays = PT.vignette(rays,ind=ind)
##    weights = weights[ind]

    #If no rays made it through
    if np.size(rays[1])==0:
        return 0,0,0,0

    #Look at resolution of this line
    try:
        cy = np.average(rays[2],weights=weights)
    except:
        pdb.set_trace()
    if findMargin is True:
        #Compute FWHM after Gaussian convolution
        fwhms = np.linspace(1,3,201)/60.**2*pi/180*12e3
        tot = np.array([convolveLSF(rays,.001,m,weights=weights)\
                        for m in fwhms/2.35])
        #Convert to arcsec, return std of margin convolution
        marg = fwhms[np.argmin(np.abs(tot/12e3*180/pi*60**2-3.))]/2.35
        return marg,rays
    fwhm = convolveLSF(rays,.001,marg,weights=weights)
    resolution = cy/fwhm
    weights = weights**.799*.83*.8*.9#Grat plates, azimuthal ribs, packing
    area = np.sum(weights)
    #print resolution
    #print area
    #print sqrt((cy/3000)**2 - lsf**2)/F * 180/pi*60**2
##    print 'Order: %i, Wave: %.2f\n'%(order,wave)
    return resolution, area, np.nanmean(rays[1]), np.nanmean(rays[2]),fwhm,\
           rays,weights


def defineSPOaperture(N,wave,offX=0.,offY=0.,gap=50.,vis=False):
    """
    Define a set of rays based on Ed's SPO module layout design.
    Radii come from his spreadsheet from 160503
    N is number of rays per SPO shell
    Use traceSPO for each individual module
    R0 = 320.443
    R1 = 811.607
    """
    #Go row by row
    ang1 = [-27.653,-16.592,-5.531,5.531,16.592,27.653]
    ang2 = [-23.174,-13.904,-4.635,4.635,13.904,23.174]
    ang3 = [-19.942,-11.965,-3.988,3.988,11.965,19.942]
    ang4 = [-16.065,-5.355,5.355,16.065]
    ang5 = [-14.317,-4.772,4.772,14.317]
    ang6 = [-12.972,-4.324,4.324,12.972]
    ang7 = [-11.756,-3.919,3.919,11.756]
    ang8 = [-10.791,-3.597,3.597,10.791]
    ang = [ang1,ang2,ang3,ang4,ang5,ang6,ang7,ang8]

    #Module radii
    rin = [320.443,382.638,444.833,507.027,569.222,631.417,693.612,755.807]
    rout = [376.243,438.438,500.633,562.827,625.022,681.217,749.412,811.607]

    #Module widths
    span = [50.159,49.839,49.614,89.363,82.476,77.572,86.892,82.053]

    #Visualization bookkeeping
    if vis is True:
        xp,yp,zp = [],[],[]
        xs,ys,zs = [],[],[]
    coords = [tran.tr.identity_matrix()]*4

    flag = 0
    for i in range(8):
        #Loop through module angles
        for a in ang[i]:
            #Trace parameters
            R = np.arange(rin[i],rout[i],.775) #Vector of shell radii
            tg = .25*np.arctan((R+.775/2)/12e3) #Primary graze angles
            L = .775/tan(tg) #Primary length
            if i==0:
                lmax = L.max()
            #L = 4*F*d/R #Vector of mirror lengths
            M = np.size(R) #Number of shells

            #Create focal length vector to implement spherical principle surface
            #This should be changed based on nominal focal lengths of
            #modules from cosine.
            focConst = 12e3**2#+rin[i]**2
            focVec = sqrt(focConst-R**2)

            #Weight vector for shell radii
            tweights = np.zeros(M*N)
            spanv = np.zeros(M)
            for k in range(M):
                #Full angular span of each shell
                spanv[k] = 2*np.arcsin(span[i]/2/R[k])
                #Geometric area in each shell - 10^2 is unit conversion
                tweights[k*N:(k+1)*N] = ((R[k]+.605)**2-R[k]**2)\
                                        * spanv[k]/2 / 100. / N
                #Radial vignetting factor
                betax = .605/2/L[k]
                vrad = max(0,(1-abs(offX)/betax))
                tweights[k*N:(k+1)*N] = tweights[k*N:(k+1)*N] * vrad
                #Azimuthal vignetting factor
                betay = .83/2/L[k]
                vaz = max(0,(1-abs(offY)/betay))
                tweights[k*N:(k+1)*N] = tweights[k*N:(k+1)*N] * vaz
                
            #Perform SPO module trace
            aa = a*np.pi/180
            if wave=='uniform':
                twave = np.random.uniform(4.2,4.2*2,size=M*N)
##                trays,tref = traceSPO(R,L,focVec,N,M,spanv,twave,\
##                                  offX=np.cos(aa)*offX-np.sin(aa)*offY,\
##                                  offY=np.cos(aa)*offY+np.sin(aa)*offX,\
##                                      )
                if vis is True:
                    trays,tref,prim,sec = traceSPO(R,L,focVec,N,M,spanv,twave,\
                                  offX=offX,\
                                  offY=offY,\
                                      ang=aa,vis=True)
                    xp = np.concatenate((xp,prim[0]))
                    yp = np.concatenate((yp,prim[1]))
                    zp = np.concatenate((zp,prim[2]))
                    xs = np.concatenate((xs,sec[0]))
                    ys = np.concatenate((ys,sec[1]))
                    zs = np.concatenate((zs,sec[2]))
                else:
                    trays,tref = traceSPO(R,L,focVec,N,M,spanv,twave,\
                                  offX=offX,\
                                  offY=offY,\
                                      ang=aa,scatter=True)
            else:
##                trays,tref = traceSPO(R,L,focVec,N,M,spanv,wave,\
##                                  offX=np.cos(aa)*offX-np.sin(aa)*offY,\
##                                  offY=np.cos(aa)*offY+np.sin(aa)*offX,\
##                                      )
                if vis is True:
                    trays,tref,prim,sec = traceSPO(R,L,focVec,N,M,spanv,wave,\
                                  offX=offX,\
                                  offY=offY,\
                                      ang=aa,vis=True)
                    xp = np.concatenate((xp,prim[0]))
                    yp = np.concatenate((yp,prim[1]))
                    zp = np.concatenate((zp,prim[2]))
                    xs = np.concatenate((xs,sec[0]))
                    ys = np.concatenate((ys,sec[1]))
                    zs = np.concatenate((zs,sec[2]))
                else:
                    trays,tref = traceSPO(R,L,focVec,N,M,spanv,wave,\
                                  offX=offX,\
                                  offY=offY,\
                                      ang=aa,scatter=True)
            tweights = tweights*tref
            
            #Rotate to appropriate angle
##            tran.transform(trays,0,0,0,0,0,aa)
            
            #Attempt to concatenate, if fail then set rays,ref to trays,tref
            try:
                rays = [np.concatenate([rays[ti],trays[ti]]) for ti in range(10)]
                weights = np.concatenate([weights,tweights])
                flags = np.append(flags,np.repeat(flag,len(tweights)))
                if wave=='uniform':
                    fwave = np.concatenate([fwave,twave])
            except:
                rays = trays
                weights = tweights
                flags = np.repeat(flag,len(tweights))
                if wave=='uniform':
                    fwave = twave
            flag = flag + 1

    #Get to plane of outermost grating
    tran.transform(rays,0,0,focVec[-1]-(L.max()+gap+95.),0,0,0,coords=coords)
    surf.flat(rays)

    if vis is True:
        pyfits.writeto('PrimaryPos.fits',np.array([xp,yp,zp]),clobber=True)
        pyfits.writeto('SecondaryPos.fits',np.array([xs,ys,zs]),clobber=True)

    if wave=='uniform':
        return rays,weights,focVec[-1],lmax,fwave,coords

    return rays,weights,focVec[-1],lmax,coords,flags

def traceSPO(R,L,focVec,N,M,spanv,wave,d=.605,t=.775,offX=0.,offY=0.,\
             vis=None,ang=None,coords=None,scatter=False):
    """Trace SPO surfaces sequentially. Collect rays from
    each SPO shell and set them to the PT rays at the end.
    Start at the inner radius, use the wafer and pore thicknesses
    to vignette and compute the next radius, loop while
    radius is less than Rout.
    """
    #Ray bookkeeping arrays
    trays = [np.zeros(M*N) for n in range(10)]
    if vis is True:
        xp,yp,zp = [],[],[]
        xs,ys,zs = [],[],[]

    #Loop through shell radii and collect rays
    ref = np.zeros(M*N)
    for i in range(M):
        #Set up source annulus
        rays = sources.subannulus(R[i],R[i]+d,spanv[i],N,zhat=-1.)
        z,n = rays[3],rays[6]
##        #Transform rays to be above xy plane
##        tran.transform(rays,0,0,-100.,0,0,0,coords=coords)
        #Apply angular offset
        tran.transform(rays,0,0,0,0,0,ang)
        #Get initial positions
        glob = None
        if vis is True:
            #Initial ray positions
            x0,y0,z0 = np.copy(rays[1]),np.copy(rays[2]),np.copy(rays[3])
##
##            #Establish 3d figure
##            fig = plt.figure('vis')
##            ax = fig.add_subplot(111,projection='3d')
        #Trace to primary
        surf.spoPrimary(rays,R[i],focVec[i])
        
        #Export primary ray positions in global reference frame
        if vis is True:
            tran.transform(rays,0,0,-focVec[i],0,0,0)
            xp = np.concatenate((xp,rays[1]))
            yp = np.concatenate((yp,rays[2]))
            zp = np.concatenate((zp,rays[3]))
            tran.transform(rays,0,0,focVec[i],0,0,0)
                           
        #Add offsets if they apply
        rays = [rays[0],rays[1],rays[2],rays[3],\
                rays[4]+offX,rays[5]+offY,\
                -np.sqrt(rays[6]**2-offX**2-offY**2),\
                rays[7],rays[8],rays[9]]
        tran.reflect(rays)
                           
        #Compute reflectivity
        inc = PT.grazeAngle(rays)#np.arcsin(l*ux+m*uy+n*uz)
        if np.size(wave)==1:
            refl = sporef(inc*180/np.pi,1239.8/wave)
        else:
            refl = np.diag(sporef(inc*180/np.pi,1239.8/wave[i*N:(i+1)*N]))
        
        #Trace to secondary
        surf.spoSecondary(rays,R[i],focVec[i])            
        tran.reflect(rays)

        #Add scatter
        if scatter is True:
            theta = np.arctan2(rays[2],rays[1])
            inplane = np.random.normal(scale=10./2.35*5e-6,size=N)
            outplane = np.random.normal(scale=1.5/2.35*5e-6,size=N)
            rays[4] = rays[4] + inplane*np.cos(theta) + outplane*np.sin(theta)
            rays[5] = rays[5] + inplane*np.sin(theta) + outplane*np.cos(theta)
            rays[6] = -np.sqrt(1.-rays[5]**2-rays[4]**2)
                           
        #Compute reflectivity
        inc = PT.grazeAngle(rays)#inc = np.arcsin(l*ux+m*uy+n*uz)
        if np.size(wave)==1:
            ref[i*N:(i+1)*N] = refl * sporef(inc*180/np.pi\
                                             ,1239.8/wave)
        else:
            ref[i*N:(i+1)*N] = refl * np.diag(sporef(inc*180/np.pi\
                                             ,1239.8/wave[i*N:(i+1)*N]))
        
        #Set plane to be at focus
        tran.transform(rays,0,0,-focVec[i],0,0,0,coords=coords)
        #Collect rays
        try:
            for t in range(1,7):
                temp = trays[t]
                temp[i*N:(i+1)*N] = rays[t]
        except:
            pdb.set_trace()

    #Export secondary ray positions in global coordinate frame    
    if vis is True:
        return trays,ref,[xp,yp,zp],[trays[1],trays[2],trays[3]]
        
    return trays,ref

def gratArray(rays,outerrad,hubdist,angle,inc,l=95.,bestFocus=None,\
              weights=None,order=0,blazeYaw=0.,wave=1.,offX=0.,\
              coords=None,vis=False):
    """Trace rays leaving SPO petal to the fanned grating array.
    Start with outermost radius and rotate grating array about
    the hub. Define outermost grating position by max ray radius
    at desired axial height.
    Rays have been traced to bottom of outermost grating.
    """
    #Visualization bookkeeping
    xg,yg,zg = [np.zeros(len(rays[1]))]*3
    
    x,y = rays[1:3]
    #Dummy rays to ensure return of reference frame
##    rays2 = sources.subannulus(220.,223.,10.*pi/180,100)
    #Put origin at bottom of outermost grating
    tran.transform(rays,outerrad,0,0,0,0,0,coords=coords)
##    PT.transform(rays2,outerrad,0,0,0,0,0)
    #Go to proper incidence angle of grating
    tran.transform(rays,0,0,0,0,0,-pi/2,coords=coords)
    tran.transform(rays,0,0,0,-pi/2-angle+inc,0,0,coords=coords)
##    PT.transform(rays2,0,0,0,0,0,-pi/2)
##    PT.transform(rays2,0,0,0,-pi/2-angle+inc,0,0)
    #Go to hub
    tran.transform(rays,0,0,0,0,0,blazeYaw,coords=coords) #Put in blaze
    tran.transform(rays,0,hubdist,0,0,0,0,coords=coords)
##    PT.transform(rays2,0,0,0,0,0,blazeYaw) #Put in blaze
##    PT.transform(rays2,0,hubdist,0,0,0,0)
    #Trace out gratings until no rays hit a grating
    #Flat
    #Indices
    #Reflect
    #Apply Grating
    #Next
    #Edit to only flat rays with a substantial incidence angle
    indg = np.abs(np.arcsin(rays[6])) > .001
    surf.flat(rays,ind=indg)
    rho = -sqrt(x**2+y**2)*np.sign(y)
    ind = np.logical_and(rho>hubdist,rho<l+hubdist)
##    pdb.set_trace()
    ind2 = np.copy(ind)
    #offx subtracted to prevent numerical vignetting...this
    #is accounted for with numerical factors, so don't want
    #any rays to be missed
    ang = l*sin(inc-offX)/hubdist*.95
    
    i = 0
    prev = np.copy(ind)
    #Loop condition needs to be rays not diffracted > 0
    while np.sum(prev)<len(rays[1]):
        i = i+1
        if np.sum(ind2)>0:
            tran.reflect(rays,ind=ind2)
            tran.radgrat(rays,160./hubdist,order,wave,ind=ind2)
        tran.transform(rays,0,0,0,ang,0,0,coords=coords)
##        PT.transform(rays2,0,0,0,ang,0,0)
        indg = np.abs(np.arcsin(rays[6])) > .001
        indg = np.logical_and(np.invert(prev),indg)
        surf.flat(rays,ind=indg)
##        pdb.set_trace()
        #Determine rays hitting new grating
        rho = -sqrt(x**2+y**2)*np.sign(y)
        ind = np.logical_and(rho>hubdist,rho<l+hubdist)
        ind2 = np.logical_and(np.invert(prev),ind) #Remove previous rays
        prev = np.logical_or(prev,ind) #Add rays hitting new grating
        #sys.stdout.write('%i \r' % i)
        #sys.stdout.flush()
    tran.reflect(rays,ind=ind2)
    tran.radgrat(rays,160./hubdist,order,wave,ind=ind2)

##    #Go to focal plane
##    PT.transform(rays,0,-hubdist,0,0,0,0)
##    PT.transform(rays,0,0,0,0,0,-blazeYaw) #Reverse blaze
##    #Currently at bottom point of innermost grating
##    pdb.set_trace()

    if vis is True:
        #Get hub position
        hub = tran.applyTPos(0,0,0,coords,inverse=True)
        pyfits.writeto('HubPos.fits',np.array(hub),clobber=True)
    
    #Get back to original outermost grating reference frame
    tran.transform(rays,0,0,0,-ang*i,0,0,coords=coords)
    tran.transform(rays,0,-hubdist,0,0,0,0,coords=coords)
    tran.transform(rays,0,0,0,0,0,-blazeYaw,coords=coords)
    tran.transform(rays,0,0,0,pi/2+angle-inc,0,0,coords=coords)
    tran.transform(rays,0,0,0,0,0,pi/2,coords=coords)
    tran.transform(rays,-outerrad,0,0,0,0,0,coords=coords)

##    PT.transform(rays2,0,0,0,-ang*i,0,0)
##    PT.transform(rays2,0,-hubdist,0,0,0,0)
##    PT.transform(rays2,0,0,0,0,0,-blazeYaw)
##    PT.transform(rays2,0,0,0,pi/2+angle-inc,0,0)
##    PT.transform(rays2,0,0,0,0,0,pi/2)
##    PT.transform(rays2,-outerrad,0,0,0,0,0)

    #Export grating ray positions
    if vis is True:
        rays2 = tran.applyT(rays,coords,inverse=True)
        pyfits.writeto('GratingPos.fits',\
                       np.array([rays2[1],rays2[2],rays2[3]]),\
                       clobber=True)
        pdb.set_trace()

    
    #Should be there
    surf.flat(rays)
    
##    PT.transform(rays,0,hubdist,0,0,0,0)
##    PT.transform(rays,0,0,0,-ang*i+pi/2+angle-inc,0,0)
##    PT.transform(rays,0,0,0,0,0,pi/2)
##    PT.flat(rays)
    
    #Find focus
    if bestFocus is None:
        return surf.focusY(rays,weights=weights)
##
##    #Focus already found, tracing diffracted line
##    PT.transform(rays,0,0,bestFocus,0,0,0)
##    PT.flat(rays)

    return None
    

def arc(inc,yaw,hubdist,wave,order,dpermm):
    """Return x and y positions of diffraction arc
    as a function of wavelength for a given order"""
    #Set up source ray
    PT.circularbeam(0.,1)
    #Transform to grating frame
    PT.transform(0,0,0,pi/2+inc,0,0)
    PT.transform(0,0,0,0,0,yaw)
    PT.flat()
    #Apply grating
    PT.reflect()
    PT.radgrat(hubdist,dpermm,order,wave)
    #Go to focus
    PT.transform(0,0,0,0,0,-yaw)
    PT.transform(0,hubdist,0,0,0,0)
    PT.transform(0,0,0,pi/2,0,0)
    PT.flat()
    #Get ray positions
    return PT.x[0],PT.y[0]

def plotArc(inc,yaw,hubdist,dpermm,order):
    """Plot arc vs. x and wavelength in a double plot"""
    #Define wavelength vector
    wave = np.linspace(.1,10.,1000)
    #Figure out coordinates for each order
    x = np.zeros(np.size(wave))
    y = np.copy(x)
    for i in range(np.size(wave)):
        x[i],y[i] = arc(inc,yaw,hubdist,wave[i],-1,dpermm)

    pdb.set_trace()
    #Wavelength scale
    scale = np.nanmean(wave/x/order)
    fn = lambda x: scale*x
    plotting.pltd2(x,y,fn)

    return x,y,scale

def effPlot(wave,w1,w2):
    """Make grating efficiency plot
    Indicate order ranges using low 1st order wave w1
    and high 1st order wave w2
    """
    color = ['b','r','g','c','k','y','m']
    for order in np.arange(-7,-1):
        ef = gratEff(order)
        plt.plot(wave,ef(wave),label=str(order),color=color[order+7])
        #plt.plot([-w1/order,-w1/order],[0.,.7],'--',color=color[order+7])
        #plt.plot([-w2/order,-w2/order],[0.,.7],'--',color=color[order+7])

def predictionPlots(wave,arcus,w1,w2):
    """Make plots of resolution and effective area
    Plot only regions covered by CCDs
    """
    eafig = plt.figure()
    plt.grid()
    resfig = plt.figure()
    plt.grid()
    for i in range(np.shape(arcus)[0]):
        order = i+1
        ind = np.logical_and(wave<=w2/order,wave>=w1/order)
        plt.figure(resfig.number)
        plt.plot(wave[ind],arcus[i,0][ind],label=str(order))
        plt.figure(eafig.number)
        plt.plot(wave[ind],arcus[i,1][ind]*4,label=str(order))

def plotPetal(x,y,tra=tr.identity_matrix(),color='blue'):
    """Make an isometric plot of a petal layout"""
    #Apply transformation matrices
    arc = [x,y,np.zeros(len(x)),np.ones(len(x))]
    arc = np.dot(tra,arc)
    #Plot arc
    plotting.isoplot(arc[0],arc[1],color=color)
    #Plot petal aperture
    appx = np.array([300.,300.,800.,800.,300.])
    appy = np.array([187.5,-187.5,-187.5,187.5,187.5])
    app = [appx,appy,np.zeros(5),np.ones(5)]
    app = np.dot(tra,app)
    plt.plot(app[0],app[1],color=color)
    #Plot zero order
    zo = np.dot(tra,[615.253,0,0,1])
    plt.plot(zo[0],zo[1],'.',color=color)
    return [arc,app,zo]

def matchArc(arc):
    """Take the raytraced arc, negate the dispersion direction, and
    compute the translation and rotation transformations to
    match the two spectra
    """
    #Create the second arc by negating the dispersion direction
    sh = np.shape(arc)
    arc = [arc[0],arc[1],np.zeros(sh[1]),np.ones(sh[1])]
    arc2 = np.copy(arc)
    arc2[1] = -arc2[1]

    #Transform second arc by pi/2
    r = tr.rotation_matrix(pi/2,[0,0,1])
    arc2 = np.dot(r,arc2)

    #Get gradients
    dx = np.diff(arc2[0])[0]
    gr = np.gradient(arc2[1],dx)
    pdb.set_trace()

    #Apply rotation to match slopes
    angle = pi/2+np.arctan(gr[-1])+np.arctan(gr[0])
    print angle
    r2 = tr.rotation_matrix(-angle,[0,0,1])
    arc3 = np.dot(r2,arc2)

    #Apply translation to match 
    pdb.set_trace()
##    arc3[0] = arc3[0] + (arc[0][-1]-arc3[0][0])
##    arc3[1] = arc3[1] + (arc[1][-1]-arc3[1][0])
    t = tr.translation_matrix([(arc[0][-1]-arc3[0][0]),\
                               (arc[1][-1]-arc3[1][0]),\
                               0.])
    arc3 = np.dot(t,arc3)
    tra = np.dot(t,(np.dot(r2,r)))

    #Get even closer
    dist = np.sqrt((arc[0][0]-arc3[0][-1])**2+(arc[1][0]-arc3[1][-1])**2)
    lever = np.sqrt((arc[0][0]-arc[0][-1])**2+(arc[1][0]-arc[1][-1])**2)
    angle2 = dist/lever
    new = tr.rotation_matrix(angle2,[0,0,1],point=[arc3[0][0],arc3[1][0],0])
    tra2 = np.dot(new,tra)

    #Return transformation
    return tra2,angle,angle2,new,tra

def plotTwoArcAlign(tra,ang,offX=0.,offY=0.,figure=None,label=None):
    """

    """
    arc = traceArcus(1,wave='uniform',order=-1,analyzeLine=False,\
                     offX=offX,offY=offY)[0]
    arc2 = traceArcus(1,wave='uniform',order=-1,analyzeLine=False,\
                      offX=np.cos(ang)*offX-np.sin(ang)*offY,\
                      offY=-np.cos(ang)*offY-np.sin(ang)*offX)[0]
    
    #Trace 
    arc = [arc[1],arc[2],np.zeros(len(arc[1])),np.ones(len(arc[1]))]
    arc2 = [arc2[1],-arc2[2],np.zeros(len(arc2[1])),np.ones(len(arc2[1]))]
    arc3 = np.dot(tra,arc2)

    plt.figure(figure)
    plotting.isoplot(arc[0],arc[1],'.',markersize=4,label=label)
    plt.plot(arc3[0],arc3[1],'.',markersize=4,label='Arc 2')
    #plt.text(620,300,'Arc 1',size=20,color='blue')
    #plt.text(620,280,'Arc 2',size=20,color='green')
    
    
    return

def makeLayout(arc,sepAngle):
    tra,angle = matchArc(arc) #Get nominal transformation for paired petal
    tra2 = np.dot(tr.rotation_matrix(angle/2,[0,0,1]),tra)
    tra1 = tr.rotation_matrix(angle/2,[0,0,1])

    #Create the second arc by negating the dispersion direction
    sh = np.shape(arc)
    arc = [arc[0],arc[1],np.zeros(sh[1]),np.ones(sh[1])]
    arc2 = np.copy(arc)
    arc2[1] = -arc2[1]

    #Apply transformations
    arc2 = np.dot(tra2,arc2)
    arc1 = np.dot(tra1,arc)
    
    #Find center of second arc
    xc,yc = fit.circle(arc2[0],arc2[1],0,0)[0]

    #Rotate each petal about this point
    tra1 = np.dot(tr.rotation_matrix(sepAngle,[0,0,1],point=[xc,yc,0]),tra1)
    tra2 = np.dot(tr.rotation_matrix(-sepAngle,[0,0,1],point=[xc,yc,0]),tra2)

    #Add translation to second petal
    tra2 = np.dot(tr.translation_matrix([10/np.sqrt(2),10/np.sqrt(2),0]),tra2)

    #Now plot the petals with the correct translations
    pet1 = plotPetal(arc[0],arc[1],tra=tra1,color='blue')
    pet2 = plotPetal(arc[0],-arc[1],tra=tra2,color='red')

    #Plot dashed line between corners closest to arcs
    x1,x2,y1,y2=pet1[1][0][3],pet2[1][0][2],pet1[1][1][3],pet2[1][1][2]
    plt.plot([pet1[1][0][3],pet2[1][0][2]],[pet1[1][1][3],pet2[1][1][2]],'k--')

def convolveLSF(rays,binsize,std,weights=None,plot=False):
    """Convolve a gaussian with the LSF determined by
    histogramming the ray LSF
    std should be supplied as a distance in mm"""
    #Bin up rays in dispersion direction
    n,b = np.histogram(rays[2],bins=\
                       np.arange(rays[2].mean()-.5,rays[2].mean()+.5,binsize),\
                       weights=weights)
    b = np.array([np.mean([b[i],b[i+1]]) for i in range(len(b)-1)])
    #Create convolution kernel
    gaussk = conv.Gaussian1DKernel(std/binsize)
    n2 = conv.convolve(n,gaussk)
    #Determine FWHM
    maxi = np.argmax(n2) #Index of maximum value
    #Find positive bound
    bp = b[maxi:]
    fwhmp = bp[np.argmin(np.abs(n2[maxi:]-n2.max()/2))]-b[maxi]
    bm = b[:maxi]
    fwhmm = b[maxi]-bm[np.argmin(np.abs(n2[:maxi]-n2.max()/2))]
    if plot is True:
        plt.plot(b,n)
        plt.plot(b,n2)
    return fwhmm+fwhmp

def collectFocalPlaneRays(z):
    tra2 = np.dot(tr.translation_matrix([40,-100,0]),\
               tr.rotation_matrix(pi/2,[0,0,1,0]))
    rot2 = tr.rotation_matrix(pi/2,[0,0,1,0])
    tra3 = np.dot(tr.translation_matrix([1000,-1000,0]),\
               tr.rotation_matrix(-pi/2,[0,0,1,0]))
    rot3 = tr.rotation_matrix(-pi/2,[0,0,1,0])
    tra4 = np.dot(tr.translation_matrix([1020,920,0]),\
               tr.rotation_matrix(pi,[0,0,1,0]))
    rot4 = tr.rotation_matrix(pi,[0,0,1,0])

    f = open('/home/rallured/Dropbox/Arcus/Raytrace/FocalPlaneLayout/160412_Rays.pkl','r')
    rays = pickle.load(f)
    f.close()

    rays2 = np.copy(rays)
    rays2 = [rays2[0],rays2[1],-rays2[2],rays2[3],\
             rays2[4],-rays2[5],rays2[6],\
             rays2[7],rays2[8],rays2[9]]
    rays3 = np.copy(rays)
    rays3 = [rays3[0],rays3[1],-rays3[2],rays3[3],\
             rays3[4],-rays3[5],rays3[6],\
             rays3[7],rays3[8],rays3[9]]
    rays4 = np.copy(rays)

    tran.itransform(rays2,40,-100,0,0,0,pi/2)
    tran.itransform(rays3,1000,1000,0,0,0,-pi/2)
    tran.itransform(rays4,1020,920,0,0,0,pi)


    #Plot to make sure
    plt.plot(rays[1],rays[2],'.')
    plt.plot(rays2[1],rays2[2],'.')
    plt.plot(rays3[1],rays3[2],'.')
    plt.plot(rays4[1],rays4[2],'.')

    #Transform everything up
    r = [rays,rays2,rays3,rays4]
    [tran.transform(ri,0,0,z,0,0,0) for ri in r]
    [surf.flat(ri) for ri in r]
    plt.figure()
    [plt.plot(ri[1],ri[2],'.') for ri in r]

def offAxis(wave=3.6,order=-1):
    """
    Investigate linearity of spot shifts due to angular shifts.
    Create array of spot shifts from nominal centroid
    """
    #Define off axis angles
    offx = np.linspace(-1*.3e-3,1*.3e-3,11)
    offy = np.copy(offx)

    #Nominal spot location
    nomx,nomy = traceArcus(3,wave=wave,order=order,calcCentroid=True)

    #Compute shifts due to pitch (x)
    pitchx = [nomx-traceArcus(3,wave=wave,order=order,\
                              calcCentroid=True,offX=ox)[0] for ox in offx]
    pitchy = [nomy-traceArcus(3,wave=wave,order=order,\
                              calcCentroid=True,offX=ox)[1] for ox in offx]
    #Compute shifts due to yaw (y)
    yawx = [nomx-traceArcus(3,wave=wave,order=order,\
                                calcCentroid=True,offY=oy)[0] for oy in offy]
    yawy = [nomy-traceArcus(3,wave=wave,order=order,\
                                calcCentroid=True,offY=oy)[1] for oy in offy]

    shiftx = [px+yx for px in pitchx for yx in yawx]
    shifty = [py+yy for py in pitchy for yy in yawy]

    #Loop through and compute spot shifts rigorously
    cent = [traceArcus(3,wave=wave,order=order,calcCentroid=True,\
                           offX=ox,offY=oy) for ox in offx for oy in offy]
    cent = np.array(cent)
    cent = cent.reshape((11,11,2))
    cent[:,:,0] = nomx - cent[:,:,0]
    cent[:,:,1] = nomy - cent[:,:,1]

    return [shiftx,shifty], cent

def saveScanStep(res,wave,offx,offy,order,filename):
    """Save the results of a scan step to a FITS file.
    Put the wavelength and off axis position as header keywords.
    Save the X,Y positions and the weight vector in a table.
    """
    if res[0] == 0:
        return
    
    #Create header
    head = pyfits.Header()
    head['WAVE'] = (wave,'[nm]') #(value,comment)
    head['OFFX'] = (offx,'[rad]')
    head['OFFY'] = (offy,'[rad]')
    head['ORDER'] = (order,'[integer]')
    hdu = pyfits.PrimaryHDU(header=head)

    #Divide by CCD+filter for Joern
    newweight = res[-1]/ccdQE(wave)

    #Temporary weight factor
    weights = res[1]#*.799*.83*.8*.9#Grat plates, azimuthal ribs, packing

    #Create table
    col1 = pyfits.Column(name='X',format='E',unit='mm',array=res[0][1])
    col2 = pyfits.Column(name='Y',format='E',unit='mm',array=res[0][2])
    col3 = pyfits.Column(name='Weight',format='E',unit='cm^2',array=weights)
    cols = pyfits.ColDefs([col1,col2,col3])
    tbhdu = pyfits.BinTableHDU.from_columns(cols)

    #Save to fits file
    hdulist = pyfits.HDUList([hdu,tbhdu])
    hdulist.writeto(filename,clobber=True)

    return

def makeOpticalLayoutPlot():
    """
    Make plots demonstrating optical prescription
    """
    os.chdir('/home/rallured/Dropbox/Arcus/Raytrace/FocalPlaneLayout')
    hub = pyfits.getdata('HubPos.fits')
    xf,yf,zf = pyfits.getdata('FocalPlanePos.fits')
    xg,yg,zg = pyfits.getdata('GratingPos.fits')
    xp,yp,zp = pyfits.getdata('PrimaryPos.fits')
    xs,ys,zs = pyfits.getdata('SecondaryPos.fits')
    
    #Set zero order and hub locations
    plt.figure('disp')
    hubdist = 11805.475
    ygrid = np.linspace(-3000,3000,101)
    plt.plot(yf,zf,'.',label='Diffraction Arc')
    plt.plot(hub[1],hub[2],'*',markersize=14,label='Hub')
    plt.plot(0,0,'D',markersize=10,label='SPO Focus')
    grat = lambda yg: np.sqrt(hubdist**2-(yg-hub[1])**2-(hub[0]**2))
    plt.plot(ygrid,grat(ygrid),\
             label='Grating Principle Surface')
    spo = np.sqrt(12e3**2-ygrid**2)
    plt.plot(ygrid,spo,label='SPO Principle Surface')
    plt.plot([-180,hub[1]],[grat(-180),hub[2]],'k--')
    plt.plot([180,hub[1]],[grat(180),hub[2]],'k--')
    plt.plot([-180,0],[spo[0],0],'k--')
    plt.plot([180,0],[spo[-1],0],'k--')
    plt.ylim([-1000,13e3])
##    plt.xlim([-400,600])
    plt.grid()
    plt.title('Arcus Prescription Projected in Dispersion Plane')
    plt.xlabel('Dispersion Distance (mm)')
    plt.ylabel('Axial Distance (mm)')
##    plt.legend(loc='lower left')

    plt.figure('cross')
    xgrid = np.linspace(-3000,3000,101)
    plt.plot(xf,zf,'.',label='Diffraction Arc')
    plt.plot(hub[0],hub[2],'*',markersize=14,label='Hub')
    plt.plot(0,0,'D',markersize=10,label='SPO Focus')
    grat = lambda xg: np.sqrt(hubdist**2-(xg-hub[0])**2-(hub[1]**2))
    plt.plot(xgrid,grat(xgrid),\
             label='Grating Principle Surface')
    spo = np.sqrt(12e3**2-xgrid**2)
    plt.plot(xgrid,spo,label='SPO Principle Surface')
    plt.plot([300,hub[0]],[grat(300),hub[2]],'k--')
    plt.plot([800,hub[0]],[grat(800),hub[2]],'k--')
    plt.plot([300,0],[spo[0],0],'k--')
    plt.plot([800,0],[spo[-1],0],'k--')
    plt.ylim([-1000,13e3])
##    plt.xlim([-100,900])
    plt.grid()
    plt.title('Arcus Prescription Projected in Cross-Dispersion Plane')
    plt.xlabel('Cross Dispersion Distance (mm)')
    plt.ylabel('Axial Distance (mm)')
##    plt.legend(loc='lower left')

##    wave = pyfits.getdata('Wavelengths.fits')
##    wave = wave - wave.min()
##    wave = wave/wave.max()
##    fit = plt.figure()
##    ax = fit.add_subplot(111,projection='3d')
##    cm = plt.get_cmap('rainbow')
##    
##    for i in range(len(wave)):
##        ax.plot([xg[i],xf[i]],[yg[i],yf[i]],[zg[i],zf[i]],color=cm(wave[i]))
##
    #Find arc semi-circle
    res = fit.circle(xf,yf,hub[0],hub[1])
    rad = res[1][1]
    
    x0 = 608.30180066485173
    x24 = 609.16280460037842
    y24 = 530.69291879569744
    plt.figure('Arc')
    plotting.isoplot(x0,0,'>',markersize=14)
    plt.plot(xf,yf,'.')
    plt.plot(0,0,'D',markersize=10)
    plt.plot(x24,y24,'<',markersize=14)
    plt.plot(hub[0],hub[1],'*',markersize=14)
    yd = np.linspace(hub[1]-rad,hub[1]+rad,1000)
    plt.plot(hub[0]+np.sqrt(rad**2-(yd-hub[1])**2),yd,'k--')
    plt.plot([hub[0],x24],[hub[1],y24],'k--')
    plt.plot([hub[0],x0],[hub[1],0],'k--')
    plt.plot([hub[0],hub[0]+rad],[hub[1],hub[1]],'k--')
    plt.title('Arcus Prescription Projected in Focal Plane')
    plt.xlabel('Cross Dispersion Distance (mm)')
    plt.ylabel('Dispersion Distance (mm)')
    plt.xlim([-100,800])
    plt.ylim([-200,750])
    plt.grid()

def depthoffocus(rays,weights):
    tran.transform(rays,0,0,-20,0,0,0)
    surf.flat(rays)
    lsf = [convolveLSF(rays,.001,marg,weights=weights)]
    for i in range(80):
        tran.transform(rays,0,0,.5,0,0,0)
        surf.flat(rays)
        lsf.append(convolveLSF(rays,.001,marg,weights=weights))
    return lsf

def joernScan(prefix,res=None):
    wave = []
    dw = .22
    wave.append(np.arange(.2,6.+dw,dw))
    wave.append(np.arange(3.,6.+dw,dw))
    wave.append(np.arange(1.5,5.5+dw,dw))
    wave.append(np.arange(1.,3.5+dw,dw))
    wave.append(np.arange(.5,2.5+dw,dw))
    wave.append(np.arange(.2,2+dw,dw))
    wave.append(np.arange(.2,1.5+dw,dw))
    wave.append(np.arange(.2,1.5+dw,dw))

    order = np.arange(0,-8,-1)
    offx = np.linspace(-11*.3e-3,11*.3e-3,7)
    offy = np.linspace(-8*.3e-3,8*.3e-3,7)

    
    if res is None:
        res = []
        for i in range(8):
            res.append([traceArcus(1,wave=w,order=order[i],offX=ox,offY=oy,\
                                   analyzeLine=False) for w in wave[i] \
                        for ox in offx for oy in offy])

    for i in range(8):
        wa = [w for w in wave[i] for ox in offx for oy in offy]
        oxa = [ox for w in wave[i] for ox in offx for oy in offy]
        oya = [oy for w in wave[i] for ox in offx for oy in offy]
        fn = [prefix+'_'+str(order[i])+'_%04i.fits' % n\
              for n in range(len(wa))]
        [saveScanStep(r,w,ox,oy,order[i],f) for r,w,ox,oy,f in \
         zip(res[i],wa,oxa,oya,fn)]
    
    return res

def joernScan2(prefix,res=None):
    wave = []
    dw = .22
    wave.append(np.array([1240./350.,1240/400.]))
    wave.append(np.array([1240./200.,1240/250.]))
    wave.append(np.array([1240./220.,1240./350.,1240/400.]))
    wave.append(np.array([1240./350.,1240/400.]))
    wave.append(np.array([1240./700.,1240/750.,1240./800,1240./1600]))

    order = np.arange(-1,-6,-1)
    offx = np.linspace(-11*.3e-3,11*.3e-3,7)
    offy = np.linspace(-8*.3e-3,8*.3e-3,7)

    
    if res is None:
        res = []
        for i in range(5):
            res.append([traceArcus(1,wave=w,order=order[i],offX=ox,offY=oy,\
                                   analyzeLine=False) for w in wave[i] \
                        for ox in offx for oy in offy])

    for i in range(5):
        wa = [w for w in wave[i] for ox in offx for oy in offy]
        oxa = [ox for w in wave[i] for ox in offx for oy in offy]
        oya = [oy for w in wave[i] for ox in offx for oy in offy]
        fn = [prefix+'_'+str(order[i])+'_%04i.fits' % n\
              for n in range(len(wa))]
        [saveScanStep(r,w,ox,oy,order[i],f) for r,w,ox,oy,f in \
         zip(res[i],wa,oxa,oya,fn)]
    
    return res
    
def twoChannelLayout():
    tr = pyfits.getdata('/home/rallured/Dropbox/Arcus/'
                        'Raytrace/FocalPlaneLayout/160817_Transformation.fits')
    arc = traceArcus(1,wave='uniform',order=-1,analyzeLine=False)[0]
    x,y = arc[1],arc[2]
    arc = np.array([x,y,np.zeros(len(x)),np.ones(len(x))])

    plotting.isoplot(arc[0],arc[1],'.')
    arc2 = np.copy(arc)
    arc2[1] = -arc2[1]
    arc2 = np.dot(tr,arc2)
    plt.plot(arc2[0],arc2[1],'.')

    #Get SPO aperture
    rays = defineSPOaperture(1,wave='uniform')[0]
    r = np.array([rays[1],rays[2],np.zeros(3494),np.ones(3494)])
    plt.plot(rays[1],rays[2],'.')
    r2 = np.dot(tr,r)
    plt.plot(r2[0],r2[1],'.')

    #Plot optical axes
    plt.plot(0,0,'*')
    oa2 = np.dot(tr,[0,0,0,1])
    plt.plot(oa2[0],oa2[1],'*')

def turnoffChannels(rays,weights,flags,off):
    """
    Create vector to turn off certain channels. Figure out
    new resolution and area.
    """
    #Compute original RMS and sum of weights
    area0 = np.sum(weights)
    rms0 = np.sqrt(np.average(rays[2]**2,weights=weights))

    if off is None:
        return area0,rms0*2.35*180/np.pi*60**2/12e3
    
    #And all of the booleans
    tf = (flags!=off[0])
    for o in off[1:]:
        tf = (tf) & (flags!=o)
    rays = tran.vignette(rays,ind=tf)
    weights = weights[tf]

    #Return RMS and sum of weights
    return np.sum(weights),np.sqrt(np.average(rays[2]**2,\
                            weights=weights))/12e3*180/np.pi*60**2*2.35
    

