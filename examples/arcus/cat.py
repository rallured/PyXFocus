import numpy as np
from numpy import sin,cos,exp,sqrt,pi,tan
import matplotlib.pyplot as plt
from scipy import interpolate
import pyfits

import traces.analyses as anal
import traces.surfaces as surf
import traces.transformations as tran
import traces.sources as sources

#CCD QE function
ccd = np.genfromtxt('/home/rallured/Dropbox/Arcus/Raytrace/160513_CCDInfo.csv',\
                    delimiter=',')
ccd = np.transpose(ccd)
detector = np.prod(ccd[1:],axis=0)
ccdQE = interpolate.interp1d(1239.8/ccd[0],detector,kind='linear',\
                             bounds_error=False,fill_value=0.)
#Define SPO reflectivity function
ref = pyfits.getdata('/home/rallured/Dropbox/Arcus/Raytrace/160512_B4C_Ir.fits')
theta = np.linspace(.2,1.3,401)
energy = np.linspace(200,2000,401)
sporef = interpolate.interp2d(theta,energy,ref,fill_value=0.)


def alignTolerances(num,azwidth=60.,axwidth=60.,f=None,catalign=np.zeros(6),\
                    returnrays=False):
    """
    Set up perfect beam and insert CAT grating
    Mess with alignment and measure spot shift at focus
    """
    #Set up converging beam
    rays = sources.convergingbeam2(12e3,-azwidth/2.,azwidth/2.,\
                                   -axwidth/2.,axwidth/2.,num,0.)
    tran.transform(rays,0,0,0,np.pi,0,0)
    tran.transform(rays,0,0,12e3,0,0,0)

    #Place CAT grating
    tran.transform(rays,*catalign)
    surf.flat(rays)
    tran.grat(rays,200.,8,1.)
    tran.itransform(rays,*catalign)

    #Go to focus
    if f is not None:
        try:
            tran.transform(rays,0,0,f,0,0,0)
            surf.flat(rays)
        except:
            pdb.set_trace()

    cx,cy = anal.centroid(rays)

    if returnrays is True:
        return rays
    
    return cx,cy

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
                                      ang=aa)
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
                                      ang=aa)
            tweights = tweights*tref
            
            #Rotate to appropriate angle
##            tran.transform(trays,0,0,0,0,0,aa)
            
            #Attempt to concatenate, if fail then set rays,ref to trays,tref
            try:
                rays = [np.concatenate([rays[ti],trays[ti]]) for ti in range(10)]
                weights = np.concatenate([weights,tweights])
                if wave=='uniform':
                    fwave = np.concatenate([fwave,twave])
            except:
                rays = trays
                weights = tweights
                if wave=='uniform':
                    fwave = twave

    #Get to plane of outermost grating
##    tran.transform(rays,0,0,focVec[-1]-(L.max()+gap+95.),0,0,0,coords=coords)
##    surf.flat(rays)

    if vis is True:
        pyfits.writeto('PrimaryPos.fits',np.array([xp,yp,zp]),clobber=True)
        pyfits.writeto('SecondaryPos.fits',np.array([xs,ys,zs]),clobber=True)

    if wave=='uniform':
        return rays,weights,focVec[-1],lmax,fwave,coords

    return rays,weights,focVec[-1],lmax,coords

def traceSPO(R,L,focVec,N,M,spanv,wave,d=.605,t=.775,offX=0.,offY=0.,\
             vis=None,ang=None,coords=None):
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
        inc = anal.grazeAngle(rays)#np.arcsin(l*ux+m*uy+n*uz)
        if np.size(wave)==1:
            refl = sporef(inc*180/np.pi,1239.8/wave)
        else:
            refl = np.diag(sporef(inc*180/np.pi,1239.8/wave[i*N:(i+1)*N]))
        
        #Trace to secondary
        surf.spoSecondary(rays,R[i],focVec[i])            
        tran.reflect(rays)
                           
        #Compute reflectivity
        inc = anal.grazeAngle(rays)#inc = np.arcsin(l*ux+m*uy+n*uz)
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

def determineTorus():
    """
    Trace the rays to the Rowland torus.
    Go to best fit plane of rays of central module,
    then rotate to blaze angle. Diffract rays, and determine
    position of diffracted focus. This will be 2.4 nm, 5th order,
    but make that variable.
    This will define three points on the Rowland circle,
    so the two torus radii can be determined.
    """
    
