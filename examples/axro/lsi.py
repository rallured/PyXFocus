import numpy as np
import matplotlib.pyplot as plt
import zernikemod as zern
import traces.surfaces as surf
import traces.transformations as tran
import traces.analyses as anal
import traces.sources as sources
import utilities.fourier as fourier

import pdb

nSiO2 = 1.4570

def intensityTriWave(coeff,L,ang):
    """Simulate the intensity observed a distance L from
    the grating. Standard Zernike coefficients, L, and
    the diffraction angle ang are used as input.
    """
    k = 2*np.pi/405.e-6 #blue wavevector
    x,y = np.meshgrid(np.linspace(-1.1,1.1,1000),np.linspace(-1.1,1.1,1000))
    m = np.sin(ang)

    coeff = np.array(coeff).astype('float')
    coeff = np.tile(coeff,(3,1))
    coeff[0][2] = -m/2.
    coeff[1][1] = m/2.*np.sqrt(3)/2
    coeff[1][2] = m/4
    coeff[2][1] = -m/2.*np.sqrt(3)/2
    coeff[2][2] = m/4

    #Construct three phases
    phi1 = zern.zernsurf(x,y-m*L,0.,0.,1.,coeff[0])
    phi2 = zern.zernsurf(x-m*L*np.sqrt(3)/2,y+m*L/2,0,0,1,coeff[1])
    phi3 = zern.zernsurf(x+m*L*np.sqrt(3)/2,y+m*L/2,0,0,1,coeff[2])

    #Transform into complex exponentials and combine
    i = np.abs(np.exp(1j*phi1*k)+np.exp(1j*phi2*k)+np.exp(1j*phi3*k))**2

    return phi1,phi2,phi3,i

def intensityQuadWave(coeff,L,ang):
    """Simulate the intensity observed a distance L from
    the grating. Standard Zernike coefficients, L, and
    the diffraction angle ang are used as input.
    """
    k = 2*np.pi/405.e-6 #blue wavevector
    x1 = np.fft.fftfreq(1000,.5/75.)
    x,y = np.meshgrid(x1,x1)
    m = 62.5*np.sin(ang)

    coeff = np.array(coeff).astype('float')
    inp = zern.zernsurf(x,y,0,0,125./2,coeff)
    dx = np.diff(x)[0][0]
    gx,gy = np.gradient(inp,dx)
    print anal.ptov(inp)
    
    coeff = np.tile(coeff,(4,1))
    coeff[0][1] = m/2./np.sqrt(2)
    coeff[0][2] = -m/2./np.sqrt(2)
    coeff[1][1] = -m/2./np.sqrt(2)
    coeff[1][2] = m/2./np.sqrt(2)
    coeff[2][1] = m/2./np.sqrt(2)
    coeff[2][2] = m/2./np.sqrt(2)
    coeff[3][1] = -m/2./np.sqrt(2)
    coeff[3][2] = -m/2./np.sqrt(2)

    #Construct three phases
    shift = np.sin(ang)*L/np.sqrt(2)
    phi1 = zern.zernsurf(x-shift,y-shift,0.,0.,125/2.,coeff[0])
    phi2 = zern.zernsurf(x+shift,y+shift,0,0,125/2.,coeff[1])
    phi3 = zern.zernsurf(x-shift,y+shift,0,0,125/2.,coeff[2])
    phi4 = zern.zernsurf(x+shift,y-shift,0,0,125/2.,coeff[3])

    #Transform into complex exponentials and combine
    i = np.abs(np.exp(1j*phi1*k)+np.exp(1j*phi2*k)+\
               np.exp(1j*phi3*k)+np.exp(1j*phi4*k))**2

    #Turn NaNs to zeros
    i[np.isnan(i)] = 0.

    return inp,gx,gy,phi1,phi2,phi3,phi4,i



def gradientDetermination(img,L,ang):
    """Implement FFT based algorithm to determine
    phase gradients from the interference pattern."""
    #Compute spatial grids
    x1 = np.fft.fftfreq(1000,.5/75.)
    x,y = np.meshgrid(x1,x1)
    
    #Compute interferogram parameters
    fc = 2*ang/405e-6/np.sqrt(2) #carrier frequency
    dx = 150./999 #pixel size
    freq = np.fft.fftfreq(1000,dx)
    fx,fy = np.meshgrid(freq,freq)
    
    #Take the 2D FFT
    f = np.fft.fft2(img)

    #Window the harmonics
    fwy = np.copy(f)
    fwy[fy<fc/2.] = 0.
    fwy[abs(fx)>fc/2.] = 0.

    fwx = np.copy(f)
    fwx[fx<fc/2.] = 0.
    fwx[abs(fy)>fc/2.] = 0.

    #Create conversion factor
    k = 2*np.pi/405e-6
    conv = -1./L/ang/np.sqrt(2)/k

    #Convert to gradients
    fiy = np.fft.ifft2(fwy)*np.exp(-1j*2*np.pi*fc*y)
    ysl = np.fft.fftshift(np.angle(fiy))*conv
    fix = np.fft.ifft2(fwx)*np.exp(-1j*2*np.pi*fc*x)
    xsl = np.fft.fftshift(np.angle(fix))*conv

    return xsl,ysl
    


def createWavefront(rad,num,coeff,rorder=None,aorder=None):
    """Bounce rays off of Zernike surface. Use flat to
    bring rays to a common plane, leaving the OPD as twice
    the figure error of the Zernike surface.
    """
    #Create set of rays
    rays = sources.circularbeam(rad,num)
    #Reflect to Zernike surface
    surf.zernsurf(rays,coeff,rad,nr=1.,rorder=rorder,aorder=aorder)
    tran.reflect(rays)
    tran.transform(rays,0,0,0,np.pi,0,0)
    surf.flat(rays,nr=1.)
    #Wavefront now has the proper Zernike form, rays pointing in
    #+z direction
    return rays

def traceWedge(rays,t=25.,wang=1.*np.pi/180,pang=45.*np.pi/180):
    """
    Make two copies of rays and trace through a wedged plate.
    Ignore multiple reflections.
    Interpolate one OPD onto the other, take difference
    modulo wavelength
    t = plate thickness (at narrow end)
    ang = wedge angle
    """
    #Make copy
    rays2 = np.copy(rays)

    #Trace first set
    ref1 = [tran.tr.identity_matrix()]*4
    pdb.set_trace()
    tran.transform(rays,0,0,300.,pang,0,0,coords=ref1)
    surf.flat(rays,nr=1.)
    tran.reflect(rays)
    tran.transform(rays,0,0,0,np.pi/2-pang,0,0,coords=ref1)
    tran.transform(rays,0,0,-300.,0,0,0,coords=ref1)
##    tran.steerY(rays,coords=ref1)
    surf.flat(rays,nr=1.)

    #Trace second set
    ref2 = [tran.tr.identity_matrix()]*4
    pdb.set_trace()
    tran.transform(rays2,0,0,300.,pang,0,0,coords=ref2)
    surf.flat(rays2,nr=1.)
    #Refract into glass and reflect
    tran.refract(rays2,1.,nSiO2)
    tran.transform(rays2,0,0,t,0,wang,0,coords=ref2)
    surf.flat(rays2,nr=nSiO2)
    tran.reflect(rays2)
    #Refract out of glass
##    tran.itransform(rays2,0,0,t,wang,0,0,coords=ref2)
    tran.transform(rays2,0,0,0,0,-wang,0,coords=ref2)
    tran.transform(rays2,0,0,-t,0,0,0,coords=ref2)
    surf.flat(rays2,nr=nSiO2)
    tran.refract(rays2,nSiO2,1.)
    #Go to focal plane
    rays2 = tran.applyT(rays2,ref2,inverse=True)
    rays2 = tran.applyT(rays2,ref1)
    surf.flat(rays2,nr=1.)

    #Both sets of rays at same plane, should have shear and tilt
    #Interpolate OPDs onto common grid
    opd1,dx,dy = anal.interpolateVec(rays,0,200,200,\
                               xr=[rays[1].min(),rays[1].max()],\
                               yr=[rays2[2].min(),rays[2].max()])
    opd2 = anal.interpolateVec(rays2,0,200,200,\
                               xr=[rays[1].min(),rays[1].max()],\
                               yr=[rays2[2].min(),rays[2].max()])[0]

    #Convert to complex phase
    opd1 = opd1/.000635*2*np.pi % (2*np.pi)
    opd2 = opd2/.000635*2*np.pi % (2*np.pi)
    opd1 = np.exp(1j*opd1)
    opd2 = np.exp(1j*opd2)

    #Compute intensity/interferogram
    return np.abs(opd1+opd2)**2
