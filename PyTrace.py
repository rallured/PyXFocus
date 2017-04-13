import numpy as np
import matplotlib.pyplot as plt
import PytranTrace as tran
import zernikemod,pdb,time
import traces.conicsolve as con
import reconstruct
import math

##Get rid of global variables. Use ray vectors as input to all functions.
##Use a flag to indicate whether arrays are on the host or the GPU.
##This will break all currently written traces. Clean up how subarrays are
##handled by coding functions to create the temporary arrays.
##OPD is only needed for visible optics simulations.
##Get OPD working first, then get rid of the global variables.

def transform(rays,tx,ty,tz,rx,ry,rz,ind=None):
    """Coordinate transformation. translations are done first,
    then Rx,Ry,Rz
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.transform(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,-tx,-ty,-tz,-rx,-ry,-rz)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.transform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return


def itransform(rays,tx,ty,tz,rx,ry,rz):
    """Inverse of coordinate transformations. -rz,-ry,-rx then
    translations.
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    tran.itransform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)
    return

def reflect(rays,ind=None):
    """Reflect rays based on surface normal
    """
    l,m,n,ux,uy,uz = rays[4:]
    if ind is not None:
        tl,tm,tn,tux,tuy,tuz = l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind]
        tran.reflect(tl,tm,tn,tux,tuy,tuz)
        l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind] = tl,tm,tn,tux,tuy,tuz
    else:
        tran.reflect(l,m,n,ux,uy,uz)
    return

def flat(rays,ind=None,nr=None):
    """Trace rays to the XY plane
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if ind is not None:
        #Temporary array
        trays = [rays[i][ind] for i in range(10)]
        #Trace
        tran.flat(*trays[1:])
        #Copy back to original
        for i in range(1,10):
            rays[i][ind] = trays[i]
    elif nr is not None:
        tran.flatopd(x,y,z,l,m,n,ux,uy,uz,opd,nr)
    else:
        tran.flat(x,y,z,l,m,n,ux,uy,uz)
    return

def refract(rays,n1,n2):
    """Refract rays based on surface normal
    and ray direction cosines from index n1
    into index n2
    """
    l,m,n,ux,uy,uz = rays[4:]
    tran.refract(l,m,n,ux,uy,uz,n1,n2)
    return

def zernsurf(rays,coeff,rad,rorder=None,aorder=None):
    """Wrapper for Zernike surface
    Coordinates are usual arctan2(y,x)
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder is None:
        rorder,aorder = zernikemod.zmodes(np.size(coeff))
    tran.tracezern(x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad)
    rho = np.sqrt(x**2+y**2)
    ind = np.where(rho<=rad)
    vignette(ind=ind)
    return

def zernphase(rays,coeff,rad,wave,rorder=None,aorder=None):
    """Wrapper for standard Zernike phase surface. Supply
    wavelength in mm, radius in mm, coeff in mm."""
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder is None:
        rorder,aorder = zernikemod.zmodes(np.size(coeff))
    tran.zernphase(opd,x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad,wave)
    rho = np.sqrt(x**2+y**2)
    ind = np.where(rho<=rad)
    rays = vignette(rays,ind=ind)
    return

#Wrapper for Zernike surface with 2 Zernike sets and one with
#arbitrary rotation angle
#Coordinates are usual arctan2(y,x)
def zernsurfrot(rays,coeff1,coeff2,rad,rot,\
                rorder1=None,aorder1=None,rorder2=None,aorder2=None):
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder1 is None:
        rorder1,aorder1 = zernikemod.zmodes(np.size(coeff1))
    if rorder2 is None:
        rorder2,aorder2 = zernikemod.zmodes(np.size(coeff2))
    tran.tracezernrot(x,y,z,l,m,n,ux,uy,uz,coeff1,np.array(rorder1),\
                      np.array(aorder1),coeff2,np.array(rorder2),\
                      np.array(aorder2),rad,rot)
    rho = np.sqrt(x**2+y**2)
    ind = np.where(rho<=rad)
    vignette(ind=ind)

def sphere(rays,rad,nr=None):
    """Wrapper for spherical surface.
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if nr is not None:
        tran.tracesphereopd(opd,x,y,z,l,m,n,ux,uy,uz,rad,nr)
    else:
        tran.tracesphere(x,y,z,l,m,n,ux,uy,uz,rad)
    return

def conic(rays,R,K,nr=None):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if nr is not None:
        tran.conicopd(opd,x,y,l,m,n,ux,uy,uz,R,K,nr)
    else:
        tran.conic(x,y,z,l,m,n,ux,uy,uz,R,K)
    return

def cyl(rays,rad):
    """Wrapper for cylindrical surface
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.tracecyl(x,y,z,l,m,n,ux,uy,uz,rad)

def cylconic(rays,rad,k):
    """Wrapper for cylindrical conics
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.cylconic(x,y,z,l,m,n,ux,uy,uz,rad,k)

def wolterprimary(rays,r0,z0):
    """Wrapper for Wolter primary surface - no vignetting
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.wolterprimary(x,y,z,l,m,n,ux,uy,uz,r0,z0)
    return

def woltersecondary(rays,r0,z0):
    """Wrapper for Wolter secondary surface - no vignetting
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.woltersecondary(x,y,z,l,m,n,ux,uy,uz,r0,z0)
    return

def wolterprimtan(rays,r0,z0):
    """Wrapper for Wolter primary surface -
    place at surface tangent point
    +z is surface normal
    +y is toward sky
    +x is azimuthal direction
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    wolterprimary(r0,z0)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

def woltersine(rays,r0,z0,amp,freq):
    """Wrapper for Wolter primary surface with sinusoid
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.woltersine(x,y,z,l,m,n,ux,uy,uz,r0,z0,amp,freq)
    return

def woltersinetan(rays,r0,z0,amp,freq):
    """Wrapper for Wolter sinusoidal surface -
    place at surface tangent point
    +z is surface normal
    +y is toward sky
    +x is azimuthal direction
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    woltersine(r0,z0,amp,freq)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

def secondaryLL(rays,r0,z0,zmax,zmin,dphi,coeff,axial,az):
    """Wrapper for L-L secondary surface
    Placed at focus
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.woltersecll(x,y,z,l,m,n,ux,uy,uz,r0,z0,zmax,zmin,dphi,coeff,axial,az)
    vignette()
    return

def primaryLL(rays,r0,z0,zmax,zmin,dphi,coeff,axial,az):
    """Wrapper for L-L primary surface
    Placed at focus
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    tran.wolterprimll(x,y,z,l,m,n,ux,uy,uz,r0,z0,zmax,zmin,dphi,coeff,axial,az)
    vignette()
    return

def primaryLLtan(rays,r0,z0,zmax,zmin,dphi,coeff,axial,az):
    """Wrapper for Wolter primary surface -
    place at surface tangent point
    +z is surface normal
    +y is toward sky
    +x is azimuthal direction
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    #Compute Wolter parameters
    alpha,p,d,e = con.woltparam(r0,z0)
    transform(0,0,0,-np.pi/2-alpha,0,0)
    #Go to Wolter focus minus gap and half mirror length
    transform(0,con.primrad(z0+75.,r0,z0),-z0-75.,0,0,0)
    #Place Wolter surface
    transform(0,0,0,0,0,-np.pi/2)
    primaryLL(r0,z0,zmax,zmin,dphi,coeff,axial,az)
    transform(0,0,0,0,0,np.pi/2)
    #Go back to original coordinate system
    transform(0,-con.primrad(z0+75.,r0,z0),z0+75.,0,0,0)
    transform(0,0,0,np.pi/2+alpha,0,0)
    return

def wsPrimary(rays,r0,z0,psi):
    """Trace a W-S primary surface
    Fortran function computes Chase parameters for an equivalent W-I
    betas, f, g, and k computed from alpha and z0
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    a,p,d,e = con.woltparam(r0,z0)
    tran.wsprimary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
    return

def wsSecondary(rays,r0,z0,psi):
    """Trace a W-S secondary surface
    Fortran function computes Chase parameters for an equivalent W-I
    betas, f, g, and k computed from alpha and z0
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    a,p,d,e = con.woltparam(r0,z0)
    tran.wssecondary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
    return

def spoCone(rays,R0,tg,ind=None):
    """Trace rays to an SPO cone with intersection radius
    R0 and slope angle tg.
    XY plane should be at SPO intersection plane
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.spocone(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,R0,tg)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.spocone(x,y,z,l,m,n,ux,uy,uz,R0,tg)
    return

def spoPrimary(rays,R0,F,d=.605,ind=None):
    """Trace rays to an SPO primary with intersection radius
    R0 and focal length F.
    XY plane should be at SPO intersection plane
    """
    #Convert F to tg
    tg = .25*np.arctan((R0+d/2)/F)
    #Call SPO wrapper
    spoCone(rays,R0,tg,ind=ind)
    return

def spoSecondary(rays,R0,F,d=.605,ind=None):
    """Trace rays to an SPO secondary with intersection radius
    R0 and focal length F.
    XY plane should be at SPO intersection plane
    """
    #Convert F to tg
    tg = .75*np.arctan((R0+d/2)/F)
    #Call SPO wrapper
    spoCone(rays,R0,tg,ind=ind)
    return

def radgrat(rays,hubdist,dpermm,order,wave,ind=None):
    """Infinite radial grating. Assumes grating in x,y plane
    with grooves converging at hubdist in positive y direction
    dpermm is nm/mm
    wave is in nm
    """
    x,y,z,l,m,n = rays[1:7]
    if ind is not None:
        tx,ty,tl,tm,tn = x[ind],y[ind],l[ind],m[ind],n[ind]
        tran.radgrat(tx,ty,tl,tm,tn,hubdist,dpermm,order,wave)
        x[ind],y[ind],l[ind],m[ind],n[ind] = tx,ty,tl,tm,tn
    else:
        tran.radgrat(x,y,l,m,n,hubdist,dpermm,order,wave)
    return

###Trying CUDA
##@cuda.jit('void(double,double,double,double,'
##          'double[:],double[:],double[:],double[:],double[:],double[:])')
##def radgratC(hubdist,dpermm,order,wave,x,y,z,l,m,n):
##    """Test
##    """
##    #Establish threadId (including any block)
##    i = cuda.grid(1)
##
##    #Handle case of threadId larger than array sizes
##    if i >= x.shape[0]:
##        return
##
##    d = dpermm * math.sqrt((hubdist-y[i])**2 + x[i]**2)
##    yaw = math.pi/2 + math.atan(x[i]/(hubdist-y[i]))
##
##    det = l[i]**2 + m[i]**2
##
##    if det<1:
##        l[i] = l[i] + math.sin(yaw)*order*wave/d
##        m[i] = m[i] - math.cos(yaw)*order*wave/d
##        n[i] = math.sqrt(1.-l[i]**2-m[i]**2)
##    

def grat(rays,d,order,wave):
    """Linear grating with groove direction in +y
    Evanescence results in position vector set to zero
    """
    x,y,z,l,m,n = rays[1:7]
    tran.grat(x,y,l,m,n,d,order,wave)
    return

def vignette(rays,ind=None):
    """Remove vignetted rays from memory
    ind is array of "good" indices, all others are removed
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if ind==None:
        mag = l**2+m**2+n**2
        ind = np.where(mag>.1) #Automatic vignetting
                            #requires position vector set to 0.
    
    return [rays[i][ind] for i in range(10)]

def lens(rays,r1,r2,thick,d,nl,reverse=False):
    """Trace lens, first surface center is coincident with xy plane
    thickness extends in positive z direction
    Rays should be traced to this plane before calling lens
    to ensure rays hit correct spherical intersection
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(rays,0.,0.,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        sphere(rays,abs(r1)) #Trace to first spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(rays,ind=ind)
    refract(rays,1.,nl) #Refract into lens
    
    transform(rays,0.,0.,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat(rays)
    if r2 != 0:
        transform(rays,0.,0.,r2,0,0,0) #Go to center of curvature
        sphere(rays,abs(r2)) #Trace to second spherical surface
    rho = np.sqrt(x**2 + y**2)/(d/2)
    ind = np.where(rho<=1.) #Remove rays with radius beyond that of lens
    vignette(rays,ind=ind)
    refract(rays,nl,1.) #Refract out of lens
    transform(rays,0.,0.,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return

def cyllens(rays,r1,r2,thick,width,height,nl,reverse=False):
    """Cylindrical lens, same principle as with standard lens
    """
    opd,x,y,z,l,m,n,ux,uy,uz
    if reverse is True:
        r1,r2 = -r2,-r1
    #Trace to first surface
    if r1 != 0:
        transform(rays,0.,0,r1,0.,0.,0.) #Go to first
                                     #center of curvature, r1>0->convex
        cyl(rays,abs(r1)) #Trace to first spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(rays,ind=ind)
    refract(rays,1.,nl) #Refract into lens

    
    transform(rays,0.,0,-r1+thick,0.,0.,0.) #Go to center of second surface
    flat(rays)
    if r2 != 0:
        transform(rays,0.,0,r2,0,0,0) #Go to center of curvature
        cyl(rays,abs(r2)) #Trace to second spherical surface
    ind = logical_and(x<width/2,y<height/2) #Remove rays outside of cylinder
    vignette(rays,ind=ind)
    refract(rays,nl,1.) #Refract out of lens
    transform(rays,0.,0,-r2,0,0,0) #Transform to xy plane tangent to
                              #center of second surface

    #Rays are left at second surface, pointing out of lens in correct direction
    return


###### SOURCES #######

def pointsource(ang,num):
    """Define point source with angular divergence
    Points in +z direction
    """
    #Radial direction cosine magnitude
    rho = np.sqrt(np.random.rand(num))*np.sin(ang)
    theta = np.random.rand(num)*2*np.pi
    l = rho*np.cos(theta)
    m = rho*np.sin(theta)
    n = np.sqrt(1.-l**2-m**2)
    x = np.repeat(0.,num)
    y = np.repeat(0.,num)
    z = np.repeat(0.,num)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return opd,x,y,z,l,m,n,ux,uy,uz

def circularbeam(rad,num):
    """Define uniform, circular beam of radius rad, pointing in +z direction
    """
    rho = np.sqrt(np.random.rand(num))*rad
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return opd,x,y,z,l,m,n,ux,uy,uz

def annulus(rin,rout,num):
    """Define annulus of rays pointing in +z direction
    """
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return opd,x,y,z,l,m,n,ux,uy,uz

def subannulus(rin,rout,dphi,num):
    """Create a subapertured annulus source in +z direction
    Annulus is centered about theta=0 which points to +x
    """
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*dphi - dphi/2.
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return opd,x,y,z,l,m,n,ux,uy,uz

def rectArray(xsize,ysize,num):
    """Creates a regular array of rays using meshgrid and linspace"""
    x,y = np.meshgrid(np.linspace(-xsize,xsize,num),\
                      np.linspace(-ysize,ysize,num))
    opd = np.repeat(0.,num**2)
    x = x.flatten()
    y = y.flatten()
    z = np.copy(opd)
    l = np.repeat(0.,num**2)
    m = np.repeat(0.,num**2)
    n = np.repeat(1.,num**2)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    return opd,x,y,z,l,m,n,ux,uy,uz

def convergingbeam(zset,rin,rout,tmin,tmax,num,lscat):
    """Converging sub-apertured annulus beam
    Place at nominal focus
    Input z position, inner and outer radius,
    min and max theta
    """
    rho = sqrt(rin**2+random.rand(num)*(rout**2-rin**2))
    theta = tmin + random.rand(num)*(tmax-tmin)
    x = rho*cos(theta)
    y = rho*sin(theta)
    z = np.repeat(zset,num)
    lscat = lscat * tan((random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -cos(arctan(rho/zset)+lscat)
    l = -sqrt(1-n**2)*cos(theta)
    m = -sqrt(1-n**2)*sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return opd,x,y,z,l,m,n,ux,uy,uz

def convergingbeam2(zset,xmin,xmax,ymin,ymax,num,lscat):
    """Rectangular converging beam
    Place at nominal focus
    Input z position and rectangular bounds
    """
    x = xmin + random.rand(num)*(xmax-xmin)
    y = ymin + random.rand(num)*(ymax-ymin)
    rho = sqrt(x**2+y**2)
    theta = arctan2(y,x)
    z = np.repeat(zset,num)
    lscat = lscat * tan((random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -cos(arctan(rho/zset)+lscat)
    l = -sqrt(1-n**2)*cos(theta)
    m = -sqrt(1-n**2)*sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return opd,x,y,z,l,m,n,ux,uy,uz

def rectbeam(xhalfwidth,yhalfwidth,num):
    """Rectangular beam pointing in +z direction
    """
    x = (np.random.rand(num)-.5)*2*xhalfwidth
    y = (np.random.rand(num)-.5)*2*yhalfwidth
    z = np.repeat(0.,num)
    n = np.repeat(1.,num)
    l = np.copy(z)
    m = np.copy(z)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return opd,x,y,z,l,m,n,ux,uy,uz
    
####### Lenses ##########
def LJ1653L2(rays,reverse=False):
    cyllens(rays,103.36,0,4.09,30.,60.,1.51501,reverse=reverse)
    return

def LJ1629L2(rays,reverse=False):
    cyllens(rays,77.52,0,4.46,30.,60.,1.51501,reverse=reverse)
    return

def AC254_400_A(rays,reverse=False):
    if reverse is False:
        lens(rays,738.5,181.55,2.,12.7*2,1.64363)
        lens(rays,181.55,-219.8,4.,12.7*2,1.51501)
    else:
        lens(rays,181.55,-219.8,4.,12.7*2,1.51501,reverse=True)
        lens(rays,738.5,181.55,2.,12.7*2,1.64363,reverse=True)
    return

def AC508_200_A(rays,reverse=False):
    if reverse is False:
        lens(rays,109.86,-93.110,8.5,50.8,1.51501)
        lens(rays,-93.110,-376.25,2.,50.8,1.64363)
    else:
        lens(rays,-93.110,-376.25,2.,50.8,1.64363,reverse=True)
        lens(rays,109.86,-93.110,8.5,50.8,1.51501,reverse=True)
    return

def cylNull(rays,reverse=False):
    if reverse is False:
        transform(rays,0,0,0,np.pi/2,0,0)
        cylconic(rays,.007626,-.575)
        refract(rays,1.,1.51501)
        transform(rays,0,0,0,-np.pi/2,0,0)
        transform(rays,0,0,50,0,0,0)
        flat(rays,)
        refract(rays,1.51501,1.)
    else:
        refract(rays,1.,1.51501)
        transform(rays,0,0,50,0,0,0)
        transform(rays,0,0,0,np.pi/2,0,0)
        cylconic(rays,-.007626,-.575)
        refract(rays,1.51501,1.)
        transform(rays,0,0,0,-np.pi/2,0,0)
    return


#######  ANALYSES #########
def centroid(rays,weights=None):
    """Compute the centroid of the rays in the xy plane
    """
    x,y = rays[1:3]
    cx = np.average(x,weights=weights)
    cy = np.average(y,weights=weights)
    return cx,cy

def rmsCentroid(rays,weights=None):
    """Compute the RMS of rays from centroid in xy plane
    """
    x,y = rays[1:3]
    cx,cy = centroid(rays,weights=weights)
    rho = (x-cx)**2 + (y-cy)**2
    return np.sqrt(np.average(rho,weights=weights))

def rmsX(rays,weights=None):
    """RMS from centroid in the X direction"""
    x = rays[1]
    cx = np.average(x,weights=weights)
    rmsx = np.sqrt(np.average((x-cx)**2,weights=weights))
    return rmsx

def rmsY(rays,weights=None):
    """RMS from centroid in the Y direction"""
    y = rays[2]
    cy = np.average(y,weights=weights)
    rmsy = np.sqrt(np.average((y-cy)**2,weights=weights))
    return rmsy

def hpd(rays,weights=None):
    """Compute HPD by taking median of radii from centroid"""
    x,y = rays[1:3]
    cx,cy = centroid(rays,weights=weights)
    rho = np.sqrt((x-cx)**2 + (y-cy)**2)
    if weights is not None:
        ind = np.argsort(rho)
        weights = weights[ind]
        rho = rho[ind]
        cdf = np.cumsum(weights)
        hpd = rho[np.argmin(np.abs(cdf-.5))] * 2.
    else:
        hpd = np.median(rho)*2.
    return hpd

def hpdY(rays,weights=None):
    """Compute HPD in y direction by taking median of radii from centroid
    Does rho need to be absolute value???
    """
    y = rays[2]
    cy = np.average(PT.y,weights=weights)
    rho = np.abs(y-cy)
    if weights is not None:
        ind = np.argsort(rho)
        weights = weights[ind]
        rho = rho[ind]
        cdf = np.cumsum(weights)
        hpd = rho[np.argmin(np.abs(cdf-.5))] * 2.
    else:
        hpd = np.median(rho)*2.
    return hpd

def findimageplane(rays,zscan,num,weights=None):
    """Scans the RMS radius vs z coordinate.
    Given the analytical best focal plane technique, this
    is largely useless. Perhaps still useful for
    investigating focal depth."""
    rms = []
    zsteps = np.linspace(-zscan,zscan,num)
    for znow in np.linspace(-zscan,zscan,num):
        #Transform to offset
        transform(rays,0,0,znow,0,0,0)
        #Trace rays to new plane
        flat(rays)
        rms.append(rmsCentroid(rays,weights=weights))
        #Return to nominal plane
        transform(rays,0,0,-znow,0,0,0)
    flat(rays)

##    plt.clf()
##    plt.plot(zsteps,rms)

    return zsteps[np.where(rms==np.min(rms))]

def findlineplane(rays,zscan,num,weights=None):
    """Scans for minimum rms radius in the Y direction.
    Largely useless given analytical technique."""
    rms = []
    zsteps = np.linspace(-zscan,zscan,num)
    for znow in zsteps:
        #Transform to offset
        transform(rays,0,0,znow,0,0,0)
        #Trace rays to new plane
        flat(rays)
        #Determine centroid
        rms.append(rmsY(rays,weights=weights))
        #Return to nominal plane
        transform(rays,0,0,-znow,0,0,0)
    flat(rays)

    return zsteps[np.where(rms==min(rms))]

def findimageplane2(rays,zscan,num):
    """Scans for minimum HEW in the Y direction.
    Largely useless given analytical technique."""
    hew = []
    zsteps = linspace(-zscan,zscan,num)
    for znow in linspace(-zscan,zscan,num):
        #Transform to offset
        transform(rays,0,0,znow,0,0,0)
        #Trace rays to new plane
        flat(rays,)
        #Append FoM array
        hew.append(hpd(rays))
        #Return to nominal plane
        transform(rays,0,0,-znow,0,0,0)
    flat(rays)

    clf()
    plot(zsteps,hew)

    return zsteps[where(hew==min(hew))]

def analyticImagePlane(rays,weights=None):
    """Find the image plane using the analytic method from
    Ron Elsner's paper
    """
    x,y,z,l,m,n = rays[1:7]
    bx = np.average(x*l/n,weights=weights)-np.average(x,weights=weights)\
         *np.average(l/n,weights=weights)
    ax = np.average((l/n)**2,weights=weights)\
         -np.average(l/n,weights=weights)**2
    by = np.average(y*m/n,weights=weights)-np.average(y,weights=weights)\
         *np.average(m/n,weights=weights)
    ay = np.average((m/n)**2,weights=weights)\
         -np.average(m/n,weights=weights)**2
    dz = -(bx+by)/(ax+ay)
    
    return dz

def analyticYPlane(rays,weights=None):
    """Find the line plane using analytic method from
    Ron Elsner's paper"""
    x,y,z,l,m,n = rays[1:7]
    by = np.average(y*m/n,weights=weights)-np.average(y,weights=weights)\
         *np.average(m/n,weights=weights)
    ay = np.average((m/n)**2,weights=weights)\
         -np.average(m/n,weights=weights)**2
    dz = -by/ay
    return dz

def grazeAngle(rays):
    """Find the graze angle of the rays with the current
    surface normal."""
    return np.arcsin(rays[4]*rays[7] +\
                     rays[5]*rays[8] +\
                     rays[6]*rays[9])

def wsPrimRad(z,psi,r0,z0):
    """Return the radius of a WS primary as a function of axial coordinate
    This is computed numerically by tracing a single ray in plane
    orthogonal to optical axis
    """
    #Set up source pointing toward +z
    rays = pointsource(0.,1)
    transform(rays,0,0,0,0,-np.pi/2,0) #Point ray to +x
    transform(rays,-r0,0,-z,0,0,0) #Go to proper axial locatio

    #Trace to WS primary
    wsPrimary(rays,r0,z0,psi)

    return rays[1][0]

def wsSecRad(z,psi,r0,z0):
    """Return the radius of a WS primary as a function of axial coordinate
    This is computed numerically by tracing a single ray in plane
    orthogonal to optical axis
    """
    #Set up source pointing toward +z
    pointsource(rays,0.,1)
    transform(rays,0,0,0,0,-np.pi/2,0) #Point ray to +x
    transform(rays,-r0,0,-z,0,0,0) #Go to proper axial location

    #Trace to WS primary
    wsSecondary(rays,r0,z0,psi)

    return rays[1][0]

def referencedWavefront(xang,yang,phase,xang2,yang2,phase2):
    """Wrapper for reconstructing referenced wavefront"""
    phaseinf = np.copy(phase)
    phaseinf[:,:] = 0.
    ind = where(logical_or(phase==100,phase2==100))
    phaseinf[ind] = 100
    xanginf = np.copy(xang)
    xanginf = xang2-xang
    yanginf = np.copy(yang)
    yanginf = yang2-yang
    xanginf[ind] = 100
    yanginf[ind] = 100

    #Reconstruct influence wavefront
    influence = reconstruct.reconstruct(xanginf,yanginf,1.e-12,phaseinf)

    #Make invalid np.pixels NaNs
    ind = where(influence==100.)
    influence[ind] = NaN

    return influence

################CUDA ROUTINES################3
##def transferToGPU():
##    """Get rays transfered to the GPU shared memory"""
##    xg = cuda.to_device(x)
##    yg = cuda.to_device(y)
##    zg = cuda.to_device(z)
##    lg = cuda.to_device(l)
##    mg = cuda.to_device(m)
##    ng = cuda.to_device(n)
##    uxg = cuda.to_device(ux)
##    uyg = cuda.to_device(uy)
##    uzg = cuda.to_device(uz)
##    return [xg,yg,zg,lg,mg,ng,uxg,uyg,uzg]
##
##def returnFromGPU(xg,yg,zg,lg,mg,ng,uxg,uyg,uzg):
##    """Get rays from GPU shared memory back to host"""
##    xg.copy_to_host(x)
##    yg.copy_to_host(y)
##    zg.copy_to_host(z)
##    lg.copy_to_host(l)
##    mg.copy_to_host(m)
##    ng.copy_to_host(n)
##    uxg.copy_to_host(ux)
##    uyg.copy_to_host(uy)
##    uzg.copy_to_host(uz)
