#Module to collect routines to trace rays to various surfaces
import numpy as np
import surfacesf as surf
import zernsurf as zern
import zernikemod
import woltsurf as wolt
import traces.transformations as tran
from traces.analyses import analyticYPlane,analyticXPlane,analyticImagePlane
import traces.conicsolve as con
import pdb

def flat(rays,ind=None,nr=None):
    """Trace rays to the XY plane
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if ind is not None:
        #Temporary array
        trays = [rays[i][ind] for i in range(10)]
        #Trace
        surf.flat(*trays[1:])
        #Copy back to original
        for i in range(1,10):
            rays[i][ind] = trays[i]
    elif nr is not None:
        surf.flatopd(x,y,z,l,m,n,ux,uy,uz,opd,nr)
    else:
        surf.flat(x,y,z,l,m,n,ux,uy,uz)
    return

def zernsurf(rays,coeff,rad,rorder=None,aorder=None,nr=None):
    """Wrapper for Zernike surface
    Coordinates are usual arctan2(y,x)
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder is None:
        rorder,aorder = zernikemod.zmodes(np.size(coeff))
    if nr is None:
        zern.tracezern(x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad)
    else:
        zern.tracezernopd(opd,x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad,nr)
##    rho = np.sqrt(x**2+y**2)
##    ind = np.where(rho<=rad)
##    vignette(ind=ind)
    return

def zernphase(rays,coeff,rad,wave,rorder=None,aorder=None):
    """Wrapper for standard Zernike phase surface. Supply
    wavelength in mm, radius in mm, coeff in mm."""
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder is None:
        rorder,aorder = zernikemod.zmodes(np.size(coeff))
    zern.zernphase(opd,x,y,z,l,m,n,ux,uy,uz,coeff,\
                   np.array(rorder),np.array(aorder),rad,wave)
##    rho = np.sqrt(x**2+y**2)
##    ind = np.where(rho<=rad)
##    rays = vignette(rays,ind=ind)
    return

def zernsurfrot(rays,coeff1,coeff2,rad,rot,\
                rorder1=None,aorder1=None,rorder2=None,aorder2=None):
    """Wrapper for Zernike surface with 2 Zernike sets and one with
    arbitrary rotation angle
    Coordinates are usual arctan2(y,x)
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if rorder1 is None:
        rorder1,aorder1 = zernikemod.zmodes(np.size(coeff1))
    if rorder2 is None:
        rorder2,aorder2 = zernikemod.zmodes(np.size(coeff2))
    zern.tracezernrot(x,y,z,l,m,n,ux,uy,uz,coeff1,np.array(rorder1),\
                      np.array(aorder1),coeff2,np.array(rorder2),\
                      np.array(aorder2),rad,rot)
##    rho = np.sqrt(x**2+y**2)
##    ind = np.where(rho<=rad)
##    vignette(ind=ind)
    return

def sphere(rays,rad,nr=None):
    """Wrapper for spherical surface.
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if nr is not None:
        surf.tracesphereopd(opd,x,y,z,l,m,n,ux,uy,uz,rad,nr)
    else:
        surf.tracesphere(x,y,z,l,m,n,ux,uy,uz,rad)
    return

def conic(rays,R,K,nr=None):
    """Wrapper for conic surface with radius of curvature R
    and conic constant K
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if nr is not None:
        surf.conicopd(opd,x,y,z,l,m,n,ux,uy,uz,R,K,nr)
    else:
        surf.conic(x,y,z,l,m,n,ux,uy,uz,R,K)
    return

def torus(rays,rin,rout):
    """Wrapper for toroidal surface. Outer radius
    is in xy plane, inner radius is orthogonal.
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    surf.torus(x,y,z,l,m,n,ux,uy,uz,rin,rout)
    return

def cyl(rays,rad,nr=None):
    """Wrapper for cylindrical surface
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    if nr is not None:
        surf.tracecylopd(opd,x,y,z,l,m,n,ux,uy,uz,rad,nr)
    else:
        surf.tracecyl(x,y,z,l,m,n,ux,uy,uz,rad)
    return

def cylconic(rays,rad,k):
    """Wrapper for cylindrical conics
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    surf.cylconic(x,y,z,l,m,n,ux,uy,uz,rad,k)

def paraxial(rays,F):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    surf.paraxial(x,y,z,l,m,n,ux,uy,uz,F)
    return

def paraxialY(rays,F):
    """
    Trace rays through an ideal, paraxial lens.
    Assume optical axis is at xy=0 in z direction
    Surface is in xy plane
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    surf.paraxialy(x,y,z,l,m,n,ux,uy,uz,F)
    return

def wolterprimary(rays,r0,z0,psi=1.):
    """Wrapper for Wolter primary surface - no vignetting
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    wolt.wolterprimary(x,y,z,l,m,n,ux,uy,uz,r0,z0,psi)
    return

def wolterprimarynode(rays,r0,z0,psi=1.):
    """Place Wolter node at current origin,
    focus at (-r0,0,-z0)
    """
    tran.transform(rays,-r0,0,-z0,0,0,0)
    wolterprimary(rays,r0,z0,psi)
    tran.itransform(rays,-r0,0,-z0,0,0,0)
    return

def woltersecondary(rays,r0,z0,psi=1.):
    """Wrapper for Wolter secondary surface - no vignetting
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    wolt.woltersecondary(x,y,z,l,m,n,ux,uy,uz,r0,z0,psi)
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
    wolt.woltersine(x,y,z,l,m,n,ux,uy,uz,r0,z0,amp,freq)
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

def secondaryLL(rays,r0,z0,psi,zmax,zmin,dphi,coeff,axial,az):
    """Wrapper for L-L secondary surface
    Placed at focus
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    wolt.woltersecll(x,y,z,l,m,n,ux,uy,uz,r0,z0,psi,\
                     zmax,zmin,dphi,coeff,axial,az)
    return

def primaryLL(rays,r0,z0,zmax,zmin,dphi,coeff,axial,az):
    """Wrapper for L-L primary surface
    Placed at focus
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    wolt.wolterprimll(x,y,z,l,m,n,ux,uy,uz,r0,z0,zmax,zmin,dphi,coeff,axial,az)
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
    wolt.wsprimary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
    return

def wsSecondary(rays,r0,z0,psi):
    """Trace a W-S secondary surface
    Fortran function computes Chase parameters for an equivalent W-I
    betas, f, g, and k computed from alpha and z0
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    a,p,d,e = con.woltparam(r0,z0)
    wolt.wssecondary(x,y,z,l,m,n,ux,uy,uz,a,z0,psi)
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
        wolt.spocone(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,R0,tg)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        wolt.spocone(x,y,z,l,m,n,ux,uy,uz,R0,tg)
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

def ellipsoidPrimary(rays,R0,F,S,psi):
    """
    Trace rays to the primary of an ellipsoid-hyperboloid
    telescope.
    Call at focus, just like with Wolter-I
    """
    #Compute ellipsoid parameters
    P,a,b,e,f = con.ellipsoidFunction(S,psi,R0,F)
    R = b**2/a

    #Move to vertex
    tran.transform(rays,0,0,F+f-P-a,0,0,0)
    #Call conic
    conic(rays,R,-e**2)
    #Go back to focus
    tran.itransform(rays,0,0,F+f-P-a,0,0,0)

    return

def ellipsoidSecondary(rays,R0,F,S,psi):
    """
    Trays rays to the secondary of an ellipsoid-hyperboloid
    telescope.
    Call at focus, just like with Wolter-I.
    Effective psi for secondary must be computed from
    ellipsoid parameters.
    """
    #Compute ellipsoid parameters
    P,a,b,e,f = con.ellipsoidFunction(S,psi,R0,F)
    psi_eff = np.arctan(R0/P)/(np.arctan(R0/F)-np.arctan(R0/P))
    #Call Wolter secondary
    woltersecondary(rays,R0,F,psi=psi_eff)
    return

def ellipsoidPrimaryLL(rays,R0,F,S,psi,zmax,zmin,dphi,coeff,axial,az):
    """
    Trace rays to the primary of an ellipsoid-hyperboloid
    telescope. Add L-L distortions.
    Call at focus, just like with Wolter-I
    """
    opd,x,y,z,l,m,n,ux,uy,uz = rays
    wolt.ellipsoidwoltll(x,y,z,l,m,n,ux,uy,uz,R0,F,psi,S,\
                         zmax,zmin,dphi,coeff,axial,az)

    return

def ellipsoidSecondaryLL(rays,R0,F,S,psi,zmax,zmin,dphi,coeff,axial,az):
    """
    Trace rays to the secondary of an ellipsoid-hyperboloid
    telescope. Add L-L distortions.
    Call at focus, just like with Wolter-I
    """
    #Compute ellipsoid parameters
    P,a,b,e,f = con.ellipsoidFunction(S,psi,R0,F)
    psi_eff = np.arctan(R0/P)/(np.arctan(R0/F)-np.arctan(R0/P))
    #Call secondary
    secondaryLL(rays,R0,F,psi_eff,zmax,zmin,dphi,coeff,axial,az)
    return

def focus(rays,fn,weights=None,nr=None,coords=None):
    dz1 = fn(rays,weights=weights)
    tran.transform(rays,0,0,dz1,0,0,0,coords=coords)
    flat(rays,nr=nr)
    dz2 = fn(rays,weights=weights)
    tran.transform(rays,0,0,dz2,0,0,0,coords=coords)
    flat(rays,nr=nr)
    
    return dz1+dz2

def focusY(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticYPlane,weights=weights,nr=nr,coords=coords)

def focusX(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticXPlane,weights=weights,nr=nr,coords=coords)

def focusI(rays,weights=None,nr=None,coords=None):
    return focus(rays,analyticImagePlane,weights=weights,nr=nr,coords=coords)
