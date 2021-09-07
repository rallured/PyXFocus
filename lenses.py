#Module to define singlet and doublet lenses
#Explicitly put in commonly used lenses such as Thorlabs doublets
import surfaces as surf
import transformMod as tran
import pdb

def singlet(rays,r1,r2,thick,nl,reverse=False):
    """Trace a spherical singlet lens. Assume reference frame
    is +z toward optical axis, xy plane tangent to first surface.
    Positive R indicates convex surface for both surfaces.
    """
    if reverse is True:
        r1,r2 = r2,r1
    #Trace to first surface
    tran.transform(rays,0,0,r1,0,0,0)
    surf.sphere(rays,r1,nr=1.)
    #Refract into material
    tran.refract(rays,1.,nl)
    #Trace to second surface
    tran.transform(rays,0,0,-r1+thick-r2,0,0,0)
    surf.sphere(rays,r2,nr=nl)
    tran.transform(rays,0,0,r2,0,0,0)
    #Refract out of surface
    tran.refract(rays,nl,1.)
    #Leave at tangent plane of surface
    surf.flat(rays,nr=1.)
    return

def singletCyl(rays,r1,r2,thick,nl,reverse=False):
    """Trace a cylindrical singlet lens. Assume reference frame
    is +z toward optical axis, xy plane tangent to first surface.
    Positive R indicates convex surface for both surfaces.
    Cylindrical axis is in y direction
    Radius of 0 is flat
    """
    if reverse is True:
        r1,r2 = r2,r1
    #Trace to first surface
    tran.transform(rays,0,0,r1,0,0,0)
    if r1==0:
        surf.flat(rays,nr=1.)
    else:
        surf.cyl(rays,r1,nr=1.)
    #Refract into material
    tran.refract(rays,1.,nl)
    #Trace to second surface
    tran.transform(rays,0,0,-r1+thick-r2,0,0,0)
    if r2==0:
        surf.flat(rays,nr=1.)
    else:
        surf.cyl(rays,r2,nr=nl)
    tran.transform(rays,0,0,r2,0,0,0)
    #Refract out of surface
    tran.refract(rays,nl,1.)
    #Leave at tangent plane of last surface
    surf.flat(rays,nr=1.)
    return

def doublet(rays,r1,r2,r3,n1,n2,t1,t2,reverse=False):
    """
    Trace rays through a cemented doublet lens,
    similar to Thorlabs Achromats
    r1,r2,t1,n1 belong to first singlet
    r2,r3,t2,n2 belong to second
    Positive r is convex for end surfaces
    r2 is same convention as r1
    """
    if reverse is True:
        r1,r3,r2 = r3,r1,-r2
        n1,n2 = n2,n1
        t1,t2 = t2,t1
    #Trace to first surface
    tran.transform(rays,0,0,r1,0,0,0)
    surf.sphere(rays,r1,nr=1.)
    #Refract into material
    tran.refract(rays,1.,n1)
    #Trace to second surface
    tran.transform(rays,0,0,-r1+t1+r2,0,0,0)
    surf.sphere(rays,r2,nr=n1)
    #Refract into second material
    tran.refract(rays,n1,n2)
    #Trace to third surface
    tran.transform(rays,0,0,-r2+t2-r3,0,0,0)
    surf.sphere(rays,r3,nr=n2)
    tran.transform(rays,0,0,r3,0,0,0)
    #Refract back to air
    tran.refract(rays,n2,1.)
    #Leave at tangent plane of last surface
    surf.flat(rays,nr=1.)
    
    return
    

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


####### Lenses ##########
def AC508_250(rays,reverse=False):
    """
    Trace rays through AC508-250.
    Assumed highly curved surface first
    """
    doublet(rays,137.1,-111.7,459.2,1.5150885,1.6437928,7.5,2.,\
            reverse=reverse)
    return

def collimator6(rays,reverse=False):
    """
    Traces through the six inch collimator from Cumberland.
    R1=1124.
    R2=9324.
    Standard orientation is collimation of point source
    Reverse is focusing of plane wave
    """
    singlet(rays,9324.,1124.,20.,1.5150885,reverse=reverse)
    return

def edmundCollimator(rays,reverse=False):
    """
    Traces through the 5" collimator from Edmund.
    Standard orientation is  collimation of point source
    Reverse is focusing of plane wave
    """
    if reverse is True:
        singlet(rays,1131.72,780.87,15.42,1.51501)
        tran.transform(rays,0,0,.1,0,0,0)
        singlet(rays,-779.37,2704.01,10.92,1.64363)
    else:
        singlet(rays,2704.01,-779.37,10.92,1.64363)
        tran.transform(rays,0,0,.1,0,0,0)
        singlet(rays,780.87,1131.72,15.42,1.51501)
    return
    

def LJ1516_L2(rays,reverse=False):
    singletCyl(rays,0,516.8,3.2,1.5150885,reverse=reverse)
    return

def LJ1144_L2(rays,reverse=False):
    singletCyl(rays,0,258.4,3.4,1.5150885,reverse=reverse)

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
