import numpy as np
import transformationsf as tran
import transformMod as tr
import pdb

def transform(rays,dx,dy,dz,rx,ry,rz,ind=None,coords=None):
    """Coordinate transformation. translations are done first,
    then Rx,Ry,Rz
    coords[0] - global to local rotation only
    coords[1] - global to local rotations and translations
    coords[2] - local to global rotations only
    coords[3] - local to global rotations and translations
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.transform(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,-dx,-dy,-dz,-rx,-ry,-rz)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.transform(x,y,z,l,m,n,ux,uy,uz,-dx,-dy,-dz,-rx,-ry,-rz)

    #Update transformation matrices
    if coords is not None:
        #Define rotation and translation matrices
        rotm = rotationM(rx,ry,rz)
        tranm = translationM(dx,dy,dz)
        rotmi = rotationM(rx,ry,rz,inverse=True)
        tranmi = translationM(-dx,-dy,-dz)
        #Dot rotation into forward transform
        coords[0] = np.dot(rotm,coords[0])
        coords[1] = np.dot(np.dot(rotm,tranm),coords[1])
        coords[2] = np.dot(coords[2],rotmi)
        coords[3] = np.dot(coords[3],np.dot(tranmi,rotmi))
    
    return


def itransform(rays,dx,dy,dz,rx,ry,rz,coords=None,ind=None):
    """Inverse of coordinate transformations. -rz,-ry,-rx then
    translations.
    """
    x,y,z,l,m,n,ux,uy,uz = rays[1:]
    #tran.itransform(x,y,z,l,m,n,ux,uy,uz,-tx,-ty,-tz,-rx,-ry,-rz)

    if ind is not None:
        tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],\
                                        l[ind],m[ind],n[ind],\
                                        ux[ind],uy[ind],uz[ind]
        tran.itransform(tx,ty,tz,tl,tm,tn,tux,tuy,tuz,-dx,-dy,-dz,-rx,-ry,-rz)
        x[ind],y[ind],z[ind],\
        l[ind],m[ind],n[ind],\
        ux[ind],uy[ind],uz[ind] = tx,ty,tz,tl,tm,tn,tux,tuy,tuz
    else:
        tran.itransform(x,y,z,l,m,n,ux,uy,uz,-dx,-dy,-dz,-rx,-ry,-rz)

    #Update transformation matrices
    if coords is not None:
        #Define rotation and translation matrices
        rotm = rotationM(rx,ry,rz,inverse=True)
        tranm = translationM(-dx,-dy,-dz)
        rotmi = rotationM(rx,ry,rz)
        tranmi = translationM(dx,dy,dz)
        #Dot rotation into forward transform
        coords[0] = np.dot(rotm,coords[0])
        coords[1] = np.dot(np.dot(tranm,rotm),coords[1])
        coords[2] = np.dot(coords[2],rotmi)
        coords[3] = np.dot(coords[3],np.dot(rotmi,tranmi))
    return

def steerY(rays,coords=None):
    """Rotate reference frame for zero mean y tilt"""
    while np.abs(np.mean(rays[5])) > 1e-6:
        transform(rays,0,0,0,-np.mean(rays[5]),0,0,coords=coords)
    return

def steerX(rays,coords=None):
    """Rotate reference frame for zero mean y tilt"""
    while np.abs(np.mean(rays[4])) > 1e-6:
        transform(rays,0,0,0,0,-np.mean(rays[4]),0,coords=coords)
    return

def pointTo(rays,x0,y0,z0,reverse=-1.):
    """
    Direct all ray direction cosines toward (x0,y0,z0)
    reverse=1. will have all rays point away from (x0,y0,z0)
    """
    R = np.sqrt((rays[1]-x0)**2 + (rays[2]-y0)**2 + (rays[3]-z0)**2)
    rays[4] = reverse*(rays[1]-x0)/R
    rays[5] = reverse*(rays[2]-y0)/R
    rays[6] = reverse*(rays[3]-z0)/R
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

def refract(rays,n1,n2):
    """Refract rays based on surface normal
    and ray direction cosines from index n1
    into index n2
    """
    l,m,n,ux,uy,uz = rays[4:]
    tran.refract(l,m,n,ux,uy,uz,n1,n2)
    return

def radgrat(rays,dpermm,order,wave,ind=None):
    """Infinite radial grating. Assumes grating in x,y plane
    with grooves converging at hubdist in positive y direction
    dpermm is nm/mm
    wave is in nm
    """
    x,y,z,l,m,n = rays[1:7]
    #Choose correct radgrat function
    if type(wave) == np.ndarray:
        fn = tran.radgratw
    else:
        fn = tran.radgrat
    if ind is not None:
        tx,ty,tl,tm,tn = x[ind],y[ind],l[ind],m[ind],n[ind]
        if np.size(wave)==1:
            tw = wave
        else:
            tw = wave[ind]
        fn(tx,ty,tl,tm,tn,tw,dpermm,order)
        x[ind],y[ind],l[ind],m[ind],n[ind] = tx,ty,tl,tm,tn
    else:
        fn(x,y,l,m,n,wave,dpermm,order)
    return

def grat(rays,d,order,wave,ind=None):
    """Linear grating with groove direction in +y
    Evanescence results in position vector set to zero
    """
    x,y,z,l,m,n = rays[1:7]

    if ind is not None:
        tx,ty,tl,tm,tn = x[ind],y[ind],l[ind],m[ind],n[ind]
        tran.grat(tx,ty,tl,tm,tn,d,order,wave)
        x[ind],y[ind],l[ind],m[ind],n[ind] = tx,ty,tl,tm,tn
    else:
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

#Transformation matrix helper functions
def rotationM(rx,ry,rz,inverse=False):
    """Return a rotation matrix, applying rotations in
    X,Y,Z order
    Negate the angle values to be consistent with transform function
    Translation translates the reference frame
    """
    if inverse is True:
        rx,ry,rz = -rx,-ry,-rz
    r1 = tr.rotation_matrix(-rx,[1,0,0])
    r2 = tr.rotation_matrix(-ry,[0,1,0])
    r3 = tr.rotation_matrix(-rz,[0,0,1])
    if inverse is True:
        return np.dot(r1,np.dot(r2,r3))
    else:
        return np.dot(r3,np.dot(r2,r1))

def translationM(tx,ty,tz):
    """
    Return a translation matrix. Negate the values in order
    to be consistent with the transform method.
    Translation translates the reference frame"""
    return tr.translation_matrix([-tx,-ty,-tz])

def applyT(rays,coords,inverse=False):
    """Apply transformation matrix to raylist.
    Only rotations to direction cosines.
    Inverse means going back to global coordinate system.
    Forward means going from global coordinate system to
    local coordinate system.
    """
    i = 0
    if inverse is True:
        i = 2
    #Extract position, wavevector, and surface normals
    on = np.ones(np.shape(rays)[1])
    pos = [rays[1],rays[2],rays[3],on]
    wave = [rays[4],rays[5],rays[6],on]
    norm = [rays[7],rays[8],rays[9],on]
    #Apply relevant transformations
    pos = np.dot(coords[i+1],pos)[:3]
    wave = np.dot(coords[i],wave)[:3]
    norm = np.dot(coords[i],norm)[:3]
    #Construct and return new raylist
    return [rays[0],\
            pos[0],pos[1],pos[2],\
            wave[0],wave[1],wave[2],\
            norm[0],norm[1],norm[2]]
    
def applyTPos(x,y,z,coords,inverse=False):
    """Apply transformation to list of points"""
    i = 0
    if inverse is True:
        i = 2
    pos = [x,y,z,np.ones(np.size(x))]
    pos = np.dot(coords[i+1],pos)[:3]
    return pos
