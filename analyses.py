import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import pdb
from utilities.imaging.fitting import circle,circleMerit
from utilities.imaging.analysis import ptov,rms
import utilities.imaging.man as man
import utilities.imaging.fitting as fit
import utilities.imaging.zernikemod as zern
import PyXFocus.reconstruct as reconstruct
import astropy.io.fits as pyfits
import PyXFocus.specialfunctions as special

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

def rmsPoint(rays,point,weights=None):
    """
    Compute the RMS of rays about a point
    """
    if np.size(point) == 10:
        rho = (rays[1]-point[1])**2 + \
              (rays[2]-point[2])**2 + \
              (rays[3]-point[3])**2
    else:
        rho = (rays[1]-point[0])**2 + \
              (rays[2]-point[1])**2 + \
              (rays[3]-point[2])**2
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

def rho(rays,weights=None,cent=False):
    """
    Compute distance from centroid for all rays
    """
    x,y = rays[1:3]
    if cent is True:
        cx,cy = centroid(rays,weights=weights)
    else:
        cx,cy = 0,0
    rho = np.sqrt((x-cx)**2+(y-cy)**2)
    
    return rho

def rhocdf(rays,weights=None,cent=True):
    """
    Compute the radial CDF of the ray distribution
    """
    r = rho(rays,weights=weights,cent=cent)
    if weights is None:
        weights = np.repeat(1,len(r))
    ind = np.argsort(r)
    weights = weights[ind]
    r = r[ind]
    cdf = np.cumsum(weights)
    cdf = cdf / cdf.max()
    
    return r,cdf

def hpd(rays,weights=None):
    """Compute HPD by taking median of radii from centroid"""
    r = rho(rays,weights=weights)
    if weights is not None:
        r,cdf = rhocdf(rays,weights=weights)
        hpd = r[np.argmin(np.abs(cdf-.75))] - \
              r[np.argmin(np.abs(cdf-.25))]
    else:
        hpd = np.median(r)*2.
    return hpd

def hpdY(rays,weights=None):
    """Compute HPD in y direction by taking median of radii from centroid
    Does rho need to be absolute value???
    """
    y = rays[2]
    cy = np.average(y,weights=weights)
    rho = np.abs(y-cy)
    if weights is not None:
        ind = np.argsort(rho)
        weights = weights[ind]
        rho = rho[ind]
        cdf = np.cumsum(weights)
        cdf = cdf / cdf.max()
        hpd = rho[np.argmin(np.abs(cdf-.75))] - \
              rho[np.argmin(np.abs(cdf-.25))]
    else:
        hpd = np.median(rho)*2.
    return hpd

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

def analyticXPlane(rays,weights=None):
    """Find the line plane using analytic method from
    Ron Elsner's paper"""
    x,y,z,l,m,n = rays[1:7]
    bx = np.average(x*l/n,weights=weights)-np.average(x,weights=weights)\
         *np.average(l/n,weights=weights)
    ax = np.average((l/n)**2,weights=weights)\
         -np.average(l/n,weights=weights)**2
    dz = -bx/ax
    return dz

#def grazeAngle(rays,flat=False):
#    """Find the graze angle of the rays with the current
#    surface normal."""
#    return np.arcsin(rays[4]*rays[7] +\
#                     rays[5]*rays[8] +\
#                     rays[6]*rays[9])

def indAngle(rays,ind = None,normal = None):
    """Find the incidence angle of the rays with either the current or a specified
    surface normal."""
    if normal is None:
        if ind is not None:
            wave,x,y,z,l,m,n,ux,uy,uz = rays
            tx,ty,tz,tl,tm,tn,tux,tuy,tuz = x[ind],y[ind],z[ind],l[ind],m[ind],n[ind],ux[ind],uy[ind],uz[ind]
            iangs = np.arccos(tl*tux + tm*tuy + tn*tuz)
        else:
            iangs = np.arccos(rays[4]*rays[7] + rays[5]*rays[8] + rays[6]*rays[9])
    else:
        if ind is not None:
            wave,x,y,z,l,m,n,ux,uy,uz = rays
            dir_cosines = np.vstack((l[ind],m[ind],n[ind]))
            iangs = np.arccos(np.dot(np.array([normal]),dir_cosines))[0]
        else:
            dir_cosines = np.vstack((rays[4],rays[5],rays[6]))
            iangs = np.arccos(np.dot(np.array([normal]),dir_cosines))[0]
    return iangs

def grazeAngle(rays,ind = None):
    """Find the graze angle of the rays with the current
    surface normal."""
    return np.pi/2 - indAngle(rays,ind = ind)

def interpolateVec(rays,I,Nx,Ny,xr=None,yr=None,method='linear',\
                   polar=False,interpVec=None):
    """
    Interpolate a ray vector onto a 2D grid based on the X and Y
    positions of the rays. Assume that the rays randomly fill a
    contiguous region of space. Automatically set up the grid
    using the min/max X/Y positions. User inputs the number of
    points in the grid as Nx,Ny, where total points is Nx*Ny.
    User can also indicate the interpolation method used by
    griddata.
    I indicates which vector to interpolate
    """
    #Unpack needed vectors
    x,y = rays[1:3]
    if interpVec is None:
        interpVec = rays[I]
    
    #Set up new grid
    if xr is None:
        xr=[x.min(),x.max()]
        yr=[y.min(),y.max()]
    gridx,gridy = np.meshgrid(np.linspace(xr[0],xr[1],Nx),\
                              np.linspace(yr[0],yr[1],Ny))
    dx = np.diff(gridx)[0][0]
    dy = np.diff(np.transpose(gridy))[0][0]

    if polar is True:
        #Convert to polar
        rho = np.sqrt(x**2+y**2)
        theta1 = np.arctan2(y,x)
        theta2 = np.arctan2(x,y)
        rhog = np.sqrt(gridx**2+gridy**2)
        azg1 = rhog * np.arctan2(gridy,gridx)
        azg2 = rhog * np.arctan2(gridx,gridy)
        res1 = griddata((rho,theta1*rho),interpVec,(rhog,azg1),method=method)
        res2 = griddata((rho,theta2*rho),interpVec,(rhog,azg2),method=method)
        res = np.nanmedian([res1,res2],axis=0)
    else:
        res = griddata((x,y),interpVec,(gridx,gridy),method=method)
    
    return res,dx,dy

def measureOPD(rays,point):
    """
    Measure distance in mm from rays to point
    Returns a vector of distances
    """
    if len(point)==10:
        return np.sqrt((rays[1]-point[1])**2+\
                       (rays[2]-point[2])**2+\
                       (rays[3]-point[3])**2)
    else:
        return np.sqrt((rays[1]-point[0])**2+\
                       (rays[2]-point[1])**2+\
                       (rays[3]-point[2])**2)

def OPDtoLegendre(x,y,opd,xo,yo,xwidth=1.,ywidth=1.):
    """
    Fit Legendre coefficients over ray XY plane to OPD.
    """
    #Perform Legendre fitting
    res = fit.legendre2d(opd,x=x/xwidth,y=y/ywidth,xo=xo,yo=yo)
    #Construct order arrays
    coeff = res[1]
    [xorder,yorder] = np.meshgrid(range(xo+1),range(yo+1))
    return coeff,xorder,yorder

def OPDtoLegendreP(x,y,opd,xo,yo,xwidth=1.,ywidth=1.):
    """
    Fit Legendre coefficients over ray XY plane to OPD.
    Return the derivatives in both axes
    """
    #Perform Legendre fitting
    res = fit.legendre2d(opd,x=x/xwidth,y=y/ywidth,xo=xo,yo=yo)
    #Construct order arrays
    coeff = res[1]
    [xorder,yorder] = np.meshgrid(range(xo+1),range(yo+1))
    #Construct the derivatives using the special function module
    cx = np.polynomial.legendre.legder(res[1],axis=0)
    cy = np.polynomial.legendre.legder(res[1],axis=1)
    #gy = np.polynomial.legendre.legval2d(x/xwidth,y/ywidth,cx)/xwidth
    #gx = np.polynomial.legendre.legval2d(x/xwidth,y/ywidth,cy)/ywidth
    gy = np.polynomial.legendre.legval2d(-y/ywidth,x/xwidth,cx)/xwidth
    gx = np.polynomial.legendre.legval2d(-y/ywidth,x/xwidth,cy)/ywidth
    
    return gx,gy

def OPDtoSqZernike(x,y,opd,N,xwidth=1.,ywidth=1.):
    """
    Fit Square Zernike coefficients over ray XY plane to OPD.
    """
    #First determine standard Zernike coefficients
    res = zern.fitvec(x/xwidth,y/ywidth,opd,N=N)
    #Convert to square Zernikes using precomputed transformation matrix
    sqcoef = np.dot(np.linalg.inv(zern.sqtost[:N,:N]),res[0])
    return sqcoef

def SqZerniketoOPD(x,y,coef,N,xwidth=1.,ywidth=1.):
    """
    Return an OPD vector based on a set of square Zernike coefficients
    """
    stcoef = np.dot(zern.sqtost[:N,:N],coef)
    x = x/xwidth
    y = y/ywidth
    zm = zern.zmatrix(np.sqrt(x**2+y**2),np.arctan2(y,x),N)
    opd = np.dot(zm,stcoef)
    return opd

def wavefront(rays,Nx,Ny,method='cubic',polar=False):
    """
    Interpolate the beam slopes and then integrate the wavefront.
    Rays assumed to be in the XY plane.
    """
    #Interpolate slopes to regular grid
    y,dx,dy = interpolateVec(rays,5,Nx,Ny,method=method,polar=polar)
    x,dx,dy = interpolateVec(rays,4,Nx,Ny,method=method,polar=polar)
        

    #Prepare arrays for integration
    #Reconstruct requires nans be replaced by 100
    #Fortran functions require arrays to be packed in
    #fortran contiguous mode
    x = man.padRect(x)
    y = man.padRect(y)
    phase = np.zeros(np.shape(x),order='F')
    phase[np.isnan(x)] = 100.
    x[np.isnan(x)] = 100.
    y[np.isnan(y)] = 100.
    y = np.array(y,order='F')
    x = np.array(x,order='F')

    #Reconstruct and remove border    
    phase = reconstruct.reconstruct(y,x,1e-12,dx,phase)
    phase[phase==100] = np.nan
    x[x==100] = np.nan
    y[y==100] = np.nan

    return phase[1:-1,1:-1],x[1:-1,1:-1],y[1:-1,1:-1]

def measurePower(rays,Nx,Ny,method='linear'):
    """Measure the radius of curvature in X and Y
    axes of OPD
    Assumes you have steered the beam to get rid of
    any tilts
    Sign of power equals transform in z to go to focus"""
    #Get l,m
    l,dx,dy = interpolateVec(rays,4,Nx,Ny,method=method)
    m = interpolateVec(rays,5,Nx,Ny,method=method)[0]
    #Get slices
    xsl = man.stripnans(l[Ny/2])
    ysl = man.stripnans(m[:,Nx/2])
    #Estimate gradients
    xpow = 1/np.gradient(xsl,dx)[Nx/2]
    ypow = 1/np.gradient(ysl,dy)[Ny/2]
    return -xpow,-ypow
    

def compareOPDandSlopes(rays,Nx,Ny,method='linear'):
    """
    Function to compare OPD gradient to ray slopes quickly
    """
    #Get opd
    opd,dx = interpolateVec(rays,0,Nx,Ny,method=method)
    #Get gradient
    gradx,grady = np.gradient(opd,dx)
    #Get slopes
    l,dx = interpolateVec(rays,4,Nx,Ny,method=method)
    m,dx = interpolateVec(rays,5,Nx,Ny,method=method)
    #Make plots
    fig = plt.figure()
    
    fig.add_subplot(321)
    plt.imshow(grady)
    plt.title('Gradx')
    plt.colorbar()
    
    fig.add_subplot(322)
    plt.imshow(gradx)
    plt.title('Grady')
    plt.colorbar()

    fig.add_subplot(323)
    plt.imshow(l)
    plt.title('l')
    plt.colorbar()

    fig.add_subplot(324)
    plt.imshow(m)
    plt.title('m')
    plt.colorbar()

    fig.add_subplot(325)
    plt.imshow(l-grady)
    plt.title('l-gradx')
    plt.colorbar()

    fig.add_subplot(326)
    plt.imshow(m-gradx)
    plt.title('m-grady')
    plt.colorbar()

    print 'X Diff:' + str(np.sqrt(np.nanmean((grady-l)**2)))
    print 'Y Diff:' + str(np.sqrt(np.nanmean((gradx-m)**2)))
    
    return

def radialGrad(x,y,hubscale,yaw,hubdist):
    """
    Compute the x and y derivatives of a radial grating phase function.
    hubscale is lambda/(nn/mm)
    Phase derivatives are:
    gx = alpha*(cos(theta)/rho)
    gy = -alpha*(sin(theta)/rho)
    """
    #Rotate coordinate system based on yaw angle
    x2 = x*np.cos(yaw)+y*np.sin(yaw)
    y2 = -x*np.sin(yaw)+y*np.cos(yaw)

    #Shift y2 axis by hubdistance
    y2 = y2+hubdist

    #Compute rho and theta
    rho = np.sqrt(x2**2+y2**2)
    theta = np.arctan2(-x2,y2)
    
    #Compute derivatives
    gx = hubscale * (np.cos(theta)/rho)
    gy = -hubscale * (np.sin(theta)/rho)

    return gx,gy


