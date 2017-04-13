import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pdb,os,pyfits
import utilities.imaging.fitting as fit
import traces.sources as sources
import traces.transformations as tran
import traces.surfaces as surf
import traces.analyses as anal

#Get distortion coefficients
os.chdir('/home/rallured/Dropbox/AXRO/Alignment/AlignmentStressCompensation/')
res = pyfits.getdata('PrimaryCoeffs.fits')
res[0][0] = 0.
res2 = pyfits.getdata('SecondaryCoeffs.fits')
res2[0][0] = 0.

def findDistortionCoefficients(filename,Nx,Ny,method='cubic'):
    """
    Fit L-L coefficients to distortion data produced by Vanessa.
    """
    #Load in data
    d = np.transpose(np.genfromtxt(filename,skip_header=1,delimiter=','))

    #Compute angle and radial perturbations
    x0,y0,z0 = d[2:5]
    t0 = np.arctan2(x0,-z0)
    r0 = 1 - np.sqrt(x0**2+z0**2)

    #Define regular grid
    tg,zg = np.meshgrid(np.linspace(t0.min(),t0.max(),Nx+2),\
                        np.linspace(y0.min(),y0.max(),Ny+2))

    #Regrid nominal node positions
    d0 = griddata((t0,y0),r0,(tg,zg),method=method)

    #Find perturbed node positions
    x1,y1,z1 = x0+d[5],y0+d[6],z0+d[7]
    t1 = np.arctan2(x1,-z1)
    r1 = 1 - np.sqrt(x1**2+z1**2)
    #Regrid distorted nodes onto original grid
    d1 = griddata((t1,y1),r1,(tg,zg),method=method)

    #Get distortion data
    #
    d = d1-d0

    #Apply fit
    res = fit.legendre2d(d,xo=10,yo=10)
    xo,yo = np.meshgrid(np.arange(11),np.arange(11))

    pdb.set_trace()
    
    #Return coefficients and order arrays    
    return res[1].flatten(),xo.flatten(),yo.flatten()

def pairRaytrace(secondaryTilt,despace):
    """Trace the distorted mirror pair. Assume no gap for now.
    Vignette rays that land outside active mirror areas."""
    #Define ray subannulus
    r1 = surf.con.primrad(8600.,1000.,8400.)
    ang = 260./1000. #arc length over radius is angular extent
    rays = sources.subannulus(1000.,r1,ang,10**3)
    tran.transform(rays,0,0,0,np.   pi,0,0) #Point in -z
    tran.transform(rays,0,0,-10000,0,0,0) #Converge from above

    #Trace to primary
    surf.primaryLL(rays,1000.,8400.,8600,8400,ang,res[0],res[2],res[1])
    #Vignette rays missing
    ind = np.logical_and(rays[3]<8600.,rays[3]>8400.)
    rays = tran.vignette(rays,ind)
    numin = float(len(rays[1]))
    #Reflect
    tran.reflect(rays)
    #Bring to midplane of despaced secondary
    
    #Apply secondary misalignment
    tran.transform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    tran.transform(rays,0,0,despace,0,secondaryTilt,0)
    tran.itransform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    #Trace to secondary
    surf.secondaryLL(rays,1000.,8400.,8400.,8200.,ang,res2[0],res2[2],res2[1])
    #Vignette rays missing
    ind = np.logical_and(rays[3]<8400.,rays[3]>8200.)
    rays = tran.vignette(rays,ind)
    numout = float(len(rays[1]))
    #Reflect
    tran.reflect(rays)
    #Reverse secondary misalignment
    tran.transform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)
    tran.itransform(rays,0,0,despace,0,secondaryTilt,0)
    tran.itransform(rays,surf.con.secrad(8300.,1000.,8400.),0,8300,0,0,0)

    #Go to focus
    surf.focusI(rays)

    #Get centroid
    cx,cy = anal.centroid(rays)

    #Return merit function
    return anal.rmsCentroid(rays)/8400.*180/np.pi*60**2,numout/numin,cx,cy
