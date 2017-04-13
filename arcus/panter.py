#Trace a Wolter I mirror pair feeding a radial grating
#Add scatter to the rays coming out of the mirror pair
#Determine how the astigmatism goes, and how it relates
#to the inherent scatter.

import numpy as np
import matplotlib.pyplot as plt
import traces.PyTrace as PT
import pdb

def setupSource(wave,N,radius=450.,focal=8000.,scatter=50e-6,\
                gwidth=25.,glength=32.,pitch=0.,yaw=0.,roll=0.):
    """Set up converging beam with some scatter.
    Idea is to place a grating at radius and focal length
    with an incidence angle of 1.5 degrees
    """
    #Setup rays over grating surface
    PT.x = np.random.uniform(low=-gwidth/2,high=gwidth/2,size=N)
    PT.y = np.random.uniform(low=-glength/2,high=glength/2,size=N)
    PT.z = np.repeat(0.,N)
    PT.l = np.zeros(N)
    PT.m = np.zeros(N)
    PT.n = np.zeros(N)
    PT.ux = np.zeros(N)
    PT.uy = np.zeros(N)
    PT.uz =  np.repeat(1.,N)
    #Transform to nominal focus to set ray cosines
    PT.transform(0,0,0,-88.5*np.pi/180,0,0)
    PT.transform(0,radius,-focal,0,0,0)
    rad = np.sqrt(PT.x**2+PT.y**2+PT.z**2)
    PT.l = -PT.x/rad
    PT.m = -PT.y/rad
    PT.n = -PT.z/rad
    #Add scatter
    PT.l = PT.l + scatter/10.*np.random.normal(size=N)
    PT.m = PT.m + scatter*np.random.normal(size=N)
    PT.n = -np.sqrt(1. - PT.l**2 - PT.m**2)
    #Go back to plane of of grating
    PT.transform(0,-radius,focal,0,0,0)
    PT.transform(0,0,0,88.5*np.pi/180,0,np.pi)
    #Place grating and diffract
    hubdist = np.sqrt(radius**2+focal**2)*np.cos(1.5*np.pi/180)
    PT.flat()
    PT.reflect()
    PT.transform(0,0,0,pitch,roll,yaw)
    PT.radgrat(hubdist,160./hubdist,1,wave)
    PT.itransform(0,0,0,pitch,roll,yaw)
    #Go to focal plane
    PT.transform(0,0,0,np.pi/2,0,0)
    PT.transform(0,0,-hubdist,0,0,0)
    PT.transform(0,0,27.27727728-.04904905,0,0,0)
    PT.flat()

    #plt.plot(PT.x,PT.y,'.')
    
    return PT.centroid()

def computeYawShift(wave,order,yaw1,yaw2):
    cx1,cy1 = setupSource(wave*order,1000,focal=8052,yaw=yaw1)
    cx2,cy2 = setupSource(wave*order,1000,focal=8052,yaw=yaw2)
    return cy2-cy1
