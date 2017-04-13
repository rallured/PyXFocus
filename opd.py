import traces.PyTrace as PT
import numpy as np

def testOPD(N):
    """Test out the OPD wavefront formalism...does it match
    wavefront slope?"""
    rays = PT.rectArray(10.,10.,N)
    opd = np.reshape(rays[0],(N,N))
    l = np.reshape(rays[4],(N,N))
    x = np.reshape(rays[1],(N,N))
    PT.transform(rays,0,0,1000,0,0,0)
    PT.sphere(rays,1000.,nr=1)
    PT.refract(rays,1.,1.51501)
    #Magnify OPD based on index change
    #opd = opd*1./1.51501
    PT.transform(rays,0,0,-1000+10,0,0,0)
    PT.flat(rays,nr=1.51501)
    return rays,opd,l,x
    PT.refract(rays,1.51501,1.)
    PT.transform(rays,0,0,1700,0,0,0)
    PT.flat(rays,nr=1.)
    
    #Compute gradient of OPD
    opd = np.reshape(rays[0],(100,100))
    x = np.reshape(rays[1],(100,100))
    dx = np.diff(x)[0][0]
    l = np.reshape(rays[4],(100,100))
    grady,gradx = np.gradient(opd,dx)

    return rays
