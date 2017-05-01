import traces.PyTrace as PT
import numpy as np
import traces.conicsolve as con
import pdb
import scipy.optimize

def traceCyl(align):    
    """Traces a cylindrical approximation to Wolter I geometry
    Assumes 1 m radius of curvature in accordance with OAB samples
    align is a 12 element array giving the transformations to be
    applied to each mirror
    Uses identical set of rays each run defined upon module import
    Adds a restraint such that any vignetting over 25% results in
    a huge merit function
    """
    
    #Set up source
    np.random.seed(5)
    a,p,d,e = con.woltparam(1000.,10000.)
    r0 = con.primrad(10025.,1000.,10000.)
    r1 = con.primrad(10125.,1000.,10000.)
    dphi = 100./1000.
    PT.subannulus(r0,r1,dphi,10**4)
    PT.transform(0,0,0,np.pi,0,0)
    PT.transform(0,0,-10500.,0,0,0)

    #Trace primary cylinder: go to tangent point at center
    #of mirror and then rotate to cone angle, then move
    #to new cylinder axis and tracecyl
    rt = con.primrad(10075.,1000.,10000.)
    PT.transform(rt,0,0,0,a,0)
    PT.transform(*align[:6])
    PT.transform(0,0,0,np.pi/2,0,0)
    PT.transform(-1000.,0,0,0,0,0)
    PT.cyl(1000.)

    #Vignette rays missing physical surface
    ind = np.logical_and(abs(PT.z)<50.,abs(PT.y)<50.)
    PT.vignette(ind=ind)

    #Reflect and reverse transformations
    PT.reflect()
    PT.transform(1000.,0,0,0,0,0)
    PT.itransform(0,0,0,np.pi/2,0,0)
    PT.itransform(*align[:6])
    PT.itransform(rt,0,0,0,a,0)

    #Trace secondary cylinder: same principle as before
    rt = con.secrad(9925.,1000.,10000.)
    PT.transform(0,0,-150.,0,0,0)
    PT.transform(rt,0,0,0,3*a,0)
    PT.transform(*align[6:])
    PT.transform(0,0,0,np.pi/2,0,0)
    PT.transform(-1000.,0,0,0,0,0)
    PT.cyl(1000.)

    #Vignette rays missing physical surface
    ind = np.logical_and(abs(PT.z)<50.,abs(PT.y)<50.)
    PT.vignette(ind=ind)

    #Reflect and reverse transformations
    PT.reflect()
    PT.transform(1000.,0,0,0,0,0)
    PT.itransform(0,0,0,np.pi/2,0,0)
    PT.itransform(*align[6:])
    PT.itransform(rt,0,0,0,3*a,0)

    #Go down to nominal focus
    PT.transform(0,0,-9925.,0,0,0)
    PT.flat()

    #Compute merit function
    nom = PT.rmsCentroid()/10**4*180./np.pi*60**2
    nom = nom + max((9500.-np.size(PT.x)),0.)
    
    return nom

#Attempt to minimize RMS from centroid via alignment of cylinders
def alignCyls():
    #Create function
    start = np.repeat(.0001,12)
    res = scipy.optimize.minimize(traceCyl,start,method='nelder-mead',\
                                  options={'disp':False})

    return res

#Determine axial sag as a function of radius of curvature
def sagVsRadius(rad,z0,z1,z2):
    """Compute axial sag needed as a function of radius of curvature
    Input is vector of radii and upper and lower axial length bounds
    """
    sag = np.zeros(np.size(rad))
    z = np.linspace(z1,z2,1000)
    for r0 in rad:
        r = con.primrad(z,r0,z0)
        r = r - np.polyval(np.polyfit(z,r,1),z)
        sag[r0==rad] = np.max(r)-np.min(r)
    return sag
