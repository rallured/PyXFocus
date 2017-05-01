import numpy as np
import matplotlib.pyplot as plt
from numpy import exp,cos,sin,pi,sqrt
import traces.transformations as tran
import traces.sources as sources
import traces.lenses as lenses
import traces.analyses as anal
import traces.surfaces as surf
import pdb
import copy
import utilities.imaging.fitting as fit
from utilities.imaging.analysis import ptov,rms

nbk7 = 1.5150885
nsf2 = 1.6437928
nsil = 1.45701735

#Get CGH coefficients
cghcoeff = np.genfromtxt('/home/rallured/Dropbox/'
                         'AXRO/Metrology/200mmCGH.txt')[2:]
cgh1m = np.genfromtxt('/home/rallured/Dropbox/AXRO/'
                      'Metrology/CylinderCGH_F9_110mmx110mm_null.txt')[2:]
#Empirically determined line focus distance
line = 199.997297
line1m = 999.95158027583807
#Empirically determined focus from field lens
foc = 264.92386077410214#275.01818439272461
foc1m = 236.25530923977206
#Empriically determined position of cylindrical field lens
cylz = 120.70707070707071
cylz1m = 173.73737373737373

def focusCyl(cylz,div=.1*pi/180,rad=220.):
    raylist = cylindricalSource(N=1000,div=div,rad=rad)
    rmsy = [backToWFS1m(rays,cylz) for rays in raylist]
    return np.mean(rmsy)


def rayBundle(N,div,az,height,rad):
    """
    Set up a diverging ray bundle on the 220 mm cylinder.
    """
    #Establish rays
    rays = sources.pointsource(div,N)
    #Go to cylindrical axis
    tran.transform(rays,0,0,rad,0,0,0)
    #Apply height offset
    tran.transform(rays,-height,0,0,0,0,0)
    #Apply azimuthal offset
    tran.transform(rays,0,0,0,-az,0,0)
    #Go back to tangent plane
    tran.transform(rays,0,0,-rad,0,0,0)
    surf.flat(rays,nr=1.)
    return rays

def cylindricalSource(height=np.linspace(-50.,50.,3),\
                      az=np.linspace(-pi/12,pi/12,5),\
                      N=100,div=pi/180.,rad=220.):
    """
    Set up source rays for focus testing. Set them up on the
    220 mm nominal cylindrical surface using coordinate transforms.
    Use extreme heights and azimuths as well as central sources.
    Use N rays per position, with div divergence
    """
    raylist = [rayBundle(N,div,a,h,rad) for a in az for h in height]
    return raylist

def investigateFocus(focErr=np.linspace(.1,1,.100),app=75.):
    """
    Trace various focus errors, fit the Legendres,
    and return the ratios of the quadratic and quartic
    Legendre coefficients
    """
    coeffs = [determineFocusCoefficients(f,app=app) for f in focErr]

    return np.transpose(coeffs)

def investigateTilt(tiltErr=np.linspace(1e-5,1e-3,100),app=75.):
    """
    Trace various tilt errors, fit the Legendres,
    and return the coefficients
    """
    coeffs = [determineTiltCoefficients(f,app=app) for f in tiltErr]
    return np.transpose(coeffs)

def investigateTip(tipErr=np.linspace(1e-5,1e-3,100),app=75.):
    """
    Trace various tip errors, fit the Legendres,
    and return the coefficients
    """
    coeffs = [determineTipCoefficients(f,app=app) for f in tipErr]
    return np.transpose(coeffs)

def investigateTwist(twistErr=np.linspace(1e-5,1e-3,100),app=75.):
    """
    Trace various twist errors, fit the Legendres,
    and return the coefficients
    """
    coeffs = [determineTwistCoefficients(f,app=app) for f in twistErr]
    return np.transpose(coeffs)

def determineFocusCoefficients(focerror,app=75.):
    d = misalignmentTerm(10000,[0,0,focerror,0,0,0],app=app)
    return fitFocus(d)

def determineTiltCoefficients(tilterror,app=75.):
    d = misalignmentTerm(10000,[0,0,0,tilterror,0,0],app=app)
    return fitTilt(d)

def determineTipCoefficients(tiperror,app=75.):
    d = misalignmentTerm(10000,[0,0,0,0,tiperror,0],app=app)
    return fitTip(d)

def determineTwistCoefficients(twisterror,app=75.):
    d = misalignmentTerm(10000,[0,0,0,0,0,twisterror],app=app)
    return fitTwist(d)
    

def fitFocus(d):
    """
    Fit up to 4th order for focus and return the
    2nd and 4th order coefficients
    """
    res,c = fit.legendre2d(d,xo=0,yo=4,xl=[1,0,0,0,0],yl=[1,0,1,2,4])
    return ptov(d),ptov(d-res),c[2,0],c[4,0],res

def fitTilt(d):
    """
    Fit the tilt term (pitch) using the 1st and 3rd
    order coefficients
    """
    res,c = fit.legendre2d(d,xo=0,yo=3,xl=[1,0,0,0],yl=[0,0,1,3])
    return ptov(d),ptov(d-res),c[1,0],c[3,0],res

def fitTip(d):
    """
    Fit the tip term (roll) using the 1st and
    2,1 coefficients
    """
    res,c = fit.legendre2d(d,xo=1,yo=2,xl=[1,0,1,1],yl=[0,0,0,2])
    return ptov(d),ptov(d-res),c[0,1],c[2,1],res

def fitTwist(d):
    """
    Fit the twist term (yaw) using the 1,1 and
    3,1 terms
    """
    res,c = fit.legendre2d(d,xo=1,yo=2,xl=[1,0,1,1,1],yl=[0,0,0,1,3])
    return ptov(d),ptov(d-res),c[1,1],c[3,1],res

def misalignmentTerm(N,align,app=75.):
    """
    Trace identical set of rays through system with and
    without misalignment of cylindrical optic.
    Return interpolated OPD difference
    """
    #Set up list of rays
    rays = traceToTestOptic(10000,app=app)
    rays2 = copy.deepcopy(rays)

    #Trace each through the system
    perfectCyl(rays)
    backToWFS(rays)
    perfectCyl(rays2,align=align)
    backToWFS(rays2)

    #Return both OPD
    opd1 = anal.interpolateVec(rays,0,100,100,method='linear')[0]
    opd2 = anal.interpolateVec(rays2,0,100,100,method='linear')[0]
    d = (opd1-opd2)/2.
    
    return d - np.nanmean(d)

def cghmisalign1m(N,cghalign):
    rays,l = traceToTestOptic1m(N,coloffset=200.)
    rays2,l = traceToTestOptic1m(N,coloffset=200.,cghalign=cghalign)

    perfectCyl1m(rays)
    perfectCyl1m(rays2)

    backToWFS1m(rays)
    backToWFS1m(rays2,cghalign=cghalign)

    pdb.set_trace()

    opd1,dx,dy = anal.interpolateVec(rays,0,200,200)
    opd2,dx,dy = anal.interpolateVec(rays2,0,200,200)
    d = (opd1-opd2)/2.

    return d - np.nanmean(d)


def backToWFS220(rays):
    """
    Trace rays from nominal test optic tangent plane back to WFS plane.
    This function can also be used with a point source to determine the
    Optimal focus positions of the field lenses.
    +z points toward CGH.
    """
    #Back to CGH
    tran.transform(rays,0,0,220+line,0,0,0)
    surf.flat(rays,nr=1.)
    #Trace back through CGH
    tran.transform(rays,0,0,0,0,1.*pi/180,0)
    tran.transform(rays,0,0,0,1.*pi/180,0,0)
    surf.flat(rays,nr=1.)
    surf.zernphase(rays,cghcoeff,80.,632.82e-6)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    tran.itransform(rays,0,0,0,1.*pi/180,0,0)
    #Go to collimator
    tran.transform(rays,0,0,100,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays,reverse=True)
    #Go to focus
    tran.transform(rays,0,0,1934.99719-100.,0,0,0)
    surf.flat(rays,nr=1.)
    #Place to AC-508-250
    lenses.AC508_250(rays,reverse=True)
    #Go to WFS location
##    tran.transform(rays,0,0,foc,0,0,0)
##    surf.flat(rays,nr=.1)
    tran.transform(rays,0,0,foc,0,0,0)
    surf.flat(rays,nr=1.)

    #Go to cylindrical field lens
    tran.transform(rays,0,0,-cylz,0,0,0)
    surf.flat(rays,nr=1.)
    tran.transform(rays,0,0,0,0,0,pi/2)
    lenses.LJ1516_L2(rays,reverse=False)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,-cylz,0,0,0)
    #Back to WFS
    surf.flat(rays,nr=1.)
    
    return anal.rmsY(rays)

def perfectCyl(rays,align=np.zeros(6)):
    """
    Trace rays from perfect cylinder with potential misalignment
    Assume rays are traced to tangent plane of nominal optic position
    +z points back toward CGH
    Leave with reference frame at tangent plane of nominal surface
    """
    #Apply misalignment
    tran.transform(rays,*align)
    #Trace cylinder
    tran.transform(rays,0,0,220.,0,0,0)
    #Get cylindrical axis in +x direction
    tran.transform(rays,0,0,0,0,0,pi/2)
    surf.cyl(rays,220.,nr=1.)
    tran.reflect(rays)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,220.,0,0,0)
    #Go back to nominal tangent plane
    tran.itransform(rays,*align)
    surf.flat(rays,nr=1.)
    
    return

def traceToTestOptic220(N,app=75.):
    """Trace a set of rays from the point source to the nominal
    test optic location
    Return the rays at the plane tangent to the nominal source position.
    """
    #Set up source
    div = app/1935.033
    rays = sources.pointsource(div,N)
    #Trace through collimator
    tran.transform(rays,0,0,1935.033,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays)
    #Trace to CGH
    tran.transform(rays,0,0,100.,0,0,0)
    #Apply proper CGH misalignment
    pdb.set_trace()
    tran.transform(rays,0,0,0,-1.*pi/180,0,0)
    #Trace through CGH
    surf.flat(rays,nr=1.)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    surf.zernphase(rays,cghcoeff,80.,632.82e-6)
    #Reverse CGH misalignment
    tran.itransform(rays,0,0,0,-1.*pi/180,0,0)
    #Go to line focus
    tran.transform(rays,0,0,0,0,1.*pi/180,0)
    surf.flat(rays,nr=1.)
    tran.transform(rays,0,0,line,0,0,0)
    surf.flat(rays,nr=1.)
    #Go to test optic
    tran.transform(rays,0,0,220.,0,0,0)
    surf.flat(rays,nr=1.)
    #Rotate reference frame so rays impinge toward -z
    tran.transform(rays,0,0,0,0,pi,0)
    
    return rays

def traceToTestOptic1m(N,app=75.,coloffset=0.,cghalign=np.zeros(6)):
    """Trace a set of rays from the point source to the nominal
    test optic location
    Return the rays at the plane tangent to the nominal source position.
    """
    #Set up source
    div = app/1935.033
    rays = sources.pointsource(div,N)
    #Trace through collimator
    tran.transform(rays,0,0,1935.033+coloffset,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays)
##    tran.transform(rays,0,0,-coloffset,0,0,0)
    #Trace to CGH
    tran.transform(rays,0,0,100.,0,0,0)
    #Apply proper CGH misalignment
    tran.transform(rays,0,0,0,-10.*pi/180,0,0)
    #Apply CGH misalignment
    tran.transform(rays,*cghalign)
    #Trace through CGH
    surf.flat(rays,nr=1.)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    surf.zernphase(rays,-cgh1m,80.,632.82e-6)
    #Reverse CGH misalignment
    tran.itransform(rays,*cghalign)
    #Go to line focus
    line = surf.focusY(rays,nr=1.)
    #Go to test optic
    tran.transform(rays,0,0,1000.,0,0,0)
    surf.flat(rays,nr=1.)
    #Go to 1m cylindrical radius of curvature
    px,py = anal.measurePower(rays,200,200)
    tran.transform(rays,0,0,1000+py,0,0,0)
    surf.flat(rays,nr=1.)
    
    return rays,line

def perfectCyl1m(rays,align=np.zeros(6)):
    """
    Trace rays from perfect cylinder with potential misalignment
    Assume rays are traced to tangent plane of nominal optic position
    +z points back toward CGH
    Leave with reference frame at tangent plane of nominal surface
    """
    #Rotate reference frame so rays impinge toward -z
    tran.transform(rays,0,0,0,0,pi,0)
    #Apply misalignment
    tran.transform(rays,*align)
    #Trace cylinder
    tran.transform(rays,0,0,1000.,0,0,0)
    #Get cylindrical axis in +x direction
    tran.transform(rays,0,0,0,0,0,pi/2)
    surf.cyl(rays,1000.,nr=1.)
    tran.reflect(rays)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,1000.,0,0,0)
    #Go back to nominal tangent plane
    tran.itransform(rays,*align)
    surf.flat(rays,nr=1.)
    
    return

def backToWFS1m(rays,cghalign=np.zeros(6)):
    """
    Trace rays from nominal test optic tangent plane back to WFS plane.
    This function can also be used with a point source to determine the
    Optimal focus positions of the field lenses.
    +z points toward CGH.
    """
    #Reverse x,z misalignments
    for i in [0,2,3,5]:
        cghalign[i] = -cghalign[i]
    #Back to CGH
    tran.transform(rays,0,0,1000+line1m,0,0,0)
    surf.flat(rays,nr=1.)
    #Trace back through CGH
    tran.transform(rays,*cghalign)
    surf.zernphase(rays,-cgh1m,80.,632.82e-6)
    tran.refract(rays,1.,nsil)
    tran.transform(rays,0,0,6.35,0,0,0)
    surf.flat(rays,nr=nsil)
    tran.refract(rays,nsil,1.)
    tran.itransform(rays,*cghalign)
    tran.transform(rays,0,0,0,-10.*pi/180,0,0)
    #Go to collimator
    tran.transform(rays,0,0,100,0,0,0)
    surf.flat(rays,nr=1.)
    lenses.collimator6(rays,reverse=True)
    #Go to focus
    tran.transform(rays,0,0,1934.90059-100.,0,0,0)
    surf.flat(rays,nr=1.)
    #Place to AC-508-250
    lenses.AC508_250(rays,reverse=True)
    #Go to WFS location
##    tran.transform(rays,0,0,foc,0,0,0)
##    surf.flat(rays,nr=.1)
    tran.transform(rays,0,0,foc1m,0,0,0)

    #Go to cylindrical field lens
    tran.transform(rays,0,0,-cylz1m,0,0,0)
    surf.flat(rays,nr=1.)
    tran.transform(rays,0,0,0,0,0,pi/2)
    lenses.LJ1144_L2(rays,reverse=False)
    tran.itransform(rays,0,0,0,0,0,pi/2)
    tran.itransform(rays,0,0,-cylz1m,0,0,0)
    #Back to WFS
    surf.flat(rays,nr=1.)
    
    return anal.rmsY(rays)

##def traceToTestOpticP(N):
##    """Trace a set of rays from the point source to the nominal
##    test optic location
##    Return the rays at the plane tangent to the nominal source position.
##    Assumes all paraxial lenses
##    """
##    #Set up source
##    rays = sources.pointsource(.038996,N)
##    #Trace through collimator
##    tran.transform(rays,0,0,1935.033,0,0,0)
##    surf.flat(rays,nr=1.)
##    surf.paraxial(rays,1935.033)
##    #Trace to CGH
##    tran.transform(rays,0,0,100.,0,0,0)
##    #Trace through CGH
##    surf.flat(rays,nr=1.)
##    surf.paraxialY(rays,200.)
##    #Go to test optic
##    tran.transform(rays,0,0,420.,0,0,0)
##    surf.flat(rays,nr=1.)
##    #Rotate reference frame so rays impinge toward -z
##    tran.transform(rays,0,0,0,0,pi,0)
##    
##    return rays
##
##def backToWFSP(rays,cylp):
##    """
##    Trace rays from nominal test optic tangent plane back to WFS plane.
##    This function can also be used with a point source to determine the
##    Optimal focus positions of the field lenses.
##    +z points toward CGH.
##    Assumes paraxial lenses
##    """
##    #Back to CGH
##    tran.transform(rays,0,0,420,0,0,0)
##    surf.flat(rays,nr=1.)
##    surf.paraxialY(rays,200.)
##    #Go to collimator
##    tran.transform(rays,0,0,100,0,0,0)
##    surf.flat(rays,nr=1.)
##    surf.paraxial(rays,1935.033)
##    #Go to focus
##    tran.transform(rays,0,0,1934.99719-100.,0,0,0)
##    surf.flat(rays,nr=1.)
##    #Place to AC-508-250
##    surf.paraxial(rays,250.)
##    #Go to WFS location
####    tran.transform(rays,0,0,foc,0,0,0)
####    surf.flat(rays,nr=.1)
##    #Determine focus location
####    return surf.focusX(rays)
##
##    
##    tran.transform(rays,0,0,271.98618296467038-132.32323232323233,0,0,0)
##    surf.flat(rays,nr=1.)
##    surf.paraxialY(rays,1000.)
##    tran.transform(rays,0,0,132.32323232323233,0,0,0)
##    surf.flat(rays,nr=1.)
##
####    #Go to cylindrical field lens
####    tran.transform(rays,0,0,-cylz,0,0,0)
####    surf.flat(rays,nr=1.)
####    tran.transform(rays,0,0,0,0,0,pi/2)
####    lenses.LJ1516_L2(rays,reverse=False)
####    tran.itransform(rays,0,0,0,0,0,pi/2)
####    tran.itransform(rays,0,0,-cylz,0,0,0)
####    #Back to WFS
####    surf.flat(rays,nr=1.)
##    
##    return anal.rmsY(rays)
##
##def testCollimatorP():
##    """Find out the true focal length of the collimator for
##    positioning of the source in raytrace.
##    """
##    #Create circular plane wave
##    rays = sources.circularbeam(25.4*3,100000)
##    #Trace through collimator
##    lenses.collimator6(rays,reverse=True)
##    #Find focus
##    foc = anal.analyticImagePlane(rays)
##    #Trace to focus and return ray spread and focal distance
##    tran.transform(rays,0,0,foc,0,0,0)
##    surf.flat(rays,nr=1.)
##    return rays
