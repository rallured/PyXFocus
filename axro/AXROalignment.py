from numpy import *
from matplotlib.pyplot import *
import traces.conicsolve as conicsolve
import traces.PyTrace as PT
import pdb
from mpl_toolkits.mplot3d import Axes3D
import time
import scipy.optimize

#Load in flat mirror deformations
foldfig = genfromtxt("/home/rallured/Dropbox/AXRO/Alignment/Simulation/"
                        "NIST/141202FoldFigCoeffs.txt")
foldsag = genfromtxt("/home/rallured/Dropbox/AXRO/Alignment/Simulation/"
                        "141202FoldSagCoeffs.txt")*1000
foldcoeffs = foldfig + foldsag
retrofig = genfromtxt("/home/rallured/Dropbox/AXRO/Alignment/Simulation/"
                      "NIST/141202RetroFigCoeffs.txt")
retrosag = genfromtxt("/home/rallured/Dropbox/AXRO/Alignment/Simulation/"
                      "141202RetroSagCoeffs.txt")*1000
retrocoeffs = retrosag + retrosag

#Load in primary deformations
pcoeff,pax,paz = genfromtxt('/home/rallured/Dropbox/AXRO/'
                        'Alignment/CoarseAlignment/150615_OP1S09Coeffs.txt')
pcoeff = pcoeff/1000.

#Set up Hartmann mask
holewidth = arcsin(.3/220)
holetheta = linspace(-arcsin(45./220),arcsin(45./220),15)
numholes = size(holetheta)

#Set up diverging beam with angular offset
#Give it divergence angle, pitch, and roll
def CDAbeam(num,div,pitch,roll,cda):
##    PT.transform(*cda)
    PT.pointsource(div,num)
    PT.transform(0,0,0,pitch,0,0)
    PT.transform(0,0,0,0,0,roll)
    PT.itransform(*cda)
    return

#Trace from primary focus to fold to primary
def primMaskTrace(fold,primary,woltVignette=True,foldrot=0.):
    #Get Wolter parameters
    alpha,p,d,e = conicsolve.woltparam(220.,8400.)
    primfoc = conicsolve.primfocus(220.,8400.)

    #Trace to fold mirror
    #translate to center of fold mirror
    PT.transform(0.,85.12,primfoc-651.57+85.12,0,0,0)
    #rotate so surface normal points in correct direction
    PT.transform(0,0,0,-3*pi/4,0,0)
    PT.transform(0,0,0,0,0,pi)
    #trace to fold flat
    PT.flat()
    #Introduce fold misalignment
    PT.transform(*fold)    
    PT.zernsurfrot(foldsag,foldfig,406./2,-174.659*pi/180+foldrot)
    PT.itransform(*fold)
    PT.reflect()
    PT.transform(0,0,0,0,0,-pi)
    PT.transform(0,0,0,pi/4,0,0)

    #Translate to optical axis mid-plane, then down to image of
    #primary focus, place primary mirror and trace
    PT.transform(0,85.12,651.57-85.12,0,0,0)
    PT.flat()
##    pdb.set_trace()
    rt = conicsolve.primrad(8475.,220.,8400.)
    PT.transform(0,-rt,75.,0,0,0)
    PT.transform(*primary)
    PT.transform(0,rt,-8475.,0,0,0)
##    PT.wolterprimary(220.,8400.)
    PT.primaryLL(220.,8400.,8525.,8425.,30.*np.pi/180.,pcoeff,pax,paz)
    if woltVignette is True:
        ind = logical_and(PT.z<8525.,PT.z>8425.)
        PT.vignette(ind=ind)
    PT.reflect()
    PT.transform(0,-rt,8475.,0,0,0)
    PT.itransform(*primary)
    PT.transform(0,rt,-8475.,0,0,0)

    #Move back up to mask plane and trace flat
    PT.transform(0,0,8400.+134.18,0,0,0)
    PT.flat()
##    pdb.set_trace()

    #Rays should now be at Hartmann mask plane

    return

def traceFromMask(N,numholes,cda,fold,retro,primary,foldrot=0.,retrorot=0.):
    #Vignette at proper hole
    h = hartmannMask()
    ind = h==N
    PT.vignette(ind=ind)

    #Continue trace up to retro and back to CDA
    PT.transform(0,-123.41,1156.48-651.57-134.18,0,0,0)
    PT.flat()
    PT.transform(0,0,0,pi,0,0)
    PT.transform(*retro)
    PT.zernsurfrot(retrosag,retrofig,378./2,-8.993*pi/180+retrorot)
    PT.itransform(*retro)
    PT.reflect()
    PT.transform(0,0,0,-pi,0,0)
    PT.transform(0,123.41,-1156.48+651.57+134.18,0,0,0)
    PT.flat()
    h = hartmannMask()
    ind = h==N
    PT.vignette(ind=ind)
    PT.transform(0,0,-134.18,0,0,0)
    rt = conicsolve.primrad(8475.,220.,8400.)
    PT.transform(0,-rt,75.,0,0,0)
    PT.transform(*primary)
    PT.transform(0,rt,-8475.,0,0,0)
    PT.wolterprimary(220.,8400.)
    ind = logical_and(PT.z<8525.,PT.z>8425.)
    PT.vignette(ind=ind)
    PT.reflect()
    PT.transform(0,-rt,8475.,0,0,0)
    PT.itransform(*primary)
    PT.transform(0,rt,-8475.,0,0,0)
    PT.transform(0,-85.12,8400.-651.57+85.12\
                 ,0,0,0)
    PT.transform(0,0,0,-pi/4,0,0)
    PT.transform(0,0,0,0,0,pi)
    PT.flat()
    PT.transform(*fold)
    PT.zernsurfrot(foldsag,foldfig,406./2,-174.659*pi/180+foldrot)
    PT.itransform(*fold)
    PT.reflect()
    PT.transform(0,0,0,0,0,-pi)
    PT.transform(0,0,0,3*pi/4,0,0)
    PT.transform(0,-85.12,-85.12-(conicsolve.primfocus(220.,8400.)-651.57)\
                 ,0,0,0)
    PT.transform(*cda)
    PT.flat()

    return


#### DOUBLE MIRROR TRACES ####
#Trace from focus to fold to Hartmann mask
def fullMaskTrace(fold,prim,sec,woltVignette=True,foldrot=0.):
    #Get Wolter parameters
    alpha,p,d,e = conicsolve.woltparam(220.,8400.)
    foc = 8400.

    #Trace to fold mirror
    #translate to center of fold mirror
    PT.transform(0.,85.12,foc-651.57+85.12,0,0,0)
    #rotate so surface normal points in correct direction
    PT.transform(0,0,0,-3*pi/4,0,0)
    PT.transform(0,0,0,0,0,pi)
    #trace to fold flat
    PT.flat()
    #Introduce fold misalignment
    PT.transform(*fold)    
    PT.zernsurfrot(foldsag,foldfig,406./2,-174.659*pi/180+foldrot)
    PT.itransform(*fold)
    PT.reflect()
    PT.transform(0,0,0,0,0,-pi)
    PT.transform(0,0,0,pi/4,0,0)

    #Translate to optical axis mid-plane, then down to image of
    #primary focus, place primary mirror and trace
    PT.transform(0,85.12,651.57-85.12,0,0,0)
    PT.flat()
    PT.transform(0,0,-8400.,0,0,0)

    #Place secondary
    #Go to tangent point, apply misalignment, place mirror, and reverse
    PT.transform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.transform(*sec)
    PT.itransform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.woltersecondary(220.,8400.)
    if woltVignette is True:
        ind = logical_and(PT.z<8375.,PT.z>8275.)
        PT.vignette(ind=ind)
    PT.reflect()
    PT.transform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.itransform(*sec) #Back at nominal secondary tangent point
    PT.itransform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    

    #Place primary
    #Go to tangent point, apply misalignment, place mirror, and reverse
    PT.transform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
    PT.transform(*prim)
    PT.itransform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
##    PT.transform(0,0,8475.,0,0,0)
##    PT.flat()
##    PT.itransform(0,0,8475.,0,0,0)
##    PT.wolterprimary(220.,8400.)
    PT.primaryLL(220.,8400.,8525.,8425.,30.*np.pi/180.,pcoeff,pax,paz)
##    pdb.set_trace()
    if woltVignette is True:
        ind = logical_and(PT.z<8525.,PT.z>8425.)
        PT.vignette(ind=ind)
    PT.reflect()
    PT.transform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
    PT.itransform(*prim)
    PT.itransform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)

    #Move back up to mask plane and trace flat
    PT.transform(0,0,8400.+134.18,0,0,0)
    PT.flat()
##    pdb.set_trace()

    #Rays should now be at Hartmann mask plane

    return

def fullFromMask(N,cda,fold,retro,prim,sec,foldrot=0.,retrorot=0.):
##    pdb.set_trace()
    #Vignette at proper hole
    h = hartmannMask()
    ind = h==N
    PT.vignette(ind=ind)

    #Continue trace up to retro and back to CDA
    PT.transform(0,-123.41,1156.48-651.57-134.18,0,0,0)
    PT.flat()
    PT.transform(0,0,0,pi,0,0)
    PT.transform(*retro)
    PT.zernsurfrot(retrosag,retrofig,378./2,-8.993*pi/180+retrorot)
    PT.itransform(*retro)
    PT.reflect()
    PT.transform(0,0,0,-pi,0,0)

    #Back to mask
    PT.transform(0,123.41,-1156.48+651.57+134.18,0,0,0)
    PT.flat()
    h = hartmannMask()
    ind = h==N
    PT.vignette(ind=ind)

    #Place Wolter surfaces
    PT.transform(0,0,-134.18-8400.,0,0,0)
    PT.transform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
    PT.transform(*prim)
    PT.itransform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
##    PT.wolterprimary(220.,8400.)
    PT.primaryLL(220.,8400.,8525.,8425.,30.*np.pi/180.,pcoeff,pax,paz)
    pdb.set_trace()
    ind = logical_and(PT.z<8525.,PT.z>8425.)
    PT.vignette(ind=ind)
    PT.transform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
    PT.itransform(*prim)
    PT.itransform(0,-conicsolve.primrad(8425.,220.,8400.),8425.,0,0,0)
    PT.reflect()

    #Wolter secondary
    PT.transform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.transform(*sec)
    PT.itransform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.woltersecondary(220.,8400.)
    ind = logical_and(PT.z<8375.,PT.z>8275.)
    PT.vignette(ind=ind)
    PT.reflect()
    PT.transform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)
    PT.itransform(*sec)
    PT.itransform(0,-conicsolve.secrad(8325.,220.,8400.),8325.,0,0,0)


    
##    PT.woltersecondary(220.,8400.)
##    ind = logical_and(PT.z<8375.,PT.z>8275.)
##    PT.vignette(ind=ind)
##    PT.reflect()

    #Back to fold
    PT.transform(0,-85.12,8400.-651.57+85.12\
                 ,0,0,0)
    PT.transform(0,0,0,-pi/4,0,0)
    PT.transform(0,0,0,0,0,pi)
    PT.flat()
    PT.transform(*fold)
    PT.zernsurfrot(foldsag,foldfig,406./2,-174.659*pi/180+foldrot)
    PT.itransform(*fold)
    PT.reflect()
    PT.transform(0,0,0,0,0,-pi)

    #Back to CDA
    PT.transform(0,0,0,3*pi/4,0,0)
    PT.transform(0,-85.12,-85.12-8400.+651.57\
                 ,0,0,0)
    PT.transform(*cda)
    PT.flat()

    return

#Return nominal pitch and roll for Hartmann hole
#Used as starting point for optimization
def hartmannStartFull(N):
    global holetheta
    a,p,d,e = conicsolve.woltparam(220.,8400.) #Pitch is 4*a

    return 4*a,holetheta[N-1]

#Return mean x and y positions as a function of pitch and roll
def traceHoleFull(num,N,div,p,r,cda,fold,prim,sec,foldrot=0.):
    CDAbeam(num,div,p,r,cda)
    fullMaskTrace(fold,prim,sec,woltVignette=False,foldrot=foldrot)
    realx,realy = hartmannPosition(N)
    res = sqrt(mean((PT.x-realx)**2+(PT.y-realy)**2))
    print str(N) + ': ' + str(res)
    return res

#Use minimization routine to aim ray bundle at proper hole
def fullAim(num,N,div,cda,fold,prim,sec,foldrot=0.):
    #Create function
    fun = lambda p: traceHoleFull(num,N,div,p[0],p[1],\
                               cda,fold,prim,sec,foldrot=foldrot)
    #Optimize function
    start = array(hartmannStartFull(N))
    if abs(start[1]) < .001:
        start[1] = .01
    print 'Begin ' + str(N)
    res = scipy.optimize.minimize(fun,start,method='nelder-mead',\
                  options={'ftol':1.e-2,'disp':True})
    print 'End ' +str(N)
##    traceHole2(num,N,numholes,div,res['x'][0],res['x'][1],cda,fold)

    return res['x']

##Return vector of pitch and roll for a Hartmann mask
def alignHartmannFull(cda,fold,prim,sec,foldrot=0.):
    global numholes
    pitch = zeros(numholes)
    roll = zeros(numholes)

    for i in range(numholes):
        pitch[i],roll[i] = fullAim(10**2,i+1,.00001*pi/180,\
                cda,fold,prim,sec,foldrot=0.)#findHole(i+1,numholes,cda,fold)

    return pitch,roll

##Do full double mirror alignment trace, making use of ray aiming results
##for efficiency
def fullAlign(num,cda,fold,retro,prim,sec,p=None,r=None,\
                foldrot=0.,retrorot=0.):
    global numholes
    if p is None:
        #Grab pitch and roll vectors
        p,r = alignHartmannFull(cda,fold,prim,sec,foldrot=foldrot)

    #Trace out holes, one by one
    xm = []
    ym = []
    xstd = []
    ystd = []
    for i in range(numholes):
        #Trace up to Hartmann mask
        CDAbeam(num,.001*pi/180,p[i],r[i],cda)
        fullMaskTrace(fold,prim,sec,woltVignette=True)
        fullFromMask(i+1,cda,fold,retro,prim,sec)

        #Print out vignetting factor
        print size(PT.x)/num

        #Evaluate mean spot position
        xm.append(mean(PT.x))
        ym.append(mean(PT.y))
        xstd.append(std(PT.x))
        ystd.append(std(PT.y))

    return array(xm),array(ym),array(xstd),array(ystd)

#Evaluate sensitivity of misalignment degree of freedom
def dofSensitivityFull(num,obj,dof,step,criteria):
    #Initial misalignment vectors
    cda = zeros(6)
    fold = zeros(6)
    retro = zeros(6)
    prim = zeros(6)
    sec = zeros(6)
    misalign = [cda,fold,retro,prim,sec]
    
    #Get nominal spot position
    x0,y0 = fullAlign(num,*misalign)

    #Increase proper dof until spot shifts breach 10 micron requirement
    merit = 0.
    figure()
    while merit < criteria:
        misalign[obj][dof] = misalign[obj][dof] + step
        try:
            x1,y1 = fullAlign(num,*misalign)
        except:
            sys.stdout.write('Hartmann Throughput Cutoff at %7.4e' %\
                             misalign[obj][dof])
            break
        dx = x1-x0
        dx = dx - mean(dx)
        dy = y1-y0
        dy = dy - mean(dy)
        merit = max(sqrt(dx**2+dy**2))
        sys.stdout.write('DoF: %7.4e Merit : %0.4f\r' %\
                         (misalign[obj][dof],merit))
        plot(dx,dy,'.')
        draw()
        sys.stdout.flush()
    
    return

#### DOUBLE MIRROR TRACES ####


#Define Hartmann mask and vignette rays
#Keep track of which hole rays hit with "hole" vector
#Need origin to be at Hartmann mask plane on optical axis
def hartmannMask():
##    #Create hole array to handle hole positions
##    hole = zeros(size(PT.x))
##    holerad = (conicsolve.primrad(8500.,220.,8400.)-220.)/2. #Hole halfwidth
##    holecent = conicsolve.primrad(8475.,220.,8400.) #Radius of center of holes
##    halfang = arcsin(50./220.)-.009 #Half angle of Hartmann mask

##    holetheta = linspace(-halfang,halfang,numholes) #Vector of Hartmann angles
##    holewidth = arcsin(holerad/220.)

    #Set holewidth and holetheta to be global variables
    global holewidth, holetheta
    
    #Loop through hole numbers
    rayang = arctan2(PT.y,PT.x) #Center of mirror is -pi/2
    i = 1
    hole = zeros(size(PT.x))
    for theta in holetheta:
        ind = logical_and(rayang < -pi/2 + theta + holewidth,\
                          rayang > -pi/2 + theta - holewidth)
        hole[ind] = i
        i = i+1

    return hole

#Set up CDA beam to trace to a given Hartmann hole
#Will use indicated beam divergence and apply
#appropriate roll to beam to find hole N
def traceHole(num,N,numholes,div,cda):
    a,p,d,e = conicsolve.woltparam(220.,8400.) #Pitch is 2*a
    halfang = arcsin(50./220.)-.009
    holetheta = linspace(-halfang,halfang,numholes)

    CDAbeam(num,div,2*a,holetheta[N-1],cda)

    return 2*a,holetheta[N-1]

#Return nominal pitch and roll for Hartmann hole
#Used as starting point for optimization
def hartmannStart(N,numholes):
    global holetheta
    a,p,d,e = conicsolve.woltparam(220.,8400.) #Pitch is 2*a

    return 2*a,holetheta[N-1]

#Return mean x and y positions as a function of pitch and roll
def traceHole2(num,N,numholes,div,p,r,cda,fold,primary,foldrot=0.):
    CDAbeam(num,div,p,r,cda)
    primMaskTrace(fold,primary,woltVignette=False,foldrot=foldrot)
    realx,realy = hartmannPosition(N)
    return mean(sqrt((PT.x-realx)**2+(PT.y-realy)**2))

#Use minimization routine to aim ray bundle at proper hole
def aimHole(num,N,numholes,div,cda,fold,primary,foldrot=0.):
    #Create function
    fun = lambda p: traceHole2(num,N,numholes,div,p[0],p[1],\
                               cda,fold,primary,foldrot=foldrot)
    #Optimize function
    start = array(hartmannStart(N,numholes))
    if abs(start[1]) < .001:
        start[1] = .01
    res = scipy.optimize.minimize(fun,start,method='nelder-mead',\
                   options={'ftol':1.e-2,'disp':False})
##    traceHole2(num,N,numholes,div,res['x'][0],res['x'][1],cda,fold)

    return res['x']

#Returns nominal position of Hartmann hole in Cartesian coordinates
def hartmannPosition(N):
    global holetheta
    holecent = conicsolve.primrad(8475.,220.,8400.) #Radius of center of holes
    thistheta = holetheta[N-1]-pi/2
    return holecent*cos(thistheta),holecent*sin(thistheta)

#Returns nominal position of Hartmann hole in polar coordinates
def hartmannAngles(N,numholes):
    holecent = conicsolve.primrad(8475.,220.,8400.) #Radius of center of holes
    halfang = arcsin(50./220.)-.009 #Half angle of Hartmann mask
    holetheta = linspace(-halfang,halfang,numholes) #Vector of Hartmann angles
    thistheta = holetheta[N-1]-pi/2
    return holecent,thistheta

#Determine chief ray to a given Hartmann hole
#Start with fairly wide divergence at nominal location
#Take mean of rays that hit the hole, then repeat with much
#smaller divergence
#Probably three iterations to converge
#Need to do this for each Hartmann hole and save ray directions
def findHole(N,numholes,cda,fold,primary):
    #Trace out first iteration
    pitch,roll = traceHole(10**3,N,numholes,.0001*pi/180,cda)
    pdb.set_trace()
    primMaskTrace(fold,woltVignette=False)
    pdb.set_trace()
    #Mean position of rays should be equal to nominal Hartmann position
    rad,roll = hartmannAngles(N,numholes)
    actrad = mean(sqrt(PT.x**2+PT.y**2))
    actroll = mean(arctan2(PT.y,PT.x))
    raddiff = rad - actrad
    rolldiff = roll - actroll
    newpitch = pitch - raddiff/(conicsolve.primfocus(220.,8400.)+134.18)
    newroll = roll - rolldiff + pi/2
    pdb.set_trace()

    while (abs(raddiff) > .01) or (abs(rolldiff) > .01/220.):     
        #Fix pitch and roll
        CDAbeam(10**3,.0001*pi/180,newpitch,newroll,cda)
        pitch = newpitch
        roll = newroll-pi/2
        primMaskTrace(fold,woltVignette=False)
        raddiff = rad - mean(sqrt(PT.x**2+PT.y**2))
        rolldiff = roll - mean(arctan2(PT.y,PT.x))
##        print 'Rolldiff: ' + str(rolldiff)
##        print 'Raddiff: ' + str(raddiff)
        newpitch = pitch - raddiff/(conicsolve.primfocus(220.,8400.)+134.18)
        newroll = roll - rolldiff + pi/2
    
    return newpitch,newroll

##Return vector of pitch and roll for a Hartmann mask
def alignHartmann(numholes,cda,fold,primary,foldrot=0.):
    pitch = zeros(numholes)
    roll = zeros(numholes)

    for i in range(numholes):
        pitch[i],roll[i] = aimHole(10**2,i+1,numholes,.00001*pi/180,\
                                   cda,fold,primary,foldrot=0.)#findHole(i+1,numholes,cda,fold)

    return pitch,roll
        
##Do full primary alignment trace, making use of ray aiming results
##for efficiency
def fullPrimary(num,numholes,cda,fold,retro,primary,p=None,r=None,\
                foldrot=0.,retrorot=0.):
    if p is None:
        #Grab pitch and roll vectors
        p,r = alignHartmann(numholes,cda,fold,primary,foldrot=foldrot)

    #Trace out holes, one by one
    xm = []
    ym = []
    for i in range(numholes):
        #Trace up to Hartmann mask
        CDAbeam(num,.001*pi/180,p[i],r[i],cda)
        primMaskTrace(fold,primary,woltVignette=True)
        traceFromMask(i+1,numholes,cda,fold,retro,primary)

        #Evaluate mean spot position
        xm.append(mean(PT.x))
        ym.append(mean(PT.y))

    return array(xm),array(ym)

#Evaluate sensitivity of misalignment degree of freedom
def dofSensitivity(num,numholes,obj,dof,step):
    #Initial misalignment vectors
    cda = zeros(6)
    fold = zeros(6)
    retro = zeros(6)
    primary = zeros(6)
    misalign = [cda,fold,retro,primary]
    
    #Get nominal spot position
    x0,y0 = fullPrimary(num,numholes,*misalign)

    #Increase proper dof until spot shifts breach 10 micron requirement
    merit = 0.
    figure()
    while merit < .01:
        misalign[obj][dof] = misalign[obj][dof] + step
        try:
            x1,y1 = fullPrimary(num,numholes,*misalign)
        except:
            sys.stdout.write('Hartmann Throughput Cutoff at %7.4f' %\
                             misalign[obj][dof])
            break
        dx = x1-x0
        dx = dx - mean(dx)
        dy = y1-y0
        dy = dy - mean(dy)
        merit = max(sqrt(dx**2+dy**2))
        sys.stdout.write('DoF: %7.4f Merit : %0.4f\r' %\
                         (misalign[obj][dof],merit))
        plot(dx,dy,'.')
        draw()
        sys.stdout.flush()
    
    return

#Evaluate sensitivity of misalignment degree of freedom
def rotSensitivity(num,numholes,step,foldrot=False,retrorot=False):
    #Initial misalignment vectors
    cda = zeros(6)
    fold = zeros(6)
    retro = zeros(6)
    primary = zeros(6)
    misalign = [cda,fold,retro,primary]
    
    #Get nominal spot position
    x0,y0 = fullPrimary(num,numholes,*misalign)

    #Increase proper dof until spot shifts breach 10 micron requirement
    merit = 0.
    figure()
    while merit < .01:
        if foldrot is not False:
            foldrot = foldrot + step
            retrorot = 0.
        else:
            retrorot = retrorot + step
            foldrot = 0.
        try:
            x1,y1 = fullPrimary(num,numholes,*misalign,\
                                foldrot=foldrot,retrorot=retrorot)
        except:
            sys.stdout.write('Hartmann Throughput Cutoff at %7.4f' %\
                             max([foldrot,retrorot])*180/pi)
            break
        dx = x1-x0
        dx = dx - mean(dx)
        dy = y1-y0
        dy = dy - mean(dy)
        merit = max(sqrt(dx**2+dy**2))
        sys.stdout.write('DoF: %7.4f Merit : %0.4f\r' %\
                         (max([foldrot,retrorot])*180/pi,merit))
        plot(dx,dy,'.')
        draw()
        sys.stdout.flush()
    
    return
