
import traces.PyTrace as PT
from numpy import *
from matplotlib.pyplot import *
import pdb

def traceBeam(gratAlign,order=1.,wave=2.48,blaze=30.):
    """Traces a beam to grating, applies misalignment,
    then traces to focal plane and returns mean
    x,y position
    x refers to spectral direction
    """
    #Compute yaw to achieve Littrow configuration
    phi0 = arccos(sin(1.5*pi/180)/cos(blaze*pi/180))
    yaw = -tan(sin(blaze*pi/180)/tan(phi0))
##    yaw = yaw - pi/2
##    pdb.set_trace()
    
    #Set up beam
    rays = PT.pointsource(0.,10) #multiple rays to avoid Python complaints

    #Rotate to nominal grating reference frame
    PT.transform(rays,0,0,0,pi/2+1.5*pi/180,0,0)

    #Apply misalignment
    PT.transform(rays,*gratAlign)

    #Trace grating
    PT.flat(rays)
    PT.transform(rays,0,0,0,0,0,yaw) #-1.73211678*pi/180
    PT.radgrat(rays,11500.,160./11500.,order,wave)
    PT.itransform(rays,0,0,0,0,0,yaw)

    #Reverse misalignment
    PT.itransform(rays,*gratAlign)

    #Trace to focal plane
    PT.transform(rays,0,11500.,0,0,0,0)
    PT.transform(rays,0,0,0,pi/2,0,0)
    PT.flat(rays)

    #Return mean positions
    return mean(rays[1]),mean(rays[2])

def dofTest(dof,step,order=1.,wave=2.48,blaze=30):
    #Get nominal beam focal plane location
    misalign = zeros(6)
    tr = lambda misalign: traceBeam(misalign,order=order,wave=wave,blaze=blaze)
    xn,yn = tr(misalign)

    #Loop until spectral spot shift goes beyond 6 microns
    while True:
        misalign[dof] = misalign[dof] + step
        x,y = tr(misalign)
        eff = y-yn
        spec = x-xn
        if isnan(spec):
            print 'Evanescence cutoff at %.2e' % misalign[dof]
            return misalign[dof],0
        if abs(eff) > 1.:
            print 'Effective area cutoff at %.2e' % misalign[dof]
            return misalign[dof],1
        if abs(spec) > 6.e-3:
            print 'Spectral cutoff at %.2e' % misalign[dof]
            return misalign[dof],2

    return misalign[dof]

def scanDoFWave(dof,step,order=1.,waves=linspace(.9,5.1,100)):
    cutoff = zeros(size(waves))
    mode = copy(cutoff)

    for w in waves:
        cutoff[waves==w],mode[waves==w] = \
            dofTest(dof,step,order=order[waves==w],wave=w)

    return cutoff,mode.astype('int')

def alignPlot():
    #Loop through DoFs and make sensitivity plots
    wave = linspace(.9,5.1,100)
    order = repeat(3.,size(wave))
    order[wave>2.55] = 1.
    order[logical_and(wave>1.275,wave<2.55)] = 2.
    
    ion()
    figure()
    pfmts = ['-.','--','-']
    col = ['b','g','r']
    trans = ['Lateral','Axial','Radial']
    ang = ['Pitch','Roll','Yaw']
    for dof in range(3):
        cutoff,mode = scanDoFWave(dof,.05,order=order,waves=wave)
        print trans[dof] + ' Cutoff: ' + str(min(cutoff))
        plot(wave,cutoff,col[dof]+pfmts[mode[0]],label=trans[dof])
        cutoff,mode = scanDoFWave(dof,-.05,order=order,waves=wave)
        plot(wave,cutoff,col[dof]+pfmts[mode[0]])
    figure()
    for dof in range(3,6):
        cutoff,mode = scanDoFWave(dof,.5/60**2*pi/180,order=order,waves=wave)
        print ang[dof-3] + ' Cutoff: ' + str(min(cutoff))
        plot(wave,cutoff*180/pi*60**2,col[dof-3]+pfmts[mode[0]],label=ang[dof-3])
        cutoff,mode = scanDoFWave(dof,-.5/60**2*pi/180,order=order,waves=wave)
        plot(wave,cutoff*180/pi*60**2,col[dof-3]+pfmts[mode[0]])
