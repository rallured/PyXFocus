import numpy as np

###### SOURCES #######

def pointsource(ang,num):
    """Define point source with angular divergence
    Points in +z direction
    Ang is half angle
    """
    #Radial direction cosine magnitude
    rho = np.sqrt(np.random.rand(num))*np.sin(ang)
    theta = np.random.rand(num)*2*np.pi
    l = rho*np.cos(theta)
    m = rho*np.sin(theta)
    n = np.sqrt(1.-l**2-m**2)
    x = np.repeat(0.,num)
    y = np.repeat(0.,num)
    z = np.repeat(0.,num)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def circularbeam(rad,num):
    """Define uniform, circular beam of radius rad, pointing in +z direction
    """
    rho = np.sqrt(np.random.rand(num))*rad
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def annulus(rin,rout,num):
    """Define annulus of rays pointing in +z direction
    """
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(1.,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def subannulus(rin,rout,dphi,num,zhat=1.):
    """Create a subapertured annulus source in +z direction
    Annulus is centered about theta=0 which points to +x
    If negz is set -1, rays will point in -z hat
    """
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*dphi - dphi/2.
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0.,num)
    l = np.repeat(0.,num)
    m = np.repeat(0.,num)
    n = np.repeat(zhat,num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def rectArray(xsize,ysize,num):
    """Creates a regular array of rays using meshgrid and linspace"""
    x,y = np.meshgrid(np.linspace(-xsize,xsize,num),\
                      np.linspace(-ysize,ysize,num))
    opd = np.repeat(0.,num**2)
    x = x.flatten()
    y = y.flatten()
    z = np.copy(opd)
    l = np.repeat(0.,num**2)
    m = np.repeat(0.,num**2)
    n = np.repeat(1.,num**2)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def convergingbeam(zset,rin,rout,tmin,tmax,num,lscat):
    """Converging sub-apertured annulus beam
    Place at nominal focus
    Input z position, inner and outer radius,
    min and max theta
    """
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = tmin + np.random.rand(num)*(tmax-tmin)
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(zset,num)
    lscat = lscat * np.tan((np.random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -np.cos(np.arctan(rho/zset)+lscat)
    l = -np.sqrt(1-n**2)*np.cos(theta)
    m = -np.sqrt(1-n**2)*np.sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def convergingbeam2(zset,xmin,xmax,ymin,ymax,num,lscat):
    """Rectangular converging beam
    Place at nominal focus
    Input z position and rectangular bounds
    """
    x = xmin + np.random.rand(num)*(xmax-xmin)
    y = ymin + np.random.rand(num)*(ymax-ymin)
    rho = np.sqrt(x**2+y**2)
    theta = np.arctan2(y,x)
    z = np.repeat(zset,num)
    lscat = lscat * np.tan((np.random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -np.cos(np.arctan(rho/zset)+lscat)
    l = -np.sqrt(1-n**2)*np.cos(theta)
    m = -np.sqrt(1-n**2)*np.sin(theta)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return [opd,x,y,z,l,m,n,ux,uy,uz]

def rectbeam(xhalfwidth,yhalfwidth,num):
    """Rectangular beam pointing in +z direction
    """
    x = (np.random.rand(num)-.5)*2*xhalfwidth
    y = (np.random.rand(num)-.5)*2*yhalfwidth
    z = np.repeat(0.,num)
    n = np.repeat(1.,num)
    l = np.copy(z)
    m = np.copy(z)
    ux = np.repeat(0.,num)
    uy = np.repeat(0.,num)
    uz = np.repeat(0.,num)
    opd = np.repeat(0.,num)
    return [opd,x,y,z,l,m,n,ux,uy,uz]
