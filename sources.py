import numpy as np


def pointsource(ang, num):
    '''
    Point source with a specified angular divergence.

    Note: Rays points in the +z direction.

    Parameters
    ----------
    ang : float
        Angular divergence of rays.
    num : int
        Number of rays to create.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    # Radial direction cosine magnitude
    rho = np.sqrt(np.random.rand(num))*np.sin(ang)
    theta = np.random.rand(num)*2*np.pi
    l = rho*np.cos(theta)
    m = rho*np.sin(theta)
    n = np.sqrt(1.-l**2-m**2)
    x = np.repeat(0., num)
    y = np.repeat(0., num)
    z = np.repeat(0., num)
    ux = np.repeat(0., num)
    uy = np.repeat(0., num)
    uz = np.repeat(0., num)
    opd = np.repeat(0., num)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def circularbeam(rad, num):
    '''
    Uniform, circular beam with specified radius.

    Note: Rays point in +z direction.

    Parameters
    ----------
    rad : int / float
        Radius of circular beam.
    num : int
        Number of rays to create.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    rho = np.sqrt(np.random.rand(num))*rad
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0., num)
    l = np.repeat(0., num)
    m = np.repeat(0., num)
    n = np.repeat(1., num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def annulus(rin, rout, num, zhat=-1.):
    '''
    Annulus of rays with specified inner and outer radii.

    Note: Default has rays pointing in -z direction.

    Parameters
    ----------
    rin : int / float
        Inner radius of annulus.
    rout : int / float
        Outer radius of annulus.
    num : int
        Number of rays to create.
    zhat : float
        Direction in which rays point. Default is zhat = -1.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*2*np.pi
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0., num)
    l = np.repeat(0., num)
    m = np.repeat(0., num)
    n = np.repeat(zhat, num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def subannulus(rin, rout, dphi, num, zhat=1.):
    '''
    Subapertured annulus source with specified inner and outer
    radii, as well as angular extent.

    Note: Annulus is centered about theta = 0, which points
          in the +x direction.

    Parameters
    ----------
    rin : int / float
        Inner radius of annulus.
    rout : int / float
        Outer radius of annulus.
    dphi : int / float
        Full angular width of subannulus.
    num : int
        Number of rays to create.
    zhat : float
        Direction in which rays point. Default is zhat = +1.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = np.random.rand(num)*dphi - dphi/2.
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(0., num)
    l = np.repeat(0., num)
    m = np.repeat(0., num)
    n = np.repeat(zhat, num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def xslit(xin, xout, num, zhat=-1.):
    '''
    Slit of rays with specified width in the x-dimension.

    Note: Rays are linearly spaced in the x-dimension.

    Parameters
    ----------
    xin : int / float
        Inner x-position of slit.
    xout : int / float
        Outer x-position of slit.
    num : int
        Number of rays to create.
    zhat : float
        Direction in which rays point. Default is zhat = -1.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    x = np.linspace(xin, xout, num)
    y = np.repeat(0., num)
    z = np.repeat(0., num)
    l = np.repeat(0., num)
    m = np.repeat(0., num)
    n = np.repeat(zhat, num)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)
    opd = np.copy(l)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def rectArray(xsize, ysize, num):
    '''
    Rectangular beam with specified width (x) and height (y) on the
    x-y plane.

    Created with numpy.meshgrid and np.linspace.

    Note: Rays point in the +z direction.

    Parameters
    ----------
    xsize : int / float
        1/2 of rectangle size in the x-dimension.
    ysize : int / float
        1/2 of rectangle size in the y-dimension.
    num : int
        Number of rays to create.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    x, y = np.meshgrid(np.linspace(-xsize, xsize, num),
                       np.linspace(-ysize, ysize, num))
    opd = np.repeat(0., num**2)
    x = x.flatten()
    y = y.flatten()
    z = np.copy(opd)
    l = np.repeat(0., num**2)
    m = np.repeat(0., num**2)
    n = np.repeat(1., num**2)
    ux = np.copy(l)
    uy = np.copy(l)
    uz = np.copy(l)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def convergingbeam(zset, rin, rout, tmin, tmax, num, lscat):
    '''
    Converging sub-apertured annulus beam with specified inner
    and outer radii, as well as angular extent.

    Note: Place rays at nominal focus. Input z position, inner
          and outer radii, and minimum and maximum theta.

    Parameters
    ----------
    zset : int / float
        Z-dimension which defines the convergence of beam.
    rin : int / float
        Inner radius of sub-apertured annulus beam.
    rout : int / float
        Outer radius of sub-apertured annulus beam.
    tmin : int / float
        Minimum angular extent of sub-apertured annulus beam.
    tmax : int/ float
        Maximum angular extent of sub-apertured annulus beam.
    num : int
        Number of rays to create.
    lscat : int / float
        Scatter in the angular convergence [arcsec].

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    rho = np.sqrt(rin**2+np.random.rand(num)*(rout**2-rin**2))
    theta = tmin + np.random.rand(num)*(tmax-tmin)
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    z = np.repeat(zset, num)
    lscat = lscat * np.tan((np.random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -np.cos(np.arctan(rho/zset)+lscat)
    l = -np.sqrt(1-n**2)*np.cos(theta)
    m = -np.sqrt(1-n**2)*np.sin(theta)
    ux = np.repeat(0., num)
    uy = np.repeat(0., num)
    uz = np.repeat(0., num)
    opd = np.repeat(0., num)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def convergingbeam2(zset, xmin, xmax, ymin, ymax, num, lscat):
    '''
    Converging rectangular beam with specified extent in both
    the x- and y-dimensions.

    Note: Place rays at nominal focus. Input z position and
          rectangular bounds.

    Parameters
    ----------
    zset : int / float
        Z-dimension which defines the convergence of beam.
    xmin : int / float
        Minimum extent in the x-dimension.
    xmax : int / float
        Maximum extent in the x-dimension.
    ymin : int / float
        Minimum extent in the y-dimension.
    ymax : int / float
        Maximum extent in the y-dimension.
    num : int
        Number of rays to create.
    lscat : int / float
        Scatter in the angular convergence [arcsec].

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    x = xmin + np.random.rand(num)*(xmax-xmin)
    y = ymin + np.random.rand(num)*(ymax-ymin)
    rho = np.sqrt(x**2+y**2)
    theta = np.arctan2(y, x)
    z = np.repeat(zset, num)
    lscat = lscat * np.tan((np.random.rand(num) - .5)*np.pi)
    lscat = lscat/60**2 * np.pi/180.
    n = -np.cos(np.arctan(rho/zset)+lscat)
    l = -np.sqrt(1-n**2)*np.cos(theta)
    m = -np.sqrt(1-n**2)*np.sin(theta)
    ux = np.repeat(0., num)
    uy = np.repeat(0., num)
    uz = np.repeat(0., num)
    opd = np.repeat(0., num)

    return [opd, x, y, z, l, m, n, ux, uy, uz]


def rectbeam(xhalfwidth, yhalfwidth, num):
    '''
    Rectangular beam with specified width (x) and height (y)
    on the x-y plane, pointing in +z direction.

    Parameters
    ----------
    xhalfwidth : int / float
        1/2 of rectangle extent in the x-dimension.
    yhalfwidth : int / float
        1/2 of rectangle extent in the y-dimension.
    num : int
        Number of rays to create.

    Returns
    -------
    rays : list
        List of ray parameters (opd, x, y, z, l, m, n, ux, uy, uz).

    '''
    x = (np.random.rand(num)-.5)*2*xhalfwidth
    y = (np.random.rand(num)-.5)*2*yhalfwidth
    z = np.repeat(0., num)
    n = np.repeat(1., num)
    l = np.copy(z)
    m = np.copy(z)
    ux = np.repeat(0., num)
    uy = np.repeat(0., num)
    uz = np.repeat(0., num)
    opd = np.repeat(0., num)

    return [opd, x, y, z, l, m, n, ux, uy, uz]
