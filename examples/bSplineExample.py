#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-------------------------------B-SPLINE EXAMPLE CODE-------------------------------
#-----------------------------------------------------------------------------------
#--------------------------------Benjamin D. Donovan--------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

#Before this example can work, one must compile the bSplineEv.f95 with f2py.
#Change directories to the location of bSplineEv.f95, then perform:
#f2py -c -m bSplineEv bSplineEv.f95

from numpy import *
import bSplineEv as b
import scipy.interpolate


def bSplineExample(xData, yData, zData, x, y):
    '''
    Evaluates the bivariate B-spline interpolation at a coordinate (x,y).
    
    Parameters
    ----------
    xData : 1D array
        The array of x-coordinates to be interpolated.
    yData : 1D array
        The array of y-coordinates to be interpolated.
    zData : 1D array
        The array of z-coordinates to be interpolated.
    x : float
        The x-coordinate at which to evaluate the interpolation.
    y : float
        The y-coordinate at which to evaluate the interpolation.
        
    Returns
    -------
    
    '''
    
    #Peform the interpolation.
    interp = scipy.interpolate.SmoothBivariateSpline(xData, yData, zData)
    
    #Get the knot vectors and the coefficient vector.
    xKnots, yKnots = interp.get_knots()
    coeffs = interp.get_coeffs()
    
    #Reshape the coefficient vector.
    coeffs = reshape(coeffs,(4,4))
    
    #Evaluate the bivariate B-spline interpolation at the coordinate (x,y).
    interpVal = b.bivariatebsplineev(x, y, xKnots, yKnots, coeffs, 0, 0)
    dxVal = b.bivariatebsplineev(x, y, xKnots, yKnots, coeffs, 1, 0)
    dyVal = b.bivariatebsplineev(x, y, xKnots, yKnots, coeffs, 0, 1)
    
    return interpVal, dxVal, dyVal