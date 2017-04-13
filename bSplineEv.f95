!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------BIVARIATE B-SPLINE EVALUATION ROUTINE-----------------------
!-----------------------------------------------------------------------------------
!--------------------------------Benjamin D. Donovan--------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------------
!BIVARIATEBSPLINEEV
!
!A subroutine to return the value of the interpolation at the coordinate (x,y).
!
!Parameters:
!-----------
!x : The x-coordinate of the point of interest.
!y : The y-coordinate of the point of interest.
!u : The knot vector for the x variable.
!      Comes from scipy.interpolate.SmoothBivariateSpline().get_knots()
!v : The knot vector for the y variable.
!      Comes from scipy.interpolate.SmoothBivariateSpline().get_knots()
!coeff : The 4x4 array of the coefficients.
!      Comes from scipy.interpolate.SmoothBivariateSpline().get_coeffs()
!dx : The derivative with respect to x. 
!dy : The derivative with respect to y.
!summ : The z-value at the coordinate (x,y).
!
!Notes:
!------
!1. For interfacing with scipy.interpolate.SmoothBivariateSpline, the array of
!   coefficients must be reshaped from a 1D array with 16 elements to a 2D with
!   4 rows and 4 columns.
!2. If you make the value of 'dx' or 'dy' anything other than 0 or 1, your value
!   from the subroutine will be zero.
!3. Information on the summation performed in this subroutine can be found at:
!   http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/surface/bspline-construct.html
!-----------------------------------------------------------------------------------

subroutine bivariateBSplineEv(x, y, u, v, coeff, dx, dy, summ)

  !Declaration of all variables used in the subroutine.
  implicit none
  real, intent(in) :: x, y
  real, dimension(0:7) :: u, v
  integer, intent(in) :: dx, dy
  real, dimension(0:3,0:3) :: coeff
  integer :: i, j
  real :: bSplineBasisDerivative, bSplineBasisFcn
  real, intent(out) :: summ
  
  !Define the summ value to be zero.
  summ = 0.
  
  !Perform the summation. Information about the summation can be found in the link above.
  do i=0,3
    do j=0,3
      if ((dx==0) .and. (dy==0)) then
        summ = summ + bSplineBasisFcn(i, 3, u, x) * bSplineBasisFcn(j, 3, v, y) * coeff(i, j)
      else if ((dx==1) .and. (dy==0)) then
        summ = summ + bSplineBasisDerivative(i, 3, u, x) * bSplineBasisFcn(j, 3, v, y) * coeff(i, j)
      else if ((dx==0) .and. (dy==1)) then
        summ = summ + bSplineBasisFcn(i, 3, u, x) * bSplineBasisDerivative(j, 3, v, y) * coeff(i, j)
      end if
    end do
  end do

end subroutine bivariateBSplineEv

!-----------------------------------------------------------------------------------
!BSPLINEBASISDERIVATIVE
!
!A recursive function to evaluate the derivative of the i-th B-spline basis function
!of degree p.
!
!Parameters:
!-----------
!i : The i-th B-spline basis function to evaluate.
!p : The degree of the i-th B-spline basis function.
!u : The array of knots.
!x : Point at which to evalute the i-th B-spline basis function derivative of
!    of degree p.
!
!Results:
!--------
!bSBFD : Value of the specified B-spline basis function derivative at x.
!-----------------------------------------------------------------------------------

recursive function bSplineBasisDerivative(i, p, u, x) result(bSBFD)

  !Declaration of variables used in the recursive function.
  implicit none
  integer :: i, p
  real, intent(in) :: x
  real, dimension(0:7) :: u
  real :: c, d, e, f, bSBFD1, bSBFD2, bSBFD, bSplineBasisFcn

  !If the degree of the B-spline basis function equals zero, the derivative of the
  !B-spline basis function is zero.
  if (p .eq. 0) then
    bSBFD = 0.

  !Otherwise, the value needs to be calculated. 
  else 
    if (u(i+p) - u(i) .eq. 0.) then
      c = 0.
      d = 0.
    else 
      c = 1/(u(i+p) - u(i))
      d = (x - u(i))/(u(i+p)-u(i))
    end if

    if (u(i+p+1) - u(i+1) .eq. 0.) then
      e = 0.
      f = 0.
    else
      e = 1/(u(i+p+1) - u(i+1))
      f = (u(i+p+1) - x)/(u(i+p+1) - u(i+1))
    end if

    !Recursive calculation of the B-spline basis function derivative.
    bSBFD1 = c * bSplineBasisFcn(i,p-1,u,x) + d * bSplineBasisDerivative(i,p-1,u,x)
    bSBFD2 = -e * bSplineBasisFcn(i+1,p-1,u,x) + f * bSplineBasisDerivative(i+1,p-1,u,x)
    bSBFD = bSBFD1 + bSBFD2
  end if

end function bSplineBasisDerivative

!-----------------------------------------------------------------------------------
!BSPLINEBASISFCN
!
!A recursive function to evaluate i-th B-spline basis function of degree p.
!
!Parameters:
!-----------
!i : The i-th B-spline basis function to evaluate.
!p : The degree of the i-th B-spline basis function.
!u : The array of knots.
!z : Point at which to evalute the i-th B-spline basis function derivative of
!    of degree p.
!
!Results:
!--------
!bSBF : Value of the specified B-spline basis function derivative at x.
!-----------------------------------------------------------------------------------
recursive function bSplineBasisFcn(i, p, u, z) result(bSBF)

  !Declaration of variables used in the recursive function. 
  implicit none
  integer :: i, p
  real, intent(in) :: z 
  real, dimension(0:7) :: u
  real :: a, b, bSBF

  !If the degree of the basis function is zero, the basis function is zero or one,
  !depending on whether or not z is between the i-th and the i-th + 1 knot.
  if (p .eq. 0) then
    if (u(i) .le. z .and. z .lt. u(i+1)) then
      bSBF = 1.
    else
      bSBF = 0.
    end if

  !If the degree of the basis function is not zero, it must be calculated.
  else
    if (u(i+p) - u(i) .eq. 0.) then
      a = 0
    else
      a = (z - u(i))/(u(i+p)-u(i))
    end if

  if (u(i+p+1) - u(i+1) .eq. 0.) then
    b = 0
  else 
    b = (u(i+p+1) - z)/(u(i+p+1) - u(i+1))
  end if

  !Recursive calculation of the B-spline basis function.
  bSBF = a * bSplineBasisFcn(i, p-1, u, z) + b * bSplineBasisFcn(i+1, p-1, u, z)
  end if

end function bSplineBasisFcn