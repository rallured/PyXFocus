include 'specialFunctions.f95'

!Function to trace set of rays to Zernike surface
!Inputs are position and cosines of rays
!Array of Zernike coefficients
!Output is position of rays at surface
!and new cosines after reflection
subroutine tracezern(x,y,z,l,m,n,ux,uy,uz,num,coeff,rorder,aorder,arrsize,rad)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize, num
  real*8, intent(in) :: coeff(arrsize),rad
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer, intent(in) :: rorder(arrsize),aorder(arrsize)
  integer :: i,c
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,rho,theta,delta,t
  real*8 :: zern(arrsize),rhoder(arrsize),thetader(arrsize)
  real*8 :: zernike,zernrhoder,zernthetader
  real*8 :: dum

  !Loop through each individual ray, trace to surface, and reflect
  !Establish convergence criteria (surface function should be < 1.e-12)
  !$omp parallel do private(t,delta,rho,theta,F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,i,c,zern,rhoder,thetader)
  do i=1,num
    t = 0.
    delta = 100.
    !count = 0
    !print *, x(i), y(i), z(i)
    !print *, l(i), m(i), n(i)
    !read *, dum
    do while (abs(delta) > 1.e-10)
      !count = count + 1
      !Convert to cylindrical coordinates
      rho = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      !Compute surface function and derivatives
      F = z(i)
      Frho = 0.
      Ftheta = 0.
      call zernset(rho/rad,theta,rorder,aorder,arrsize,zern,rhoder,thetader)
      !print *, zern
      !read *, dum
      do c=1,arrsize
        !F = F - coeff(c)*zernike(rho/rad,theta,rorder(c),aorder(c))
        !Frho = Frho - coeff(c)*zernrhoder(rho/rad,theta,rorder(c),aorder(c))/rad
        !Ftheta = Ftheta - coeff(c)*zernthetader(rho/rad,theta,rorder(c),aorder(c))
        F = F - coeff(c)*zern(c)
        Frho = Frho - coeff(c)*rhoder(c)/rad
        Ftheta = Ftheta - coeff(c)*thetader(c)
        !print *, rho/rad, " ", theta, " ", rorder(c), " ", aorder(c)
        !print *, zern(c), rhoder(c)
        !print *, coeff(c)
        !print *, F
        !read *, dum
      end do
      !print *, rho/rad
      !print *, theta
      !print *, zernrhoder(rho/rad,theta,rorder(c),aorder(c))
      !print *, zernthetader(rho/rad,theta,rorder(c),aorder(c))
      !Convert back to cartesian coordinates
      !print *, Frho, Ftheta
      !read *, dum
      Frhox = (x(i)/rho) * Frho
      Frhoy = (y(i)/rho) * Frho
      Fthetax = (-y(i)/rho) * Ftheta/rho
      Fthetay = (x(i)/rho) * Ftheta/rho
      !if (rho==0) then
      !  Frhox = Frho
      !  Frhoy = 0.
      !  Fthetax = 0.
      !  Fthetay = 0.
      !end if
      Fx = Frhox + Fthetax
      Fy = Frhoy + Fthetay
      Fz = 1.
      !print *, Fx, Fy, Fz
      !print *, l(i), m(i), n(i)
      !read *, dum
      !Compute delta and new ray position
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      !print *, 'Fp:' , Fp
      !read *, dum
      delta = -F/Fp
      !print *, F, Fp, delta
      !read *, dum
      x(i) = x(i) + l(i)*delta
      y(i) = y(i) + m(i)*delta
      z(i) = z(i) + n(i)*delta
      !Keep track of total length ray is traced
      t = t + delta
    end do
    !print *, "Number of iterations: ", count
    !We have converged, do any vignetting and compute surface normal
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
  end do
  !$omp end parallel do

end subroutine tracezern

!Function to trace set of rays to Zernike surface
!Inputs are position and cosines of rays
!Array of Zernike coefficients
!Output is position of rays at surface
!and new cosines after reflection
subroutine tracezernOPD(opd,x,y,z,l,m,n,ux,uy,uz,num,coeff,rorder,aorder,arrsize,rad,nr)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize, num
  real*8, intent(in) :: coeff(arrsize),rad,nr
  real*8, intent(inout) :: opd(num),x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer, intent(in) :: rorder(arrsize),aorder(arrsize)
  integer :: i,c
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,rho,theta,delta,t
  real*8 :: zern(arrsize),rhoder(arrsize),thetader(arrsize)
  real*8 :: zernike,zernrhoder,zernthetader
  real*8 :: dum

  !Loop through each individual ray, trace to surface, and reflect
  !Establish convergence criteria (surface function should be < 1.e-12)
  !$omp parallel do private(t,delta,rho,theta,F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,i,c,zern,rhoder,thetader)
  do i=1,num
    t = 0.
    delta = 100.
    !count = 0
    !print *, x(i), y(i), z(i)
    !print *, l(i), m(i), n(i)
    !read *, dum
    do while (abs(delta) > 1.e-10)
      !count = count + 1
      !Convert to cylindrical coordinates
      rho = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      !Compute surface function and derivatives
      F = z(i)
      Frho = 0.
      Ftheta = 0.
      call zernset(rho/rad,theta,rorder,aorder,arrsize,zern,rhoder,thetader)
      !print *, zern
      !read *, dum
      do c=1,arrsize
        !F = F - coeff(c)*zernike(rho/rad,theta,rorder(c),aorder(c))
        !Frho = Frho - coeff(c)*zernrhoder(rho/rad,theta,rorder(c),aorder(c))/rad
        !Ftheta = Ftheta - coeff(c)*zernthetader(rho/rad,theta,rorder(c),aorder(c))
        F = F - coeff(c)*zern(c)
        Frho = Frho - coeff(c)*rhoder(c)/rad
        Ftheta = Ftheta - coeff(c)*thetader(c)
        !print *, rho/rad, " ", theta, " ", rorder(c), " ", aorder(c)
        !print *, zern(c), rhoder(c)
        !print *, coeff(c)
        !print *, F
        !read *, dum
      end do
      !print *, rho/rad
      !print *, theta
      !print *, zernrhoder(rho/rad,theta,rorder(c),aorder(c))
      !print *, zernthetader(rho/rad,theta,rorder(c),aorder(c))
      !Convert back to cartesian coordinates
      !print *, Frho, Ftheta
      !read *, dum
      Frhox = (x(i)/rho) * Frho
      Frhoy = (y(i)/rho) * Frho
      Fthetax = (-y(i)/rho) * Ftheta/rho
      Fthetay = (x(i)/rho) * Ftheta/rho
      !if (rho==0) then
      !  Frhox = Frho
      !  Frhoy = 0.
      !  Fthetax = 0.
      !  Fthetay = 0.
      !end if
      Fx = Frhox + Fthetax
      Fy = Frhoy + Fthetay
      Fz = 1.
      !print *, Fx, Fy, Fz
      !print *, l(i), m(i), n(i)
      !read *, dum
      !Compute delta and new ray position
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      !print *, 'Fp:' , Fp
      !read *, dum
      delta = -F/Fp
      !print *, F, Fp, delta
      !read *, dum
      x(i) = x(i) + l(i)*delta
      y(i) = y(i) + m(i)*delta
      z(i) = z(i) + n(i)*delta
      !Keep track of total length ray is traced
      t = t + delta
    end do
    !print *, "Number of iterations: ", count
    !We have converged, do any vignetting and compute surface normal
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Iterate OPD
    opd(i) = opd(i) + t*nr
  end do
  !$omp end parallel do

end subroutine tracezernOPD

!Trace a Zernike standard phase surface as in ZEMAX
subroutine zernphase(opd,x,y,z,l,m,n,ux,uy,uz,num,coeff,rorder,aorder,arrsize,rad,wave)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize, num
  real*8, intent(in) :: coeff(arrsize),rad,wave
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num),opd(num)
  integer, intent(in) :: rorder(arrsize),aorder(arrsize)
  real*8 :: zern(arrsize),rhoder(arrsize),thetader(arrsize)
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy
  real*8 :: pi,theta,rho,dum
  integer :: i,c

  pi = acos(-1.)

  do i=1,num
    !Compute the set of Zernike polynomials
    rho = sqrt(x(i)**2+y(i)**2)
    theta = atan2(y(i),x(i))
    call zernset(rho/rad,theta,rorder,aorder,arrsize,zern,rhoder,thetader)
    !Sum derivatives
    F = 0.
    Frho = 0.
    Ftheta = 0.
    do c=1,arrsize
      F = F + coeff(c)*zern(c)
      Frho = Frho + coeff(c)*rhoder(c)/rad
      Ftheta = Ftheta + coeff(c)*thetader(c)
    end do
    !Convert to Cartesian coordinates
    Frhox = (x(i)/rho) * Frho
    Frhoy = (y(i)/rho) * Frho
    Fthetax = (-y(i)/rho) * Ftheta/rho
    Fthetay = (x(i)/rho) * Ftheta/rho
    Fx = Frhox + Fthetax
    Fy = Frhoy + Fthetay
    !print *, rho,theta,Frho, Ftheta,Fx,Fy,wave
    !read *, dum
    !Add appropriate perturbations to ray wavevector
    l(i) = l(i) + Fx*wave!/2/pi
    m(i) = m(i) + Fy*wave!/2/pi
    n(i) = sign(sqrt(1.-l(i)**2-m(i)**2),n(i))
    opd(i) = opd(i) + F*wave
  end do

end subroutine zernphase

!Function to trace set of rays to Zernike surface
!Inputs are position and cosines of rays
!Array of Zernike coefficients
!Output is position of rays at surface
!and new cosines after reflection
subroutine tracezernrot(x,y,z,l,m,n,ux,uy,uz,num,coeff1,rorder1,aorder1,arrsize1,coeff2,rorder2,aorder2,arrsize2,rad,rot)
  !Declarations
  implicit none
  integer, intent(in) :: arrsize1, arrsize2, num
  real*8, intent(in) :: coeff1(arrsize1),coeff2(arrsize2),rad,rot
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer, intent(in) :: rorder1(arrsize1),aorder1(arrsize1)
  integer, intent(in) :: rorder2(arrsize2),aorder2(arrsize2)
  integer :: i,c
  real*8 :: F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,rho,theta,delta,t
  real*8 :: zern1(arrsize1),rhoder1(arrsize1),thetader1(arrsize1)
  real*8 :: zern2(arrsize2),rhoder2(arrsize2),thetader2(arrsize2)
  real*8 :: zernike,zernrhoder,zernthetader
  real*8 :: dum

  !Loop through each individual ray, trace to surface, and reflect
  !Establish convergence criteria (surface function should be < 1.e-12)
  !$omp parallel do private(t,delta,rho,theta,F,Frho,Ftheta,Frhox,Frhoy,Fthetax,Fthetay,Fx,Fy,Fz,Fp,i,c,zern1,zern2)
  do i=1,num
    t = 0.
    delta = 100.
    !count = 0
    !print *, x(i), y(i), z(i)
    !print *, l(i), m(i), n(i)
    !read *, dum
    do while (abs(delta) > 1.e-10)
      !count = count + 1
      !Convert to cylindrical coordinates
      rho = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      !Compute surface function and derivatives
      F = z(i)
      Frho = 0.
      Ftheta = 0.
      call zernset(rho/rad,theta,rorder1,aorder1,arrsize1,zern1,rhoder1,thetader1)
      !print *, zern1(1)
      !read *, dum
      do c=1,arrsize1
        !F = F - coeff(c)*zernike(rho/rad,theta,rorder(c),aorder(c))
        !Frho = Frho - coeff(c)*zernrhoder(rho/rad,theta,rorder(c),aorder(c))/rad
        !Ftheta = Ftheta - coeff(c)*zernthetader(rho/rad,theta,rorder(c),aorder(c))
        F = F - coeff1(c)*zern1(c)
        Frho = Frho - coeff1(c)*rhoder1(c)/rad
        Ftheta = Ftheta - coeff1(c)*thetader1(c)
        !print *, rho/rad, " ", theta, " ", rorder(c), " ", aorder(c)
        !print *, zern(c), rhoder(c)
        !print *, coeff(c)
        !print *, F
        !read *, dum
      end do
      call zernset(rho/rad,theta+rot,rorder2,aorder2,arrsize2,zern2,rhoder2,thetader2)
      !print *, zern2(1)
      do c=1,arrsize2
        F = F - coeff2(c)*zern2(c)
        Frho = Frho - coeff2(c)*rhoder2(c)/rad
        Ftheta = Ftheta - coeff2(c)*thetader2(c)
      end do
      !print *, rho/rad
      !print *, theta
      !print *, zernrhoder(rho/rad,theta,rorder(c),aorder(c))
      !print *, zernthetader(rho/rad,theta,rorder(c),aorder(c))
      !Convert back to cartesian coordinates
      !print *, Frho, Ftheta
      !read *, dum
      Frhox = (x(i)/rho) * Frho
      Frhoy = (y(i)/rho) * Frho
      Fthetax = (-y(i)/rho) * Ftheta/rho
      Fthetay = (x(i)/rho) * Ftheta/rho
      !if (rho==0) then
      !  Frhox = Frho
      !  Frhoy = 0.
      !  Fthetax = 0.
      !  Fthetay = 0.
      !end if
      Fx = Frhox + Fthetax
      Fy = Frhoy + Fthetay
      Fz = 1.
      !print *, Fx, Fy, Fz
      !print *, l(i), m(i), n(i)
      !read *, dum
      !Compute delta and new ray position
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      !print *, 'Fp:' , Fp
      !read *, dum
      delta = -F/Fp
      !print *, F, Fp, delta
      !read *, dum
      x(i) = x(i) + l(i)*delta
      y(i) = y(i) + m(i)*delta
      z(i) = z(i) + n(i)*delta
      !Keep track of total length ray is traced
      t = t + delta
    end do
    !print *, "Number of iterations: ", count
    !We have converged, do any vignetting and compute surface normal
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
  end do
  !$omp end parallel do

end subroutine tracezernrot
