!This function rotates a vector in a right handed fashion
!axis is 1,2,3 for x,y,z axis rotation
subroutine rotatevector(x,y,z,theta,axis)
  !Declarations
  real*8, intent(inout) :: x,y,z,theta
  integer, intent(in) :: axis
  real*8 :: output(3)
  !real*8, dimension(3) :: rotatevector

  if (axis==1) then
    output(1) = x
    output(2) = cos(theta)*y-sin(theta)*z
    output(3) = sin(theta)*y+cos(theta)*z
  else if (axis==2) then
    output(1) = cos(theta)*x+sin(theta)*z
    output(2) = y
    output(3) = -sin(theta)*x+cos(theta)*z
  else
    output(1) = cos(theta)*x-sin(theta)*y
    output(2) = sin(theta)*x+cos(theta)*y
    output(3) = z
  end if

  x = output(1)
  y = output(2)
  z = output(3)

end subroutine rotatevector

!This function rotates a vector in a right handed fashion
!axis is given by input ux,uy,uz
subroutine rotateaxis(x,y,z,theta,ux,uy,uz)
  !Declarations
  real*8, intent(inout) :: x,y,z
  real*8, intent(inout) :: theta,ux,uy,uz
  real*8 :: output(3),s,c,mag
  !real*8, dimension(3) :: rotatevector

  !Ensure axis is normalized
  mag = sqrt(ux**2 + uy**2 + uz**2)
  ux = ux/mag
  uy = uy/mag
  uz = uz/mag

  c = cos(theta)
  s = sin(theta)
  output(1) = (c+ux**2*(1-c))*x + (ux*uy*(1-c)-uz*s)*y + (ux*uz*(1-c)+uy*s)*z
  output(2) = (uy*ux*(1-c)+uz*s)*x + (c+uy**2*(1-c))*y + (uy*uz*(1-c)-ux*s)*z
  output(3) = (uz*ux*(1-c)-uy*s)*x + (uz*uy*(1-c)+ux*s)*y + (c+uz**2*(1-c))*z

  x = output(1)
  y = output(2)
  z = output(3)

end subroutine rotateaxis

!Function to reflect about local surface normal
!Assume ray cosines point into surface, surface normal points out of surface
!If i = -incident, then reflected ray is 2(i.n)n-i
subroutine reflect(l,m,n,ux,uy,uz,num)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8 :: dot
  integer :: i

  !Loop through rays and reflect about normal
  !$omp parallel do private(i,dot)
  do i=1,num
    !Compute dot of incident with normal
    dot = ux(i)*l(i) + uy(i)*m(i) + uz(i)*n(i)
    !Compute reflected direction
    l(i) = l(i) - 2*dot*ux(i)
    m(i) = m(i) - 2*dot*uy(i)
    n(i) = n(i) - 2*dot*uz(i)
  end do
  !$omp end parallel do

end subroutine reflect

!Refract from one index to another
subroutine refract(l,m,n,ux,uy,uz,num,n1,n2)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8, intent(in) :: n1,n2
  real*8, intent(inout) :: l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i
  real*8 :: dot,t1,t2,alpha,cx,cy,cz

  !Loop through rays
  !$omp parallel do private(dot,t1,t2,cx,cy,cz,alpha)
  do i=1,num
    !Ensure normal vector is pointing into second index (dot product should be positive)
    dot = l(i)*ux(i) + m(i)*uy(i) + n(i)*uz(i)
    if (dot < 0) then
      ux(i) = -ux(i)
      uy(i) = -uy(i)
      uz(i) = -uz(i)
      dot = -dot
    end if
    !If wavevector is equal to surface normal, do nothing
    !print *, l(i),m(i),n(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dummy
    !if ((l(i) .eq. ux(i)) .and. (m(i) .eq. uy(i)) .and. (n(i) .eq. uz(i))) then
    if (dot==1) then
      cycle
    end if
    !print *, "Passed conditional"
    !read *, dummy
    !Compute Snell's law
    t1 = acos(dot)
    t2 = asin((n1/n2)*sin(t1))
    !Compute cross product
    cx = uy(i)*n(i)-m(i)*uz(i)
    cy = l(i)*uz(i)-ux(i)*n(i)
    cz = ux(i)*m(i)-l(i)*uy(i)
    !Rotate about cross vector an angle t2-t1
    call rotateaxis(l(i),m(i),n(i),t2-t1,cx,cy,cz)
    !Normalize
    alpha = sqrt(l(i)**2 + m(i)**2 + n(i)**2)
    l(i) = l(i)/alpha
    m(i) = m(i)/alpha
    n(i) = n(i)/alpha
  end do
  !$omp end parallel do

end subroutine refract

!Coordinate system transform, translations are done first, then rotations in XYZ order
!Rotations act to rotate a surface via the right hand rule
subroutine transform(x,y,z,l,m,n,ux,uy,uz,num,tx,ty,tz,rx,ry,rz)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: tx,ty,tz,rx,ry,rz
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i

  !Loop through rays
  !$omp parallel do
  do i=1,num
    !Perform translation
    x(i) = x(i) + tx
    y(i) = y(i) + ty
    z(i) = z(i) + tz
    !Perform x rotation
    call rotatevector(x(i),y(i),z(i),rx,1)
    call rotatevector(l(i),m(i),n(i),rx,1)
    call rotatevector(ux(i),uy(i),uz(i),rx,1)
    !Perform y rotation
    call rotatevector(x(i),y(i),z(i),ry,2)
    call rotatevector(l(i),m(i),n(i),ry,2)
    call rotatevector(ux(i),uy(i),uz(i),ry,2)
    !Perform z rotation
    call rotatevector(x(i),y(i),z(i),rz,3)
    call rotatevector(l(i),m(i),n(i),rz,3)
    call rotatevector(ux(i),uy(i),uz(i),rz,3)
  end do
  !$omp end parallel do

end subroutine transform

!Coordinate system transform, rotations in ZYX order, then translations
!Rotations act to rotate a surface via the right hand rule
!This is meant to invert the transform routine
subroutine itransform(x,y,z,l,m,n,ux,uy,uz,num,tx,ty,tz,rx,ry,rz)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(inout) :: tx,ty,tz,rx,ry,rz
  real*8, intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  integer :: i

  !Loop through rays
  !$omp parallel do
  do i=1,num
    !Perform z rotation
    call rotatevector(x(i),y(i),z(i),-rz,3)
    call rotatevector(l(i),m(i),n(i),-rz,3)
    call rotatevector(ux(i),uy(i),uz(i),-rz,3)
    !Perform y rotation
    call rotatevector(x(i),y(i),z(i),-ry,2)
    call rotatevector(l(i),m(i),n(i),-ry,2)
    call rotatevector(ux(i),uy(i),uz(i),-ry,2)
    !Perform x rotation
    call rotatevector(x(i),y(i),z(i),-rx,1)
    call rotatevector(l(i),m(i),n(i),-rx,1)
    call rotatevector(ux(i),uy(i),uz(i),-rx,1)
    !Perform translation
    x(i) = x(i) - tx
    y(i) = y(i) - ty
    z(i) = z(i) - tz
  end do
  !$omp end parallel do

end subroutine itransform

!Radially grooved grating diffraction
!Assumes grating in x y plane, with grooves converging at
!hubdist in positive y direction
subroutine radgrat(x,y,l,m,n,wave,num,dpermm,order)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num)
  real*8, intent(in) :: dpermm,wave,order
  integer :: i
  real*8 :: d, yaw, pi, dum, det, sn

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Save sign of n
    sn = n(i) / abs(n(i))
    !Compute local d spacing in nm
    d = dpermm * sqrt(y(i)**2 + x(i)**2)
    !Compute local yaw
    yaw = pi/2 - atan(-x(i)/abs(y(i)))
    !print *, x(i),y(i),d,yaw
    !print *, l(i),m(i),n(i)

    !Evanescence?
    !det = l(i)**2+m(i)**2
    !Compute new direction cosines - evanescence will result in NaNs
    l(i) = l(i) + sin(yaw)*order*wave/d
    m(i) = m(i) - cos(yaw)*order*wave/d
    n(i) = sn*sqrt(1. - l(i)**2 - m(i)**2)

  end do

end subroutine radgrat

!Radially grooved grating diffraction
!Assumes grating in x y plane, with grooves converging at 
!hubdist in positive y direction
!Placed with origin at center of grating
subroutine radgratcenter(x,y,l,m,n,wave,num,dpermm,order,hubdist)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num)
  real*8, intent(in) :: dpermm,wave,order,hubdist
  integer :: i
  real*8 :: d, yaw, pi, dum, det, sn

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Save sign of n
    sn = n(i) / abs(n(i))
    !Compute local d spacing in nm
    d = dpermm * sqrt((-y(i)+hubdist)**2 + x(i)**2)
    !Compute local yaw
    yaw = pi/2 - atan(x(i)/abs((hubdist-y(i))))
    !print *, x(i),y(i),d,yaw
    !print *, l(i),m(i),n(i)

    !Evanescence?
    !det = l(i)**2+m(i)**2
    !Compute new direction cosines - evanescence will result in NaNs
    l(i) = l(i) + sin(yaw)*order*wave/d
    m(i) = m(i) - cos(yaw)*order*wave/d
    n(i) = sn*sqrt(1. - l(i)**2 - m(i)**2)
    
  end do

end subroutine radgratcenter

!Radially grooved grating diffraction with wavelength vector
!Assumes grating in x y plane, with grooves converging at
!hubdist in positive y direction
subroutine radgratW(x,y,l,m,n,wave,num,dpermm,order)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num),wave(num)
  real*8, intent(inout) :: l(num),m(num),n(num)
  real*8, intent(in) :: dpermm,order
  integer :: i
  real*8 :: d, yaw, pi, dum, det, sn

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  do i=1,num
    !Compute local d spacing in nm
    d = dpermm * sqrt(y(i)**2 + x(i)**2)
    !Compute local yaw
    sn = n(i) / abs(n(i))
    yaw = pi/2 - atan(x(i)/abs(y(i)))
    !print *, x(i),y(i),d,yaw
    !print *, l(i),m(i),n(i)

    !Evanescence?
    !det = l(i)**2+m(i)**2
    !Compute new direction cosines - evanescence will result in NaNs
    l(i) = l(i) + sin(yaw)*order*wave(i)/d
    m(i) = m(i) - cos(yaw)*order*wave(i)/d
    n(i) = sn*sqrt(1. - l(i)**2 - m(i)**2)
  end do

end subroutine radgratW

!Linear grating with groove period d
!Wavelength wave
!Groove direction assumed in y direction
subroutine grat(x,y,l,m,n,num,d,order,wave)
  !Declarations
  integer, intent(in) :: num
  real*8, intent(in) :: x(num),y(num)
  real*8, intent(inout) :: l(num),m(num),n(num),wave(num),order(num)
  real*8, intent(in) :: d
  integer :: i
  real*8 :: pi, dum

  pi = acos(-1.)

  !Loop through rays, compute new diffracted ray direction
  !$omp parallel do
  do i=1,num
    !Save sign of n
    sn = n(i) / abs(n(i))
    !Compute new direction cosines
    l(i) = l(i) - order(i)*wave(i)/d
    n(i) = sn*sqrt(1 - l(i)**2 - m(i)**2)
    !Evanescence?
    if ((l(i)**2+m(i)**2)>1) then
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    end if
  end do
  !$omp end parallel do

end subroutine grat
