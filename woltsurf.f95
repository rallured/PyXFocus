include 'specialFunctions.f95'

!This function traces to a Wolter I primary mirror
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine wolterprimary(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,psi
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 2*(1+2*psi)/(1+psi) * alpha
  thetap = 2*psi/(1+psi) * alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  Fz = 2*p
  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fp)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-8)
      F = 2*p*z(i) + p**2 + 4*e**2*p*d/(e**2-1) - x(i)**2 - y(i)**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wolterprimary

!This function traces to a Wolter I primary mirror
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine wolterprimaryopd(opd,x,y,z,l,m,n,ux,uy,uz,num,r0,z0,psi,nr)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: opd(num),x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,psi,nr
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 2*(1+2*psi)/(1+psi) * alpha
  thetap = 2*psi/(1+psi) * alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  Fz = 2*p
  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fp)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      F = 2*p*z(i) + p**2 + 4*e**2*p*d/(e**2-1) - x(i)**2 - y(i)**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      opd(i) = opd(i) + nr*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wolterprimaryopd

!This function traces to a Wolter I secondary mirror
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine woltersecondary(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,psi
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 2*(1+2*psi)/(1+psi) * alpha
  thetap = 2*psi/(1+psi) * alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-8)
      F = e**2*(d+z(i))**2 - z(i)**2 - x(i)**2 - y(i)**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fz = 2*e**2*(d+z(i)) - 2*z(i)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F, Fx, Fy, Fz
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersecondary

!This function traces to a Wolter I primary mirror with sinusoidal perturbation
!Defined by Van Speybroeck prescription
!For WFS test, use flat to get rays close so they find correct intersection
!Surface should be placed at common focus with z+ pointing toward mirrors
subroutine woltersine(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,amp,freq)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,amp,freq
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad
  integer :: i

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      rad = sqrt(x(i)**2+y(i)**2) + amp*sin(2*acos(-1.)*freq*z(i))
      F = 2*p*z(i) + p**2 + 4*e**2*p*d/(e**2-1) - rad**2
      Fx = -2.*x(i)
      Fy = -2.*y(i)
      Fz = 2.*p - 2*rad*amp*2*acos(-1.)*freq*cos(2*acos(-1.)*freq*z(i))
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersine

!Construct a paraboloid as in Wolter but with Legendre-Legendre
!deformations. Define Legendre and Legendre derivative functions.
!Pass in coeff, axial order, and azimuthal order as in Zemax implementation
subroutine wolterprimLL(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,zmax,zmin,dphi,coeff,axial,az,cnum)
  !Declarations
  implicit none
  integer, intent(in) :: num,cnum
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,zmax,zmin,dphi,coeff(cnum)
  integer, intent(in) :: axial(cnum),az(cnum)
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad,add,addx,addy,addz,G
  integer :: i,a
  real*8 :: pi,legendre,legendrep,zarg,targ,ang

  !Compute Van Speybroeck parameters
  pi = acos(-1.)
  alpha = .25*atan(r0/z0)
  thetah = 3.*alpha
  thetap = alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad,ang,zarg,targ,add,addx,addy,addz,a,G)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      ang = atan2(y(i),x(i))
      zarg = (z(i)-((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
      targ = 2*ang/dphi
      
      !Compute Legendre additive terms
      add = 0.
      addx = 0.
      addy = 0.
      addz = 0.
      do a=1,cnum
        add = add + coeff(a)*legendre(zarg,axial(a))*legendre(targ,az(a))
        addx = addx - coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(y(i)/(y(i)**2+x(i)**2))
        addy = addy + coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(x(i)/(y(i)**2+x(i)**2))
        addz = addz + coeff(a)*legendrep(zarg,axial(a))*legendre(targ,az(a))*2/(zmax-zmin)
      end do
      G = sqrt(x(i)**2+y(i)**2) + add
      F = -(G**2 - p**2 - 2*p*z(i) - 4*e**2*p*d/(e**2-1))
      Fx = -2*G*(x(i)/sqrt(x(i)**2+y(i)**2)+addx)
      Fy = -2*G*(y(i)/sqrt(x(i)**2+y(i)**2)+addy)
      Fz = 2*p - 2*G*addz
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Set rays outside mirror definition to NaN
    !if (abs(zarg)>1 .or. abs(targ) > 1) then
    !  x(i) = 0.
    !  y(i) = 0.
    !  z(i) = 0.
    !end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wolterprimLL

!Construct a hyperboloid as in Wolter but with Legendre-Legendre
!deformations. Define Legendre and Legendre derivative functions.
!Pass in coeff, axial order, and azimuthal order as in Zemax implementation
subroutine woltersecLL(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,psi,zmax,zmin,dphi,coeff,axial,az,cnum)
  !Declarations
  implicit none
  integer, intent(in) :: num,cnum
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,psi,zmax,zmin,dphi,coeff(cnum)
  integer, intent(in) :: axial(cnum),az(cnum)
  real*8 :: alpha,thetah,thetap,p,d,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad,add,addx,addy,addz,G
  integer :: i,a
  real*8 :: pi,legendre,legendrep,zarg,targ,ang

  !Compute Van Speybroeck parameters
  alpha = .25*atan(r0/z0)
  thetah = 2*(1+2*psi)/(1+psi) * alpha
  thetap = 2*psi/(1+psi) * alpha
  p = z0*tan(4*alpha)*tan(thetap)
  d = z0*tan(4*alpha)*tan(4*alpha-thetah)
  e = cos(4*alpha)*(1+tan(4*alpha)*tan(thetah))
  !print *, “Parameters ok”,p
  !read *, dum

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad,ang,zarg,targ,add,addx,addy,addz,a,G)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-7)
      ang = atan2(y(i),x(i))
      zarg = (z(i)-((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
      targ = 2*ang/dphi

      !print *, x(i),y(i),z(i)
      !read *, dum
      
      !Compute Legendre additive terms
      add = 0.
      addx = 0.
      addy = 0.
      addz = 0.
      do a=1,cnum
        add = add + coeff(a)*legendre(zarg,axial(a))*legendre(targ,az(a))
        !print *, add
        !read *, dum
        addx = addx - coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(y(i)/(y(i)**2+x(i)**2))
        !print *, addx
        addy = addy + coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(x(i)/(y(i)**2+x(i)**2))
        !print *, addy
        addz = addz + coeff(a)*legendrep(zarg,axial(a))*legendre(targ,az(a))*2/(zmax-zmin)
        !print *, addz
        !read *, dum
      end do
      G = sqrt(x(i)**2+y(i)**2) + add
      F = -(G**2 - e**2*(d+z(i))**2 + z(i)**2)
      Fx = -2*G*(x(i)/sqrt(x(i)**2+y(i)**2)+addx)
      Fy = -2*G*(y(i)/sqrt(x(i)**2+y(i)**2)+addy)
      Fz = 2*e**2*(d+z(i)) - 2*z(i) - 2*G*addz
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, G, F, Fx, Fy, Fz
      !print *, l(i), m(i), n(i), delt
      !print *, i
      !read *, dum
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Set rays outside mirror definition to NaN
    !if (abs(zarg)>1 .or. abs(targ) > 1) then
    !  x(i) = 0.
    !  y(i) = 0.
    !  z(i) = 0.
    !end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine woltersecLL

!This function traces to a Wolter-Schwarzschild primary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wsprimary(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi
  real*8 :: k,kterm,dbdx,dbdy,beta,betas,ff,g,r
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb,xi,yi,zi
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,flag,r,xi,yi,zi)
  do i=1,num
    delt = 100.
    c = 0
    !Initial ray position
    xi = x(i)
    yi = y(i)
    zi = z(i)
    do while(abs(delt)>1.e-8)
      beta = asin(sqrt(x(i)**2 + y(i)**2)/ff)
      flag = 0
      if (beta<=betas) then
        beta = betas
        flag = 1
        kterm = 0.
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
      end if
      F = -z(i) - ff*sin(betas/2)**2 + &
          ff**2*sin(beta)**2/(4*ff*sin(betas/2)**2) + &
          g*cos(beta/2)**4*(kterm)**(1-k)
      Fb = ff**2*sin(beta)*cos(beta)/(2*ff*sin(betas/2)**2) - &
           2*g*cos(beta/2)**3*sin(beta/2)*(kterm)**(1-k) + &
           g*(1-k)*cos(beta/2)*sin(beta/2)*(kterm)**(-k)*(1/k)
      Fz = -1.
      if (flag==1) then
        r = sqrt(x(i)**2 + y(i)**2)
        Fb = ff**2*sin(betas)*cos(betas)/(2*ff*sin(betas/2)**2) + &
              g*(1-k)*cos(betas/2)*sin(betas/2)*(1/k)
        F = F + (r - ff*sin(betas))*z(i)/(r**2+z(i)**2)*Fb
        Fz = Fz + (r-ff*sin(betas))*(r**2-z(i)**2)/(r**2+z(i)**2)**2*Fb
        !print *, Fb, F, Fz
        !read *, dum
      end if
      dbdx = x(i)/sqrt(1-(x(i)**2+y(i)**2)/ff**2)/ff/sqrt(x(i)**2+y(i)**2)
      dbdy = y(i)/sqrt(1-(x(i)**2+y(i)**2)/ff**2)/ff/sqrt(x(i)**2+y(i)**2)
      Fx = Fb * dbdx
      Fy = Fb * dbdy
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !if (isnan(delt)) then
      !  print *, c,x(i),y(i),z(i)
      !  print *, F, Fx, Fy, Fz
      !  print *, kterm, Fb, flag,k,tan(beta/2)**2
      !  print *, betas,ff
      !  read *, dum
      !end if
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 25 .or. isnan(delt)) then
        delt = 0.
        x(i) = xi
        y(i) = yi
        z(i) = zi
      end if
      c = c + 1
      !read *, dum
    end do
    if (c < 26) then
      Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
      ux(i) = -Fx/Fp
      uy(i) = -Fy/Fp
      uz(i) = -Fz/Fp
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wsprimary

!This function traces to a Wolter-Schwarzschild secondary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wssecondary(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi
  real*8 :: k,kterm,dbdx,dbdy,dbdz,dadb,beta,betas,ff,g,d,a
  real*8 :: gam,dbdzs,dadbs
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb,xi,yi,zi
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,dbdz,flag,c,a,dadbs,dbdzs,gam,dadb,xi,yi,zi)
  do i=1,num
    delt = 100.
    c = 0
    !Initial ray position
    xi = x(i)
    yi = y(i)
    zi = z(i)
    do while(abs(delt)>1.e-8)
      beta = atan2(sqrt(x(i)**2 + y(i)**2),z(i))
      flag = 0
      if (beta<=betas) then
        beta = betas
        kterm = 0
        a = 1/ff
        flag = 1
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
        a = (1-cos(beta))/(1-cos(betas))/ff + &
          (1+cos(beta))/(2*g)*(kterm)**(1+k)
      end if
      F = -z(i) + cos(beta)/a
      !Add corrective term to F if beta was < betas
      if (flag==1) then
        Fb = 0.
        dadbs = sin(betas)/ff/(1-cos(betas)) + &
                (k+1)*(cos(betas)+1)*tan(betas/2)/cos(betas/2)**2/2/g/k
        dbdzs = -sin(betas)**2/sqrt(x(i)**2+y(i)**2)
        gam = (-ff*sin(betas)-ff**2*cos(betas)*dadbs)*dbdzs
        F = F + gam*(z(i)-sqrt(x(i)**2+y(i)**2)/tan(betas))
        Fx = -2./tan(betas)*x(i)/sqrt(x(i)**2+y(i)**2)
        Fy = -2./tan(betas)*y(i)/sqrt(x(i)**2+y(i)**2)
        Fz = gam - 1.
        !print *, x(i), y(i), z(i)
        !print *, F, Fx, Fy, Fz
        !print *, delt
        !read *, dum
      !Otherwise, business as usual
      else
        dadb = sin(beta)/ff/(1-cos(betas)) - &
               sin(beta)/(2*g)*(kterm)**(1+k) + &
               (k+1)*(cos(beta)+1)*tan(beta/2)*kterm**k/2/g/k/(cos(beta/2)**2)
        Fb = -sin(beta)/a - cos(beta)/a**2*dadb
        dbdx = x(i)*z(i)/(x(i)**2+y(i)**2+z(i)**2)/sqrt(x(i)**2+y(i)**2)
        dbdy = y(i)*z(i)/(x(i)**2+y(i)**2+z(i)**2)/sqrt(x(i)**2+y(i)**2)
        dbdz = -sqrt(x(i)**2+y(i)**2)/(x(i)**2+y(i)**2+z(i)**2)
        Fx = Fb * dbdx
        Fy = Fb * dbdy
        Fz = -1. + Fb*dbdz
      end if
      !We have derivatives, now compute the iteration
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, x(i), y(i), z(i)
      !print *, F, Fx, Fy, Fz
      !print *, delt
      !read *, dum
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 25 .or. isnan(delt)) then
        delt = 0.
        x(i) = xi
        y(i) = yi
        z(i) = zi
      end if
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
      c = c+1
    end do
    if (c < 26) then
      Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
      ux(i) = Fx/Fp
      uy(i) = Fy/Fp
      uz(i) = Fz/Fp
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !print *, F, Fx, Fy, Fz
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wssecondary

!Intersection with SPO Cone
subroutine spoCone(x,y,z,l,m,n,ux,uy,uz,num,R0,tg)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: R0,tg
!  real*8 :: k,kterm,dbdx,dbdy,dbdz,dadb,beta,betas,ff,g,d,a
!  real*8 :: gam,dbdzs,dadbs
!  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb
  integer :: i
  real*8 :: A,B,C,sl,det,t1,t2

  !Loop through rays and trace to mirror
  sl = tan(tg)
  !$omp parallel do private(i,A,B,C,det,t1,t2)
  do i=1,num
    !Solve quadratic equation for ray advancement distance
    A = n(i)**2*sl**2 - m(i)**2 - l(i)**2
    B = 2*n(i)*sl*R0 + 2*z(i)*sl**2*n(i) - 2*x(i)*l(i) - 2*y(i)*m(i)
    C = R0**2 + 2*sl*R0*z(i) + z(i)**2*sl**2 - x(i)**2 - y(i)**2
    det = B**2 - 4*A*C
    if (det .ge. 0) then
      t1 = (-B + sqrt(det))/(2*A)
      t2 = (-B - sqrt(det))/(2*A)
      if (abs(t2) < abs(t1)) then
        t1 = t2
      end if
      !Set new ray position
      x(i) = x(i) + t1*l(i)
      y(i) = y(i) + t1*m(i)
      z(i) = z(i) + t1*n(i)
      !Set up surface normal
      ux(i) = -x(i)/sqrt(x(i)**2+y(i)**2)*cos(tg)
      uy(i) = -y(i)/sqrt(x(i)**2+y(i)**2)*cos(tg)
      uz(i) = sin(tg)!*abs(z(i))/z(i)
    else
      l(i) = 0.
      m(i) = 0.
      n(i) = 0.
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !print *, F, Fx, Fy, Fz
    !read *, dum
  end do
  !$omp end parallel do

end subroutine spoCone

!Construct an ellipsoid primary but with Legendre-Legendre
!deformations. Define Legendre and Legendre derivative functions.
!Pass in coeff, axial order, and azimuthal order as in Zemax implementation
subroutine ellipsoidWoltLL(x,y,z,l,m,n,ux,uy,uz,num,r0,z0,psi,S,zmax,zmin,dphi,coeff,axial,az,cnum)
  !Declarations
  implicit none
  integer, intent(in) :: num,cnum
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: r0,z0,psi,S,zmax,zmin,dphi,coeff(cnum)
  integer, intent(in) :: axial(cnum),az(cnum)
  real*8 :: P,ff,bq,cq,aa,bb,e
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,rad,add,addx,addy,addz,G
  integer :: i,a
  real*8 :: pi,legendre,legendrep,zarg,targ,ang,dummy,zfoc

  !Compute telescope parameters
  pi = acos(-1.)
  P = R0/sin((psi*asin(R0/z0)-asin(R0/S))/(1+psi))
  ff = (S+P)/2.
  bq = -(R0**2+(ff-P)**2+ff**2)
  cq = ff**2*(ff-P)**2
  aa = sqrt((-bq+sqrt(bq**2-4*cq))/2.)
  bb = sqrt(aa**2-ff**2)
  e = ff/aa
  zfoc = ff-P+z0

  !Loop through rays and trace to mirror
  !$omp parallel do private(delt,F,Fx,Fy,Fz,Fp,rad,ang,zarg,targ,add,addx,addy,addz,a,G)
  do i=1,num
    delt = 100.
    do while(abs(delt)>1.e-10)
      ang = atan2(y(i),x(i))
      zarg = (z(i)-((zmax+zmin)/2.)) / ((zmax-zmin)/2.)
      targ = 2*ang/dphi
      
      !Compute Legendre additive terms
      add = 0.
      addx = 0.
      addy = 0.
      addz = 0.
      do a=1,cnum
        add = add + coeff(a)*legendre(zarg,axial(a))*legendre(targ,az(a))
        addx = addx - coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(y(i)/(y(i)**2+x(i)**2))
        addy = addy + coeff(a)*legendre(zarg,axial(a))*legendrep(targ,az(a))*(2/dphi)*(x(i)/(y(i)**2+x(i)**2))
        addz = addz + coeff(a)*legendrep(zarg,axial(a))*legendre(targ,az(a))*2/(zmax-zmin)
      end do
      G = sqrt(x(i)**2+y(i)**2) + add
      F = (z(i)-zfoc)**2/aa**2 + G**2/bb**2 - 1.
      Fx = 2*G/bb**2*(x(i)/sqrt(x(i)**2+y(i)**2)+addx)
      Fy = 2*G/bb**2*(y(i)/sqrt(x(i)**2+y(i)**2)+addy)
      Fz = 2*(z(i)-zfoc)/aa**2 + (2*G/bb**2)*(addz)
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, 'Position: ',x(i), y(i), z(i)
      !print *, 'Function: ', F, Fx, Fy, Fz
      !print *, 'Change: ',delt, z0,aa
      !print *, z(i)-zfoc
      !read *, dummy
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
    end do
    Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
    ux(i) = Fx/Fp
    uy(i) = Fy/Fp
    uz(i) = Fz/Fp
    !Set rays outside mirror definition to NaN
    !if (abs(zarg)>1 .or. abs(targ) > 1) then
    !  x(i) = 0.
    !  y(i) = 0.
    !  z(i) = 0.
    !end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine ellipsoidWoltLL

!This function traces to a Wolter-Schwarzschild primary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wsprimaryBack(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi,thick)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi,thick
  real*8 :: k,kterm,dbdx,dbdy,beta,betas,ff,g,r,theta,x2,y2
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb,xi,yi,zi
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,flag,r,x2,y2,theta,xi,yi,zi)
  do i=1,num
    delt = 100.
    c = 0
    !Initial ray position
    xi = x(i)
    yi = y(i)
    zi = z(i)
    do while(abs(delt)>1.e-8)
      !Adjust x and y positions
      r = sqrt(x(i)**2+y(i)**2)
      theta = atan2(y(i),x(i))
      x2 = (r-thick)*cos(theta)
      y2 = (r-thick)*sin(theta)
      beta = asin(sqrt(x2**2 + y2**2)/ff)
      flag = 0
      if (beta<=betas) then
        beta = betas
        flag = 1
        kterm = 0.
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
      end if
      F = -z(i) - ff*sin(betas/2)**2 + &
          ff**2*sin(beta)**2/(4*ff*sin(betas/2)**2) + &
          g*cos(beta/2)**4*(kterm)**(1-k)
      Fb = ff**2*sin(beta)*cos(beta)/(2*ff*sin(betas/2)**2) - &
           2*g*cos(beta/2)**3*sin(beta/2)*(kterm)**(1-k) + &
           g*(1-k)*cos(beta/2)*sin(beta/2)*(kterm)**(-k)*(1/k)
      Fz = -1.
      if (flag==1) then
        r = sqrt(x2**2 + y2**2)
        Fb = ff**2*sin(betas)*cos(betas)/(2*ff*sin(betas/2)**2) + &
              g*(1-k)*cos(betas/2)*sin(betas/2)*(1/k)
        F = F + (r - ff*sin(betas))*z(i)/(r**2+z(i)**2)*Fb
        Fz = Fz + (r-ff*sin(betas))*(r**2-z(i)**2)/(r**2+z(i)**2)**2*Fb
        !print *, Fb, F, Fz
        !read *, dum
      end if
      dbdx = x2/sqrt(1-(x2**2+y2**2)/ff**2)/ff/sqrt(x2**2+y2**2)
      dbdy = y2/sqrt(1-(x2**2+y2**2)/ff**2)/ff/sqrt(x2**2+y2**2)
      Fx = Fb * dbdx
      Fy = Fb * dbdy
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, x(i),y(i),z(i)
      !print *, F, Fx, Fy, Fz
      !print *, kterm, Fb, flag,k,tan(beta/2)**2
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 25 .or. isnan(delt)) then
        delt = 0.
        x(i) = xi
        y(i) = yi
        z(i) = zi
      end if
      c = c + 1
      !read *, dum
    end do
    if (c < 26) then
      Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
      ux(i) = -Fx/Fp
      uy(i) = -Fy/Fp
      uz(i) = -Fz/Fp
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wsprimaryBack

!This function traces to a Wolter-Schwarzschild secondary mirror
!Defined by Van Speybroeck prescription
!Surface should be placed at common focus with z+ pointing toward mirrors
!If ray is within inner radius of mirror (defined by betas), it will be
!traced to minimum z position
!Code in Python wrapper must vignette such rays
subroutine wssecondaryBack(x,y,z,l,m,n,ux,uy,uz,num,alpha,z0,psi,thick)
  !Declarations
  implicit none
  integer, intent(in) :: num
  real*8 , intent(inout) :: x(num),y(num),z(num),l(num),m(num),n(num),ux(num),uy(num),uz(num)
  real*8, intent(in) :: alpha,z0,psi,thick
  real*8 :: k,kterm,dbdx,dbdy,dbdz,dadb,beta,betas,ff,g,d,a,r,theta,x2,y2
  real*8 :: gam,dbdzs,dadbs
  real*8 :: F,Fx,Fy,Fz,Fp,delt,dum,Fb,xi,yi,zi
  integer :: i, flag, c

  !Compute Chase parameters
  betas = 4*alpha
  ff = z0/cos(betas)
  g = ff / psi
  k = tan(betas/2)**2

  !Loop through rays and trace to mirror
  !$omp parallel do private(i,delt,F,Fx,Fy,Fz,Fp,Fb,kterm,beta,dbdx,dbdy,dbdz,flag,c,a,dadbs,dbdzs,gam,dadb,r,theta,x2,y2,xi,yi,zi)
  do i=1,num
    delt = 100.
    c = 0
    !Initial ray position
    xi = x(i)
    yi = y(i)
    zi = z(i)
    do while(abs(delt)>1.e-8)
      !Adjust x and y positions
      r = sqrt(x(i)**2 + y(i)**2)
      theta = atan2(y(i),x(i))
      x2 = (r-thick)*cos(theta)
      y2 = (r-thick)*sin(theta)
      beta = atan2(sqrt(x2**2 + y2**2),z(i))
      flag = 0
      if (beta<=betas) then
        beta = betas
        kterm = 0
        a = 1/ff
        flag = 1
      else
        kterm = (1/k)*tan(beta/2)**2 - 1
        a = (1-cos(beta))/(1-cos(betas))/ff + &
          (1+cos(beta))/(2*g)*(kterm)**(1+k)
      end if
      F = -z(i) + cos(beta)/a
      !Add corrective term to F if beta was < betas
      if (flag==1) then
        Fb = 0.
        dadbs = sin(betas)/ff/(1-cos(betas)) + &
                (k+1)*(cos(betas)+1)*tan(betas/2)/cos(betas/2)**2/2/g/k
        dbdzs = -sin(betas)**2/sqrt(x2**2+y2**2)
        gam = (-ff*sin(betas)-ff**2*cos(betas)*dadbs)*dbdzs
        F = F + gam*(z(i)-sqrt(x2**2+y2**2)/tan(betas))
        Fx = -2./tan(betas)*x2/sqrt(x2**2+y2**2)
        Fy = -2./tan(betas)*y2/sqrt(x2**2+y2**2)
        Fz = gam - 1.
        !print *, x(i), y(i), z(i)
        !print *, F, Fx, Fy, Fz
        !print *, delt
        !read *, dum
      !Otherwise, business as usual
      else
        dadb = sin(beta)/ff/(1-cos(betas)) - &
               sin(beta)/(2*g)*(kterm)**(1+k) + &
               (k+1)*(cos(beta)+1)*tan(beta/2)*kterm**k/2/g/k/(cos(beta/2)**2)
        Fb = -sin(beta)/a - cos(beta)/a**2*dadb
        dbdx = x2*z(i)/(x2**2+y2**2+z(i)**2)/sqrt(x2**2+y2**2)
        dbdy = y2*z(i)/(x2**2+y2**2+z(i)**2)/sqrt(x2**2+y2**2)
        dbdz = -sqrt(x2**2+y2**2)/(x2**2+y2**2+z(i)**2)
        Fx = Fb * dbdx
        Fy = Fb * dbdy
        Fz = -1. + Fb*dbdz
      end if
      !We have derivatives, now compute the iteration
      Fp = Fx*l(i) + Fy*m(i) + Fz*n(i)
      delt = -F/Fp
      !print *, x(i), y(i), z(i)
      !print *, F, Fx, Fy, Fz
      !print *, delt
      !read *, dum
      x(i) = x(i) + l(i)*delt
      y(i) = y(i) + m(i)*delt
      z(i) = z(i) + n(i)*delt
      if (c > 25 .or. isnan(delt)) then
        delt = 0.
        x(i) = xi
        y(i) = yi
        z(i) = zi
      end if
      !print *, x(i),y(i),z(i)
      !print *, F
      !print * ,delt
      !read *, dum
      c = c+1
    end do
    if (c < 26) then
      Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
      ux(i) = Fx/Fp
      uy(i) = Fy/Fp
      uz(i) = Fz/Fp
    end if
    !print *, x(i),y(i),z(i)
    !print *, ux(i),uy(i),uz(i)
    !print *, F, Fx, Fy, Fz
    !read *, dum
  end do
  !$omp end parallel do

end subroutine wssecondaryBack
