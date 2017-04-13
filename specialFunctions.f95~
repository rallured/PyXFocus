!Factorial Function
real*8 function factorial(n)
    implicit none
    integer, intent(in) :: n
    real*8 :: Ans
    integer :: i

    Ans = 1
    do i=1,n
        Ans = Ans * dble(i)
    end do

    factorial = Ans
end function factorial

!Function to compute radial Zernike polynomial at radius rho of order n,m
real*8 function radialpoly(rho,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho
    integer, intent(in) :: n,m
    real*8 :: output
    real*8 :: const, factorial
    integer :: j

    !n,m are assumed to be valid (m>=n, n-abs(m) is even)
    !Following Python module convention, negative m indicates sinusoidal Zernike
    output = 0.
    if (rho<=1) then
        do j = 0,(n-abs(m))/2
            const = (-1)**j*factorial(n-j)/factorial(j)/factorial((n+m)/2-j)/factorial((n-m)/2-j)
            output = output + const*rho**(n-2*j)
        end do
    end if

    radialpoly = output

end function radialpoly

!This function makes use of the q recursive method to compute a Zernike
!Radial polynomial
recursive function fastradpoly(rho,n,m) result(output)
  !Declarations
  implicit none
  real*8, intent(in) :: rho,n,m
  !integer, intent(in) :: n,m
  real*8 :: output,r2,r4,h1,h2,h3,tempm

  !Handle rho=0
  if (rho==0) then
    if (m==0) then
      output = (-1.)**(n/2)
      return
    else
      output = 0.
      return
    end if
  end if

  !Three possibilities, m=n, m=n-2, m<n-2
  if (m==n) then
    output = rho**n
    return
  else if (m==n-2) then
    output = n*rho**n - (n-1)*rho**(n-2)
    return
  else
    tempm = n-4
    r4 = rho**n
    r2 = n*rho**n - (n-1)*rho**(n-2)
    do while(tempm >= m)
      h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
      h2 = h3*(n+tempm+4)*(n-tempm-2)/4./(tempm+3) + (tempm+2)
      h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.
      output = h1*r4 + (h2 + h3/rho**2)*r2
      if (tempm == m) then
        return
      end if
      r4 = r2
      r2 = output
      tempm = tempm - 2
    end do
  end if

end function fastradpoly

!This function makes use of the q recursive method to compute a Zernike
!Radial polynomial DERIVATIVE
recursive function fastradder(rho,n,m) result(output)
  !Declarations
  implicit none
  real*8, intent(in) :: rho
  real*8, intent(in) :: n,m
  real*8 :: output,r2,r4,h1,h2,h3,fastradpoly,tempm

  !Handle rho=0
  if (rho==0) then
    if (m==1) then
      output = (-1.)**((n-1)/2)*(n+1)/2
      return
    else
      output = 0.
      return
    end if
  end if  

  !Three possibilities, m=n, m=n-2, m<n-2
  if (m==n) then
    output = n*rho**(n-1)
    return
  else if (m==n-2) then
    output = n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
    return
  else
    tempm = n-4
    r4 = n*rho**(n-1)
    r2 = n*(n*rho**(n-1)) - (n-1)*(n-2)*rho**(n-3)
    do while(tempm >= m)
      h3 = -4*(tempm+2)*(tempm+1)/(n+tempm+2)/(n-tempm)
      h2 = h3*(n+tempm+4)*(n-tempm-2)/4./(tempm+3) + (tempm+2)
      h1 = .5*(tempm+4)*(tempm+3) - (tempm+4)*h2 + h3*(n+tempm+6)*(n-tempm-4)/8.
      !Might need to literally copy and paste fastradpoly here to save on function call...
      output = h1*r4 + (h2 + h3/rho**2)*r2 - 2*h3*fastradpoly(rho,n,tempm+2)/rho**3
      if (tempm == m) then
        return
      end if
      r4 = r2
      r2 = output
      tempm = tempm - 2
    end do
  end if

end function fastradder

!This function will return a vector of Zernike polynomial evaluations for
!a point (rho,theta) and a number N of standard Zernike polynomials
!Start at n=0, and loop up in radial order
!For each radial order, start at Znn, and recursively calculate down
!Adding azimuthal dependence is easy, just note what overall index you're at
!and add sinmt if odd and cosmt if even, where Z1 = 1, so Z2=2rcost
!Still need special handling of rho=0
subroutine zernset(rho,theta,rorder,aorder,znum,polyout,derrho,dertheta)
  !Declarations
  implicit none
  integer, intent(in) :: znum
  integer, intent(in) :: rorder(znum),aorder(znum)
  real*8, intent(in) :: rho,theta
  real*8, intent(out) :: polyout(znum),derrho(znum),dertheta(znum)
  real*8, allocatable :: rnm(:,:),rprime(:,:)
  integer :: i,j,radnum,tznum
  real*8 :: n,m,mm,h1,h2,h3,norm
  real*8 :: fastradpoly,fastradder,zernthetader,zernrhoder,zernike

  !Allocate temp output array for full sets of radial polynomials
  !Will slice off unused terms at the end
  tznum = 1
  radnum = 1
  do while(tznum < znum)
    tznum = tznum + (radnum + 1)
    radnum = radnum + 1
  end do
  !print *, tznum, " ", radnum
  !Allocate radnum by radnum arrays for radial polynomials
  !This is much larger than necessary, but for practical purposes this should be ok
  allocate(rnm(radnum,radnum))
  allocate(rprime(radnum,radnum))

  !Compute radial polynomials and derivatives
  do i=1,radnum
    n = dble(i)-1
    do j=int(n)+1,1,-2
      m = dble(j)-1
      !Handle case rho=0
      if (rho==0) then
        if (m==1) then
          rprime(i,j) = (-1.)**((n-1)/2)*(n+1)/2
          rnm(i,j) = 0.
        else if (m==0) then
          rprime(i,j) = 0.
          rnm(i,j) = (-1.)**(n/2)
        else
          rprime(i,j) = 0.
          rnm(i,j) = 0.
        end if
      !Handle general case
      else
        if (n==m) then
          rnm(i,j) = rho**n
          rprime(i,j) = n*rho**(n-1)
        else if (m==n-2) then
          rnm(i,j) = n*rnm(i,i) - (n-1)*rnm(i-2,i-2)
          rprime(i,j) = n*rprime(i,i) - (n-1)*rprime(i-2,i-2)
        else
          h3 = -4*(m+2)*(m+1)/(n+m+2)/(n-m)
          h2 = h3*(n+m+4)*(n-m-2)/4./(m+3) + (m+2)
          h1 = .5*(m+4)*(m+3) - (m+4)*h2 + h3*(n+m+6)*(n-m-4)/8.
          rnm(i,j) = h1*rnm(i,j+4) + (h2+h3/rho**2)*rnm(i,j+2)
          rprime(i,j) = h1*rprime(i,j+4) + (h2+h3/rho**2)*rprime(i,j+2) - 2*h3/rho**3*rnm(i,j+2)
        end if
      end if
      !print *, i-1, " ", j-1
      !print *, rnm(i,j), " ", fastradpoly(rho,dble(i-1),dble(j-1))
      !print *, rprime(i,j), " ", fastradder(rho,dble(i-1),dble(j-1))
    end do
  end do

  !Radial terms are computed. Now construct Zernike, Zernderrho, and Zerndertheta
  !Use standard order up to znum
  do i=1,znum
    !Rorder and Aorder are passed to this function
    !Simply construct if statement on m, and compute all three!
    n = rorder(i)
    mm = aorder(i)
    m = abs(mm)
    norm = sqrt(2*(n+1))
    if (mm<0) then
        polyout(i) = norm * rnm(int(n)+1,int(m)+1) * sin(m*theta)
        derrho(i) = norm * rprime(int(n)+1,int(m)+1) * sin(m*theta)
        dertheta(i) = norm * rnm(int(n)+1,int(m)+1) * cos(m*theta) * m
    else if (mm>0) then
        polyout(i) = norm * rnm(int(n)+1,int(m)+1) * cos(m*theta)
        derrho(i) = norm * rprime(int(n)+1,int(m)+1) * cos(m*theta)
        dertheta(i) = -norm * rnm(int(n)+1,int(m)+1) * sin(m*theta) * m
    else
        polyout(i) = norm*sqrt(0.5)*rnm(int(n)+1,int(m)+1)
        derrho(i) = norm*sqrt(0.5)*rprime(int(n)+1,int(m)+1)
        dertheta(i) = 0.
    end if
    !print *, dertheta(i), " ",zernthetader(rho,theta,int(n),int(mm))
  end do

end subroutine zernset

!Function to compute full Zernike polynomial at radius rho and angle theta of order n,m
!Following Python Zernike mod conventions
real*8 function zernike(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: fastradpoly,rnm,znm,norm

    !Compute radial polynomial
    rnm = fastradpoly(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * sin(m*theta)
    else if (m>0) then
        znm = norm * rnm * cos(m*theta)
    else
        znm = norm*sqrt(0.5)*rnm
    end if

    !Return Zernike polynomial
    zernike = znm

end function zernike

!Function to compute radial polynomial derivative
real*8 function radialder(rho,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho
    integer, intent(in) :: n,m
    real*8 :: output
    real*8 :: const, factorial
    integer :: j
    !n,m are assumed to be valid (m>=n, n-abs(m) is even)
    !Following Python module convention, negative m indicates sinusoidal Zernike
    output = 0.
    if (rho<=1) then
        do j = 0,(n-abs(m))/2
            const = (-1)**j*factorial(n-j)/factorial(j)/factorial((n+m)/2-j)/factorial((n-m)/2-j)
            output = output + const*(n-2*j)*rho**(n-2*j-1)
        end do
    end if

    radialder = output

end function radialder

!Function to compute partial Zernike derivative wrt rho
real*8 function zernrhoder(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: fastradder, rnm, norm, znm

    !Compute radial polynomial
    rnm = fastradder(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * sin(abs(m)*theta)
    else if (m>0) then
        znm = norm * rnm * cos(m*theta)
    else
        znm = norm*sqrt(0.5)*rnm
    end if

    !Return Zernike polynomial
    zernrhoder = znm

end function zernrhoder

!Function to compute partial Zernike derivative wrt theta
real*8 function zernthetader(rho,theta,n,m)
    !Declarations
    implicit none
    real*8, intent(in) :: rho,theta
    integer, intent(in) :: n,m
    real*8 :: rnm, znm, fastradpoly, norm

    !Compute radial polynomial
    rnm = fastradpoly(rho,dble(n),abs(dble(m)))

    !Compute full Zernike polynomial
    norm = sqrt(2*(dble(n)+1))
    if (m<0) then
        znm = norm * rnm * abs(m) * cos(abs(m)*theta)
    else if (m>0) then
        znm = -norm * rnm * m * sin(m*theta)
    else
        znm = 0.
    end if

    !Return Zernike polynomial
    zernthetader = znm

end function zernthetader

!This function computes a Legendre polynomial of order n
real*8 function legendre(x,n)
  !Declarations
  implicit none
  real*8, intent(in) :: x
  integer, intent(in) :: n
  integer :: i
  real*8 :: factorial,x2

  if (abs(x) > 1.) then
    x2 = x/abs(x)
  else
    x2 = x
  end if

  legendre = 0.
  if (n==0) then
    legendre = 1.
  else
    do i=0,floor(real(n)/2)
      legendre = legendre + (-1)**(i)*factorial(2*n-2*i)/factorial(i)/factorial(n-i)/factorial(n-2*i)/2**n*x2**(n-2*i)
    end do
  end if

end function legendre

!This function computes a Legendre polynomial first derivative of order n
real*8 function legendrep(x,n)
  !Declarations
  implicit none
  real*8, intent(in) :: x
  integer, intent(in) :: n
  integer :: i
  real*8 :: factorial

  legendrep = 0.
  if (n==0) then
    legendrep = 0.
  else if (n==1) then
    legendrep = 1.
  else if (x==0. .and. mod(n,2)==0) then
    legendrep = 0.
  else
    do i=0,floor(real(n)/2)
      legendrep = legendrep + (-1)**(i)*factorial(2*n-2*i)/factorial(i)/factorial(n-i)/factorial(n-2*i)/2**n*(n-2*i)*x**(n-2*i-1)
    end do
  end if  

  if (abs(x) > 1.) then
    legendrep = 0.
  end if

end function legendrep
