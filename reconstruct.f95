subroutine reconstruct(xang,yang,xdim,ydim,criteria,h,phase,phasec,maxiter)
    !Declarations
    integer, intent(in) :: xdim,ydim,maxiter
    real*8, dimension(xdim,ydim), intent(inout) :: xang,yang,phase
    !real*8, dimension(xdim,ydim) :: xang,yang,phaseo
    real*8, intent(in) :: criteria,h
    real*8, dimension(xdim,ydim), intent(out) :: phasec
    real*8 :: pi, w, bk, psum, rms, goodpix
    real*8 :: xplus,xneg,yplus,yneg,pxplus,pxneg,pyplus,pyneg
    integer :: xi,yi,compute,nannum,counter
    real*8 :: dum

    pi = 3.1415926535897931
    w = 2/(1+sin(pi/(sqrt(dble(xdim)*dble(ydim))+1)))

    !Set phasec to phase
    phasec = phase
    !phaseo = phase
    !xang = xang2
    !yang = yang2

    !Enter reconstruction loop, only loop through non-nan elements
    counter = 0
    do
        !Loop through phase array and update in Gauss-Seidel fashion
        rms = 0. !Reset RMS squared norm accumulator
        nannum = 0 !Reset RMS counter
        do xi = 2,xdim-1
            do yi = 2,ydim-1
                compute = 1
                !If NaN, continue to next element
                if (phasec(xi,yi)==100.) then
                    compute = 0
                end if

                if (compute==1) then
                    !Make array of nearest neighbors
                    yplus = yang(xi,yi+1)
                    if (yplus==100.) then
                        yplus=-yang(xi,yi)
                    end if
                    yneg = yang(xi,yi-1)
                    if (yneg==100.) then
                        yneg=-yang(xi,yi)
                    end if
                    xplus = xang(xi+1,yi)
                    if (xplus==100.) then
                        xplus=-xang(xi,yi)
                    end if
                    xneg = xang(xi-1,yi)
                    if (xneg==100.) then
                        xneg=-xang(xi,yi)
                    end if

                    !Compute slope term
                    bk = .5*(yplus-yneg+xplus-xneg)*h

                    goodpix = 4.
                    !Compute nearest neighbor sum
                    pyplus = phasec(xi,yi+1)
                    if (pyplus==100.) then
                        pyplus=0.
                        !pyplus = phasec(xi,yi)
                        goodpix = goodpix - 1
                    end if
                    pyneg = phasec(xi,yi-1)
                    if (pyneg==100.) then
                        pyneg=0.
                        !pyneg = phasec(xi,yi)
                        goodpix = goodpix - 1
                    end if
                    pxplus = phasec(xi+1,yi)
                    if (pxplus==100.) then
                        pxplus=0.
                        !pxplus = phasec(xi,yi)
                        goodpix = goodpix - 1
                    end if
                    pxneg = phasec(xi-1,yi)
                    if (pxneg==100.) then
                        pxneg=0.
                        !pxneg = phasec(xi,yi)
                        goodpix = goodpix - 1
                    end if
                    psum = pyplus + pyneg + pxplus + pxneg

                    !If lenslet has less than 2 good neighbors, set to NaN to discard
                    if (goodpix == 0.) then
                        phasec(xi,yi) = 100.
                        phase(xi,yi) = 100.
                        xang(xi,yi) = 100.
                        yang(xi,yi) = 100.
                        compute = 0
                    end if
                end if

                if (compute==1) then
                    !Correct phase
                    phasec(xi,yi) = phasec(xi,yi) + w*((psum+bk)/goodpix-phasec(xi,yi))

                    !Accumulate squared norms
                    nannum = nannum + 1
                    rms = rms + (phasec(xi,yi)-phase(xi,yi))**2
                end if

            end do

        end do

        !Determine rms improvement
        print *, rms, nannum
        rms = sqrt(rms/nannum)
        print *, rms/nannum
        phase = phasec !Copy new phase array to old phase array

        if (rms < criteria) then
            exit
        end if

        counter = counter + 1
        if (counter > maxiter) then
            exit
        end if

    end do

end subroutine reconstruct

!This function bins a list of x,y intercepts and l,m cosines
!and constructs the arrays that reconstruct requires
!Inputs are 1D x,y,l,m vectors of size num
!2D arrays of xang,yang, and phase with size xdim*ydim, init to 0
!double binsize indicating physical size of lenslets
subroutine southwellbin(x,y,l,m,num,binsize,xang,yang,phase,xdim,ydim)
  !Declarations
  implicit none
  integer, intent(in) :: num,xdim,ydim
  real*8, intent(in) :: x(num),y(num),l(num),m(num),binsize
  real*8, intent(out) :: xang(xdim,ydim),yang(xdim,ydim),phase(xdim,ydim)
  integer :: i,j,xb,yb,accum(xdim,ydim)

  !accum(:,:) = 0
  do i=1,xdim
    do j=1,ydim
      accum(i,j) = 0
    end do
  end do

  !Loop through rays
  !$omp parallel do private(xb,yb)
  do i=1,num
    !Figure out proper bin
    if (mod(xdim,2)==0) then
      xb = floor(x(i)/binsize) + xdim/2 + 1 + 1
    else
      xb = floor((x(i)+binsize/2)/binsize) + (xdim-1)/2 + 1
    end if
    if (mod(ydim,2)==0) then
      yb = floor(y(i)/binsize) + ydim/2 + 1 + 1
    else
      yb = floor((y(i)+binsize/2)/binsize) + (ydim-1)/2 + 1
    end if
    !Add l and m to proper array
    xang(xb,yb) = xang(xb,yb) + l(i)
    yang(xb,yb) = yang(xb,yb) + m(i)
    !Accumulate total number of rays in WFS bin
    accum(xb,yb) = accum(xb,yb) + 1
  end do
  !$omp end parallel do

  !Normalize xang and yang
  !Set invalid pixels to 100.
  !$omp parallel do
  do i=1,xdim
    do j=1,ydim
      if (accum(i,j)==0) then
        phase(i,j) = 100.
        xang(i,j) = 100.
        yang(i,j) = 100.
      else
        xang(i,j) = tan(asin(xang(i,j)/accum(i,j)))
        yang(i,j) = tan(asin(yang(i,j)/accum(i,j)))
      end if
    end do
  end do
  !$omp end parallel do

end subroutine southwellbin
