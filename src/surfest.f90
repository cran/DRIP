!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
! This Fortran subroutine estimates 2-dimensional continuous surfaces by conventional   ! 
! local linear kernel smoothing. Estimated can be used to obtain residuals or noisy     !
! level.                                                                                !
!                                                                                       !
! Creator: Yicheng Kang                                                                 !
! Date: Sep 06, 2015                                                                    !
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine surfest(n, data, bandw, llkbw, sigma, fitted, resid, optb1)

  implicit none

  integer :: n, i, j, k, l, bandw, i1, j1, optb

  double precision :: llkbw(1:bandw)

  integer, parameter :: ext_size = 20
  
  double precision :: data(0:n, 0:n), z(0:n, 0:n), z1(0:(n+ 2*ext_size), 0:(n+2*ext_size)), &
       fbhat(0:(n+ext_size), 0:(n+ext_size)), sigma, resid(0:n, 0:n), optb1, &
       fitted(0:n, 0:n), cv(1:bandw), r00,   &
       r10, r01, r11, r02, r20, det1, det2, det3, det, ra,  &
       ttemp1, ttemp2, ttemp, ker, mincv, fhlin, r1(0:(n+ext_size), 0:(n+ext_size)), &
       x, y, x1, y1, bb, tolerance

  external:: ker, extend

  !! Assign value to parameters.

  tolerance = 1D-8

  !! Read in data

  do i = 0, n
     do j = 0, n

        z(i, j) = data(i, j)

     end do
  end do

  !! Extend to avoid boundary problems

  call extend(n, ext_size, z, z1) ! This cannot be done if bandwidth >=20 or sample size <= 20. 


  !! Start bandwidth selection by cross validation. !!!!

  do l = 1, bandw

     k = int(dble(n) * llkbw(l))
     ra = dble(k)/dble(n)

     do i = ext_size, n + ext_size
        do j = ext_size, n + ext_size

           x = dble(i)/dble(n)
           y = dble(j)/dble(n)
           fhlin = 0D0
           r00 = 0D0
           r10 = 0D0
           r01 = 0D0
           r11 = 0D0
           r20 = 0D0
           r02 = 0D0

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 x1 = dble(i1)/dble(n)
                 y1 = dble(j1)/dble(n)
                 ttemp1 = x1 - x
                 ttemp2 = y1 - y
                 ttemp = dsqrt(ttemp1**2 + ttemp2**2)

                 if (ttemp .LE. ra .and. ttemp > 0D0) then

                    !! Start to fit leave-one-out local plane.

                    bb = ker(ttemp1/ra, ttemp2/ra)
                    r00 = r00 + bb
                    r10 = r10 + ttemp1 * bb
                    r01 = r01 + ttemp2 * bb
                    r11 = r11 + ttemp1 * ttemp2 * bb
                    r20 = r20 + ttemp1**2 * bb
                    r02 = r02 + ttemp2**2 * bb

                 end if

              end do
           end do

           det = r00 * r20 * r02 + r10 * r11 * r01 + r01 * r10 * r11 - &
                r01 * r20 * r01 - r10 * r10 * r02 - r00 * r11 * r11

           det1 = r20 * r02 - r11**2
           det2 = r01 * r11 - r10 * r02
           det3 = r10 * r11 - r01 * r20

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 x1 = dble(i1)/dble(n)
                 y1 = dble(j1)/dble(n)
                 ttemp1 = x1 - x
                 ttemp2 = y1 - y
                 ttemp = dsqrt(ttemp1**2 + ttemp2**2)

                 if (ttemp .LE. ra .and. ttemp > 0D0) then

                    fhlin = fhlin + z1(i1, j1) * &
                         ker(ttemp1/ra, ttemp2/ra) * &
                         (det1 + det2 * ttemp1 + det3 * ttemp2)

                 end if

              end do
           end do

           fhlin = fhlin/det
           fbhat(i, j) = fhlin
           r1(i, j) = z1(i, j) - fbhat(i, j) ! Leave-one-out residuals!

        end do
     end do

     ! Calculate cross validation.

     cv(l) = 0D0

     do i = ext_size, n + ext_size
        do j = ext_size, n + ext_size

           cv(l) = cv(l) + r1(i, j)**2

        end do
     end do

  end do

  mincv = minval(cv)
  optb = 0

  do l = 1, bandw

     if ( abs(cv(l) - mincv) <= tolerance) then

        optb = int(llkbw(l) * n) ! Optimal bandwidth.
        optb1 = llkbw(l)

     end if

  end do

  ra = dble(optb)/dble(n)

  !! End of bandwidth selection


  !! Start to estimate surface with optimal bandwidth.

  sigma = 0D0

  do i = ext_size, ext_size + n
     do j = ext_size, ext_size + n

        fitted(i - ext_size, j - ext_size) = 0D0
        fhlin = 0D0
        r00 = 0D0
        r10 = 0D0
        r01 = 0D0
        r11 = 0D0
        r20 = 0D0
        r02 = 0D0

        do i1 = i - optb, i + optb
           do j1 = j - optb, j + optb

              ttemp1 = dble(i1 - i)/dble(n)
              ttemp2 = dble(j1 - j)/dble(n)

              if ( ttemp1**2 + ttemp2**2 <= ra**2 ) then

                 bb = ker(ttemp1/ra, ttemp2/ra)
                 r00 = r00 + bb
                 r10 = r10 + ttemp1 * bb
                 r01 = r01 + ttemp2 * bb
                 r11 = r11 + ttemp1 * ttemp2 * bb
                 r20 = r20 + ttemp1**2 * bb
                 r02 = r02 + ttemp2**2 * bb

              end if

           end do
        end do

        det = r00 * r20 * r02 + r10 * r11 * r01 + r01 * r10 * r11 - &
             r01 * r20 * r01 - r10 * r10 * r02 - r00 * r11 * r11

        det1 = r20 * r02 - r11**2
        det2 = r01 * r11 - r10 * r02
        det3 = r10 * r11 - r01 * r20


        do i1 = i -optb, i + optb
           do j1 = j -optb, j + optb

              ttemp1 = dble(i1 - i)/dble(n)
              ttemp2 = dble(j1 - j)/dble(n)

              if ( ttemp1**2 + ttemp2**2 <= ra**2 ) then

                 fhlin = fhlin + z1(i1, j1) * &
                      ker(ttemp1/ra, ttemp2/ra) * &
                      (det1 + det2 * ttemp1 + det3 * ttemp2)

              end if

           end do
        end do

        fhlin = fhlin/det
        fitted(i - ext_size, j - ext_size) = fhlin

        !! Residuals obtained from conventional local linear kernel smoothing.

        resid(i - ext_size, j - ext_size) = z1(i, j) - fitted(i - ext_size, j - ext_size) 
        sigma = sigma + resid(i - ext_size, j - ext_size)**2

     end do
  end do

  !! Estimate sigma.

  sigma = sigma/dble((n + 1)**2)
  sigma = sqrt(sigma)

end subroutine surfest
