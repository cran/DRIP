!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
! This is Fortran subroutine for computing leave-one-out cross validation (CV) of  !
! jump-preserving (JP) local linear kernel (LLK) smoothing procedure (equation (6))!
! proposed in the paper                                                            !
! Qiu, P. (2009), `Jump-preserving surface reconstruction from noisy data'.        !
!                                                                                  !
! Creator: Yicheng Kang                                                            !
! Date: April 9 2013                                                               !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine JP_LLK_CV(n, obsImg, nband, bandwidth, cv)

  implicit none

  integer :: i, j, n, itemp, i1, j1, k, nband, iband, bandwidth(1:nband)
  
  double precision :: z(0:n, 0:n), z1(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), cv(1:nband), temp, &
       e1, e2, obsImg(0:n, 0:n), r00, r20, bhat, chat, bb, ra, gradperp, &
       fhat1, fhat2, ttemp1, fhat, ttemp2, temp1, temp2, ker, rplus00, rplus10, &
       rplus01, rplus11, rplus20, rplus02, ZKplus, XZKplus, YZKplus, rminus00, &
       rminus10, rminus01, rminus11, rminus20, rminus02, ZKminus, XZKminus, &
       YZKminus, bhat1, bhat2, chat1, chat2, det1, det2

  external :: extend, ker

  ! Read in the observed image.

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)

     end do
  end do

  do iband = 1, nband

     k = bandwidth(iband)
     ra = dble(k)/dble(n)

     ! Extend to avoid boundary problems.

     call extend(n, k, z, z1(0:(n+2*k), 0:(n+2*k)))


     ! Calculate CV.


     cv(iband) = 0D0

     ! Calculate first derivatives by leave-one out llk smoothing.

     r00 = 0D0
     r20 = 0D0

     do i = -k, k
        do j = -k, k

           itemp = i**2 + j**2

           if ((itemp <= k**2) .and. (itemp > 0)) then

              temp = ker(dble(i)/dble(k), dble(j)/dble(k))
              r00 = r00 + temp
              r20 = r20 + (dble(i)/dble(n))**2 * temp

           end if

        end do
     end do

     do i = k, n + k
        do j = k, n + k

           bhat = 0D0
           chat = 0D0

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 itemp = (i1 - i)**2 + (j1 - j)**2

                 if ((itemp <= k**2) .and. (itemp > 0)) then

                    bb = ker(dble(i1 - i)/dble(k), dble(j1 - j)/dble(k)) &
                         * z1(i1, j1)
                    bhat = bhat +  dble(i1 - i)/dble(n) * bb
                    chat = chat +  dble(j1 - j)/dble(n) * bb

                 end if

              end do
           end do

           bhat = bhat/r20
           chat = chat/r20

           ! Start to estimate surface by leave-one-out piece-wise LLK smoothing.

           ZKplus = 0D0
           XZKplus = 0D0
           YZKplus = 0D0

           rplus00 = 0D0
           rplus10 = 0D0
           rplus01 = 0D0
           rplus11 = 0D0
           rplus20 = 0D0
           rplus02 = 0D0

           ZKminus = 0D0
           XZKminus = 0D0
           YZKminus = 0D0

           rminus00 = 0D0
           rminus10 = 0D0
           rminus01 = 0D0
           rminus11 = 0D0
           rminus20 = 0D0
           rminus02 = 0D0

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 itemp = (i1 - i)**2 + (j1 - j)**2

                 if ((itemp <= k**2) .and. (itemp > 0)) then

                    bb = ker(ttemp1/ra, ttemp2/ra)
                    gradperp = ttemp1 * bhat + ttemp2 * chat

                    if ( gradperp >= 0D0 ) then

                       rplus00 = rplus00 + bb
                       rplus10 = rplus10 + ttemp1 * bb
                       rplus01 = rplus01 + ttemp2 * bb
                       rplus11 = rplus11 + ttemp1 * ttemp2 * bb
                       rplus20 = rplus20 + ttemp1**2 * bb
                       rplus02 = rplus02 + ttemp2**2 * bb

                       ZKplus = ZKplus + z1(i1, j1) * bb
                       XZKplus = XZKplus + ttemp1 * z1(i1, j1) * bb
                       YZKplus = YZKplus + ttemp2 * z1(i1, j1) * bb

                    else

                       rminus00 = rminus00 + bb
                       rminus10 = rminus10 + ttemp1 * bb
                       rminus01 = rminus01 + ttemp2 * bb
                       rminus11 = rminus11 + ttemp1 * ttemp2 * bb
                       rminus20 = rminus20 + ttemp1**2 * bb
                       rminus02 = rminus02 + ttemp2**2 * bb

                       ZKminus = ZKminus + z1(i1, j1) * bb
                       XZKminus = XZKminus + ttemp1 * z1(i1, j1) * bb
                       YZKminus = YZKminus + ttemp2 * z1(i1, j1) * bb

                    end if

                 end if

              end do
           end do

           det1 = rplus00 * rplus20 * rplus02 + 2D0 * rplus10 * rplus01 * rplus11 - rplus01**2 * rplus20 -&
                rplus11**2 * rplus00 - rplus10**2 * rplus02
           det2 = rminus00 * rminus20 * rminus02 + 2D0 * rminus10 * rminus01 * rminus11 - rminus01**2 * rminus20 -&
                rminus11**2 * rminus00 - rminus10**2 * rminus02

           fhat1 = ((rplus02 * rplus20 - rplus11**2) * ZKplus + (rplus01 * rplus11 - rplus10 * rplus02) * XZKplus + &
                (rplus10 * rplus11 - rplus01 * rplus20) * YZKplus)/det1
           fhat2 = ((rminus02 * rminus20 - rminus11**2) * ZKminus + (rminus01 * rminus11 - rminus10 * rminus02) * XZKminus + &
                (rminus10 * rminus11 - rminus01 * rminus20) * YZKminus)/det2

           bhat1 = ((rplus11 * rplus01 - rplus10 * rplus02) * ZKplus + (rplus00 * rplus02 - rplus01**2) * XZKplus + &
                (rplus01 * rplus10 - rplus00 * rplus11) * YZKplus)/det1
           bhat2 = ((rminus11 * rminus01 - rminus10 * rminus02) * ZKminus + (rminus00 * rminus02 - rminus01**2) * XZKminus + &
                (rminus01 * rminus10 - rminus00 * rminus11) * YZKminus)/det2

           chat1 = ((rplus10 * rplus11 - rplus01 * rplus20) * ZKplus + (rplus01 * rplus10 - rplus00 * rplus11) * XZKplus + &
                (rplus00 * rplus20 - rplus10**2) * YZKplus)/det1
           chat2 = ((rminus10 * rminus11 - rminus01 * rminus20) * ZKminus + (rminus01 * rminus10 - rminus00 * rminus11) * &
                XZKminus + (rminus00 * rminus20 - rminus10**2) * YZKminus)/det2

           ! Start to calculate weighted residual mean squares.

           e1 = 0D0  
           e2 = 0D0
           temp1 = 0D0
           temp2 = 0D0

           do i1 = i -k, i + k
              do j1 = j -k, j + k

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 itemp = (i1 - i)**2 + (j1 - j)**2

                 if ((itemp <= k**2) .and. (itemp > 0)) then

                    bb = ker(ttemp1/ra, ttemp2/ra)
                    gradperp = ttemp1 * bhat + ttemp2 * chat

                    if (gradperp >= 0D0 ) then

                       e1 = e1 + (z1(i1, j1) - fhat1- ttemp1 * bhat1 - ttemp2 * chat1)**2 * bb

                    else

                       e2 = e2 + (z1(i1, j1) - fhat2 - ttemp1 * bhat2 - ttemp2 * chat2)**2 * bb

                    end if

                 end if

              end do
           end do

           e1 = e1/temp1
           e2 = e2/temp2

           if ( e1 >= e2 ) then

              fhat = fhat2

           else

              fhat = fhat1

           end if

           ! Residuals obtained

           cv(iband) = cv(iband) + (z1(i, j) - fhat)**2

        end do
     end do

     cv(iband) = cv(iband)/dble((n + 1)**2)

  end do

end subroutine JP_LLK_CV
