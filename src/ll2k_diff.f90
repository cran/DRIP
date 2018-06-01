!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!  This is Fortran source code of applying step edge detector LL2K    !
!  to an image. Step edge detector LL2K is proposed in the paper      !
!  Kang, Y., and Qiu, P.,`Jump detection in blurred regression        !
!  surfaces'.                                                         !
!  Creator : Yicheng Kang                                             !
!  Date: April 29, 2013                                               !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ll2k_diff(n, obsImg, bandwidth, diff)

  implicit none

  integer :: i, j, n, i1, j1, k, bandwidth
  
  double precision :: z(0:n, 0:n), z1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), x, y, temp, &
       x1, y1, ker, r00, r20, bhat, chat, bb, ra, fhat1, fhat2, ker1, &
       gradperp, u00, u01, u10, u20, u02, u11, obsImg(0:n, 0:n), &
       den1, den2, v00, v01, v10, v20, v02, v11, Au1, Au2, Au3, Av1, &
       Av2, Av3, ttemp1, ttemp2, diff(0:n, 0:n), dist, m1hat

  external :: extend, ker, ker1


  ! Assign values to parameters.
  
  k = bandwidth ! The chosen bandwidth
  ra = dble(k)/dble(n)

  ! Start to read in data.

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)

     end do
  end do

  ! Extend image to avoid boundary problems.

  Call extend(n, k, z, z1(0:(n+2*k), 0:(n+2*k)))

  ! Calculate first derivatives                       

  r00 = 0D0
  r20 = 0D0

  do i = -k, k
     do j = -k, k

        if (i**2 + j**2 <= k**2) then

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

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 bb = ker(dble(i1 - i)/dble(k), dble(j1 - j)/dble(k)) &
                      * z1(i1, j1)
                 bhat = bhat +  dble(i1 - i)/dble(n) * bb
                 chat = chat +  dble(j1 - j)/dble(n) * bb

              end if

           end do
        end do

        bhat = bhat/r20
        chat = chat/r20
        m1hat = sqrt(bhat**2 + chat**2)

        ! To estimate the line which is perpendicular to gradient.

        x = dble(i - k)/dble(n)
        y = dble(j - k)/dble(n)
        fhat1 = 0D0
        fhat2 = 0D0
        den1 = 0D0
        den2 = 0D0
        u00 = 0D0
        u10 = 0D0
        u01 = 0D0
        u11 = 0D0
        u20 = 0D0
        u02 = 0D0
        v00 = 0D0
        v10 = 0D0
        v01 = 0D0
        v11 = 0D0
        v20 = 0D0
        v02 = 0D0

        do i1 = i -k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 -j)**2 <= k**2) then

                  x1 = dble(i1 - k)/dble(n)
                  y1 = dble(j1 - k)/dble(n)
                  ttemp1 = x1 - x
                  ttemp2 = y1 - y
                  gradperp = bhat * ttemp1 + chat * ttemp2
                  dist = abs(gradperp)/m1hat
                  bb = ker(ttemp1/ra, ttemp2/ra) * ker1(dist/ra)

                  if (gradperp >= 0D0) then

                     u00 = u00 + bb
                     u10 = u10 + ttemp1 * bb
                     u01 = u01 + ttemp2 * bb
                     u11 = u11 + ttemp1 * ttemp2 * bb
                     u20 = u20 + ttemp1**2 * bb
                     u02 = u02 + ttemp2**2 * bb

                  else
                     
                     v00 = v00 + bb
                     v10 = v10 + ttemp1 * bb
                     v01 = v01 + ttemp2 * bb
                     v11 = v11 + ttemp1 * ttemp2 * bb
                     v20 = v20 + ttemp1**2 * bb
                     v02 = v02 + ttemp2**2 * bb

                  end if

               end if

            end do
         end do

         Au1 = u20 * u02 - u11*2
         Au2 = u01 * u11 - u10 * u02
         Au3 = u10 * u11 - u01 * u20
         Av1 = v20 * v02 - v11*2
         Av2 = v01 * v11 - v10 * v02
         Av3 = v10 * v11 - v01 * v20            

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 x1 = dble(i1 - k)/dble(n)
                 y1 = dble(j1 - k)/dble(n)
                 ttemp1 = x1 - x
                 ttemp2 = y1 - y
                 gradperp = bhat * ttemp1 + chat * ttemp2
                 dist = abs(gradperp)/m1hat

                 ! Start to fit two local constants.

                 if (gradperp >= 0D0) then

                    bb = (Au1 + ttemp1 * Au2 + ttemp2 * Au3) * &
                         ker(ttemp1/ra, ttemp2/ra) * ker1(dist/ra)
                    den1 = den1 + bb
                    fhat1 = fhat1 + z1(i1, j1) * bb
 
                 else

                    bb = (Av1 + ttemp1 * Av2 + ttemp2 * Av3) * &
                         ker(ttemp1/ra, ttemp2/ra) * ker1(dist/ra)
                    den2 = den2 + bb
                    fhat2 = fhat2 + z1(i1, j1) * bb

                 end if

              end if

           end do
        end do

        if (abs(den1) > 0D0) then

           fhat1 = fhat1/den1

        else

           fhat1 = z1(i, j)

        end if

        if (abs(den2) > 0D0) then

           fhat2 = fhat2/den2

        else

           fhat2 = z1(i, j)

        end if

        diff(i - k, j - k) = abs(fhat1 - fhat2)
 !       stdz = sqrt(num1/(den1**2) + num2/(den2**2)) * sigma

     end do
  end do

end subroutine ll2k_diff
