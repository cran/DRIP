!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!  This is Fortran source code of applying jump detector LC2K         !
!  to an image. Jump detector LC2K is proposed in the paper           !
!  Kang, Y., and Qiu, P.,`Jump detection in blurred regression        !
!  surfaces'.                                                         !
!  Creator : Yicheng Kang                                             !
!  Date: April 21, 2013                                               !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lc2k_diff(n, obsImg, bandwidth, diff)

  implicit none

  integer :: i, j, n, i1, j1, k, bandwidth
  
  double precision :: z(0:n, 0:n), z1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), x, y, temp, &
       x1, y1, ker, r00, r20, bhat, chat, bb, ra, fhat1, fhat2, ker1, &
       gradperp, diff(0:n, 0:n), obsImg(0:n, 0:n), den1, den2, m1hat, &
       dist

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

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              x1 = dble(i1 - k)/dble(n)
              y1 = dble(j1 - k)/dble(n)

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 gradperp = bhat * (x1 - x) + chat * (y1 - y)
                 dist = abs(gradperp)/m1hat
                 bb = ker((x1 - x)/ra, (y1 - y)/ra) * ker1(dist/ra)

                 ! Start to fit two local constants.

                 if (gradperp >= 0D0) then

                    fhat1 = fhat1 + z1(i1, j1) * bb
                    den1 = den1 + bb

                 else

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

     end do
  end do

end subroutine lc2k_diff
