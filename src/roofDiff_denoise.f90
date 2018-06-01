!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
! Fortran source code detects roof/valley edges when NO blur!
! is present using the edge detector proposed in the paper  !
! Kang, Y. and Qiu, P. 'Blind Image Deblurring Using Jump   !
!                       Regression Analysis'.               !
! Creator: Yicheng Kang                                     !
! Date: May 6, 2013                                         !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine roofDiff_denoise(n, obsImg, bandwidth, diff)

  implicit none

  integer :: n, bandwidth

  integer :: k, i, i1, j, j1
  
  double precision :: z(0:n, 0:n), obsImg(0:n, 0:n), &
       ker, temp, G1plus, ttemp2, ra, detn, detp, ln01, &
       G2plus, G3plus, z1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), r00, r20, r22, &
       G1minus, G2minus, G3minus, r40, eta1, eta2, H1plus, &
       H2plus, H3plus, H1minus, H2minus, eta3, det, KZ, &
       H3minus, dhat, ehat, fhat, bb, lp00, lp10, ln00, &
       bhat1, bhat2, chat1, chat2, hassperp, ttemp1, lp01, & 
       lp11, lp20, lp02, diff(0:n, 0:n), ln11, ln20, ln02, &
       ln10, X2KZ, Y2KZ, XYKZ, aa

  external :: extend, ker

!!!!! Assign values to certain parameters

  k = bandwidth
  ra = dble(k)/dble(n)

!!!!! Read in noisy observations

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)

     end do
  end do

  ! Extend to avoid boundary problems.

  call extend(n, k, z, z1(0:(n+2*k), 0:(n+2*k)))

  ! Calculate first and second derivatives.

  r00 = 0D0
  r20 = 0D0
  r22 = 0D0
  r40 = 0D0

  do i = -k, k
     do j = -k, k

        temp = ker(dble(i)/dble(k), dble(j)/dble(k))
        r00 = r00 + temp
        r20 = r20 + (dble(i)/dble(n))**2 * temp
        r40 = r40 + (dble(i)/dble(n))**4 * temp
        r22 = r22 + (dble(i)/dble(n))**2 * (dble(j)/dble(n))**2 &
             * temp

     end do
  end do

  eta1 = r20 * (r22 - r40)
  eta2 = r00 * r40 - r20**2
  eta3 = r20**2 - r00 * r22
  det = r00 * r40**2 + 2 * r22 * r20**2 - 2 * r20**2 * r40 - r00 * r22**2


  do i = k, n + k
     do j = k, n + k

        KZ = 0D0
        X2KZ = 0D0
        Y2KZ = 0D0
        XYKZ = 0D0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              bb = ker(dble(i1 - i)/dble(k), dble(j1 - j)/dble(k)) &
                   * z1(i1, j1)
              KZ = KZ + bb
              X2KZ = X2KZ + (dble(i1-i)/dble(n))**2 * bb
              Y2KZ = Y2KZ + (dble(i1-i)/dble(n))**2 * bb
              XYKZ = XYKZ + (dble(i1-i)/dble(n)) * (dble(j1-j)/dble(n)) * bb

           end do
        end do

        dhat = (eta1 * KZ + eta2 * X2KZ + eta3 * Y2KZ)/det * 2D0
        ehat = (eta1 * KZ + eta3 * X2KZ + eta2 * Y2KZ)/det * 2D0
        fhat = XYKZ/r22

        ! Start to detect roof edge by using LLK detector. First, every neiborhood
        ! is divided into two halves along the direction perpendicular to (fxx, fxy).


        bhat1 = 0D0
        bhat2 = 0D0

        lp00 = 0D0
        lp10 = 0D0
        lp01 = 0D0
        lp11 = 0D0
        lp20 = 0D0
        lp02 = 0D0
        ln00 = 0D0
        ln10 = 0D0
        ln01 = 0D0
        ln11 = 0D0
        ln20 = 0D0
        ln02 = 0D0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 hassperp = dhat * ttemp1 + fhat * ttemp2
                 aa = ker(ttemp1/ra, ttemp2/ra)

                 ! Start to fit two local linears with kernels.  

                 if (hassperp >= 0D0) then

                    lp00 = lp00 + aa
                    lp10 = lp10 + ttemp1 * aa
                    lp01 = lp01 + ttemp2 * aa
                    lp11 = lp11 + ttemp1 * ttemp2 * aa
                    lp20 = lp20 + ttemp1**2 * aa
                    lp02 = lp02 + ttemp2**2 * aa

                 else

                    ln00 = ln00 + aa
                    ln10 = ln10 + ttemp1 * aa
                    ln01 = ln01 + ttemp2 * aa
                    ln11 = ln11 + ttemp1 * ttemp2 * aa
                    ln20 = ln20 + ttemp1**2 * aa
                    ln02 = ln02 + ttemp2**2 * aa

                 end if
              end if

           end do
        end do

        detp = lp00 * lp20 * lp02 + lp10 * lp11 * lp01 + lp01 * &
             lp10 * lp11 - lp01 * lp20 * lp01 - lp10 * lp10 * &
             lp02 - lp00 * lp11 * lp11

        G1plus = lp11 * lp01 - lp10 * lp02
        G2plus = lp00 * lp02 - lp01**2
        G3plus = lp01 * lp10 - lp00 * lp11

        detn = ln00 * ln20 * ln02 + ln10 * ln11 * ln01 + ln01 * &
             ln10 * ln11 - ln01 * ln20 * ln01 - ln10 * ln10 * &
             ln02 - ln00 * ln11 * ln11

        G1minus = ln11 * ln01 - ln10 * ln02
        G2minus = ln00 * ln02 - ln01**2
        G3minus = ln01 * ln10 - ln00 * ln11

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 hassperp = dhat * ttemp1 + fhat * ttemp2
                 bb = ker(ttemp1/ra, ttemp2/ra)

                 if (hassperp >= 0D0) then

                    bhat1 = bhat1 + z1(i1, j1) * bb * &
                         (G1plus + G2plus * ttemp1 + G3plus &
                         * ttemp2)

                 else

                    bhat2 = bhat2 + z1(i1, j1) * bb * &
                         (G1minus + G2minus * ttemp1 + G3minus &
                         * ttemp2)

                 end if

              end if

           end do
        end do

        if (abs(detp) > 0D0) then

           bhat1 = bhat1/detp

        else

           bhat1 = z1(i, j)

        end if

        if (abs(detn) > 0D0) then

           bhat2 = bhat2/detn

        else

           bhat2 = z1(i, j)

        end if


        ! Start to detect roof edge by using LLK detector. First, every neiborhood
        ! is divided into two halves along the direction perpendicular to (fyx, fyy).

        chat1 = 0D0
        chat2 = 0D0

        lp00 = 0D0
        lp10 = 0D0
        lp01 = 0D0
        lp11 = 0D0
        lp20 = 0D0
        lp02 = 0D0
        ln00 = 0D0
        ln10 = 0D0
        ln01 = 0D0
        ln11 = 0D0
        ln20 = 0D0
        ln02 = 0D0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 hassperp = fhat * ttemp1 + ehat * ttemp2
                 aa = ker(ttemp1/ra, ttemp2/ra)

                 ! Start to fit two local linears with kernels.  

                 if (hassperp >= 0D0) then

                    lp00 = lp00 + aa
                    lp10 = lp10 + ttemp1 * aa
                    lp01 = lp01 + ttemp2 * aa
                    lp11 = lp11 + ttemp1 * ttemp2 * aa
                    lp20 = lp20 + ttemp1**2 * aa
                    lp02 = lp02 + ttemp2**2 * aa

                 else

                    ln00 = ln00 + aa
                    ln10 = ln10 + ttemp1 * aa
                    ln01 = ln01 + ttemp2 * aa
                    ln11 = ln11 + ttemp1 * ttemp2 * aa
                    ln20 = ln20 + ttemp1**2 * aa
                    ln02 = ln02 + ttemp2**2 * aa

                 end if
              end if

           end do
        end do

        detp = lp00 * lp20 * lp02 + lp10 * lp11 * lp01 + lp01 * &
             lp10 * lp11 - lp01 * lp20 * lp01 - lp10 * lp10 * &
             lp02 - lp00 * lp11 * lp11

        H1plus = lp11 * lp10 - lp01 * lp20
        H2plus = lp01 * lp10 - lp00 * lp11
        H3plus = lp00 * lp20 - lp10**2

        detn = ln00 * ln20 * ln02 + ln10 * ln11 * ln01 + ln01 * &
             ln10 * ln11 - ln01 * ln20 * ln01 - ln10 * ln10 * &
             ln02 - ln00 * ln11 * ln11

        H1minus = ln11 * ln10 - ln01 * ln20
        H2minus = ln01 * ln10 - ln00 * ln11
        H3minus = ln00 * ln20 - ln10**2

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 hassperp = fhat * ttemp1 + ehat * ttemp2
                 bb = ker(ttemp1/ra, ttemp2/ra)

                 if (hassperp >= 0D0) then

                    chat1 = chat1 + z1(i1, j1) * bb * &
                         (H1plus + H2plus * ttemp1 + H3plus &
                         * ttemp2)

                 else

                    chat2 = chat2 + z1(i1, j1) * bb * &
                         (H1minus + H2minus * ttemp1 + H3minus &
                         * ttemp2)

                 end if

              end if

           end do
        end do

        if (abs(detp) > 0D0) then

           chat1 = chat1/detp

        else

           chat1 = z1(i, j)

        end if

        if (abs(detn) > 0D0) then

           chat2 = chat2/detn

        else

           chat2 = z1(i, j)

        end if



        diff(i - k, j - k) = max(abs(bhat1 - bhat2), &
             abs(chat1 - chat2))

     end do
  end do

end subroutine roofDiff_denoise

