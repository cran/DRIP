!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! This is a Fortran subroutine to select bandwidth by computing         !
! bootstrap mse at edges and cv at continuity points. This selection    !
! procedure is proposed in the paper                                    !
! P. Qiu and Y. Kang "Blind Image Deblurring Using Jump Regression      !
! Analysis", Statistica Sinica.                                         !
!                                                                       !
! Creator: Yicheng Kang                                                 !
! Date: Sep 06 2015                                                     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deblur_3stage_bandwidth(n, obsImg, nband, bandwidth, edge1, edge2, &
     nboot, msecv)

  implicit none

  integer :: n, nband, bandwidth(1:nband)

  integer :: i, j, i1, j1, step(0:n, 0:n), i2, i3, k1, step1(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), &
       nboot, roof(0:n, 0:n), nstep, nroof, roof1(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), u, v, &
       nedge, bandw, edge1(0:n, 0:n), edge2(0:n, 0:n)
  
  double precision :: z(0:n, 0:n), z2(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), x, y, mse(1:nband), sigma, &
       fbhat(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), resid(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))),&
       temp, x1, y1, cv(1:nband), ra, w00, w10, w01, w11, w20, w02, ttemp1, ttemp2, ker, det, temp1, temp2, &
       temp3, fhat(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), temp11, sigmaxx, sigmayy, sigmaxy, del, &
       xbar, ybar, lambda1, prin, ttemp, ker1, temp22, zb(0:n, 0:n), zb2(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), &
       fboot(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), llkbw(1:7), aa, msecv(1:nband), optb1, obsImg(0:n, 0:n), &
       u_temp, v_temp

  external :: extend, extend1, ker, ker1, surfest, rndstart, rndend
  double precision :: myrunif
  external :: myrunif

  !! Read in data. Initialize fbhat and residuals.

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)
        step(i, j) = edge1(i, j)
        roof(i, j) = edge2(i, j)
        fbhat(i, j) = 0D0
        resid(i, j) = 0D0

     end do
  end do

  !! Use a preliminary local linear kernel smoothing to obtain residuals.

  sigma = 0D0
  optb1 = 0D0
  bandw = 7 ! The length of llkbw

  do i = 1, 7

     llkbw(i) = dble(i)/dble(n) !search bdwdth for conventnl LLK.

  end do

  call surfest(n, z(0:n, 0:n), bandw, llkbw, sigma, fbhat(0:n, 0:n), resid(0:n, 0:n), optb1)

  !! Edge detection.

  nedge = 0 ! Number of edges.

  do i = 0, n
     do j = 0, n

        nedge = nedge + step(i, j) + roof(i, j)

     end do
  end do

  !! Search bandwidth for deblurring. For each bandwidth, 
  !! compute  bootstrap mse at edges and cv at continuities.

  do i2 = 1, nband

     k1 = bandwidth(i2) ! seach bandwdith from bandwidth(1:nband).
     cv(i2) = 0D0
     mse(i2) = 0D0

     !! In a nbhd with bdwdth k1 of a given point, if the number of detected
     !! edge pixels is smaller than (k1-1)/2, then the surface is estimated as
     !! usual. Otherwise, fit a PC line through the detected edge pixels and then
     !! estimate surface locally.

     call extend(n, k1, z, z2(0:(n+2*k1), 0:(n+2*k1)))
     call extend1(n, k1, step, step1(0:(n+2*k1), 0:(n+2*k1)))
     call extend1(n, k1, roof, roof1(0:(n+2*k1), 0:(n+2*k1)))

     ra = dble(k1)/dble(n)

     !! These parameters are used in conventional surface estimation

     w00 = 0D0
     w10 = 0D0
     w01 = 0D0
     w11 = 0D0
     w20 = 0D0
     w02 = 0D0

     do i = -k1, k1
        do j = -k1, k1

           ttemp1 = dble(i)/dble(n)
           ttemp2 = dble(j)/dble(n)
           temp = ker(ttemp1/ra, ttemp2/ra)

           w00 = w00 + temp
           w10 = w10 + temp * ttemp1
           w01 = w01 + temp * ttemp2
           w11 = w11 + temp * ttemp1 * ttemp2
           w20 = w20 + temp * ttemp1**2
           w02 = w02 + temp * ttemp2**2

        end do
     end do

     det = w00 * w20 * w02 + w10 * w11 * w01 + w01 * w10 * w11 - &
          w01 * w20 * w01 - w10 * w10 * w02 - w00 * w11 * w11

     temp1 = w20 * w02 - w11**2
     temp2 = w01 * w11 - w10 * w02
     temp3 = w10 * w11 - w01 * w20

     !! Start estimating surface

     do i = k1, n + k1
        do j = k1, n + k1

           x = dble(i - k1)/dble(n)
           y = dble(j - k1)/dble(n)

           fhat(i, j)=0D0
           temp11 = 0D0

           !! nstep/nroof denotes the total number of step/roof edge pixels in the neighborhood
           !! If both nstep and nroof <= (k1-1)/2, then we estimate the surface by conventional 
           !! local linear kernel estimator. Otherwise, we do deblurring.

           nstep = 0
           nroof = 0

           do i1 =i - k1, i + k1
              do j1 = j - k1, j + k1

                 x1 = dble(i1 - k1)/dble(n)
                 y1 = dble(j1 - k1)/dble(n)

                 if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                      (step1(i1, j1) .EQ. 1)) then

                    nstep = nstep + 1

                 end if

                 if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                      (roof1(i1, j1) .EQ. 1)) then

                    nroof = nroof + 1

                 end if

              end do
           end do

           if (nstep .LE. int(k1/2) .and. nroof <= int(k1/2)) then

              do i1 = i - k1, i + k1
                 do j1 = j - k1, j + k1

                    x1 = dble(i1 - k1)/dble(n)
                    y1 = dble(j1 - k1)/dble(n)

                    ttemp1 = x1 - x
                    ttemp2 = y1 - y
                    temp = ker(ttemp1/ra, ttemp2/ra)

                    fhat(i, j) = fhat(i, j) + z2(i1, j1) * temp * &
                         (temp1 + temp2 * ttemp1 + temp3 * ttemp2)

                 end do
              end do

              if (abs(det) > 0D0) then

                 fhat(i, j) = fhat(i, j)/det

              else

                 fhat(i, j) = z2(i, j)

              end if

              !! End for conventional smoothing.

           else

              if (nstep > int(k1/2)) then

                 !! Fit a principal component line through the detected step edge pixels
                 !! and ignore roof edges.

                 sigmaxx = 0D0
                 sigmayy = 0D0
                 sigmaxy = 0D0
                 xbar = 0D0
                 ybar = 0D0

                 do i1 = i - k1, i + k1
                    do j1 = j - k1, j + k1

                       x1 = dble(i1 - k1)/dble(n)
                       y1 = dble(j1 - k1)/dble(n)

                       if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                            (step1(i1, j1) .EQ. 1)) then

                          xbar = xbar + x1
                          ybar = ybar + y1
                          sigmaxx = sigmaxx + x1 * x1
                          sigmayy = sigmayy + y1 * y1
                          sigmaxy = sigmaxy + x1 * y1

                       end if

                    end do
                 end do

                 xbar=xbar/dble(nstep)
                 ybar=ybar/dble(nstep)
                 sigmaxx = sigmaxx/dble(nstep) - xbar**2
                 sigmayy = sigmayy/dble(nstep) - ybar**2
                 sigmaxy = sigmaxy/dble(nstep) - xbar * ybar

                 del = (sigmaxx - sigmayy)**2 + 4D0 * sigmaxy**2
                 lambda1 = (sigmayy - sigmaxx - sqrt(del))/2D0

                 !!Two cases for principal component line.

                 if (sigmaxx > 0D0) then

                    prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)
                    ttemp1 = abs(prin)/sqrt(sigmaxy**2 + lambda1**2) + ra

                    if (prin .GE. 0D0) then

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  &
                                  (temp22 .GE. 0D0)) then 

                                ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i, j)

                       end if

                    else

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  & 
                                  (temp22 .LT. 0D0)) then 

                                ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa


                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i,j)

                       end if

                    end if

                 else !In this case, principal component line x=xbar

                    prin = x - xbar
                    ttemp1 = abs(prin) + ra

                    if (prin .GE. 0D0) then

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = x1 - xbar
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  &
                                  (temp22 .GE. 0D0)) then 

                                ttemp2 = abs(temp22)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i, j)

                       end if

                    else

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = x1 - xbar
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  & 
                                  (temp22 .LT. 0D0)) then 

                                ttemp2 = abs(temp22)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i,j)

                       end if

                    end if

                 end if

                 !! End for deblurring around step edges.       

              else

                 !! Fit a principal component line through detected roof edge pixels. 

                 sigmaxx = 0D0
                 sigmayy = 0D0
                 sigmaxy = 0D0
                 xbar = 0D0
                 ybar = 0D0

                 do i1 = i - k1, i + k1
                    do j1 = j - k1, j + k1

                       x1 = dble(i1 - k1)/dble(n)
                       y1 = dble(j1 - k1)/dble(n)

                       if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                            (roof1(i1, j1) .EQ. 1)) then

                          xbar = xbar + x1
                          ybar = ybar + y1
                          sigmaxx = sigmaxx + x1 * x1
                          sigmayy = sigmayy + y1 * y1
                          sigmaxy = sigmaxy + x1 * y1

                       end if

                    end do
                 end do

                 xbar=xbar/dble(nroof)
                 ybar=ybar/dble(nroof)
                 sigmaxx = sigmaxx/dble(nroof) - xbar**2
                 sigmayy = sigmayy/dble(nroof) - ybar**2
                 sigmaxy = sigmaxy/dble(nroof) - xbar * ybar

                 del = (sigmaxx - sigmayy)**2 + 4D0 * sigmaxy**2
                 lambda1 = (sigmayy - sigmaxx - sqrt(del))/2D0

                 !! Two cases for principal component line.

                 if (sigmaxx > 0D0) then

                    prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)
                    ttemp1 = abs(prin)/sqrt(sigmaxy**2 + lambda1**2) + ra

                    if (prin .GE. 0D0) then

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  &
                                  (temp22 .GE. 0D0)) then 

                                ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i, j)

                       end if

                    else

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  & 
                                  (temp22 .LT. 0D0)) then 

                                ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i,j)

                       end if

                    end if

                 else ! In this case, the PC line through roof edge pixels x=xbar

                    prin = x - xbar
                    ttemp1 = abs(prin) + ra

                    if (prin .GE. 0D0) then

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = x1 - xbar
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  &
                                  (temp22 .GE. 0D0)) then 

                                ttemp2 = abs(temp22)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i, j)

                       end if

                    else

                       do i1 = i - k1, i + k1
                          do j1 = j - k1, j + k1

                             x1 = dble(i1 - k1)/dble(n)
                             y1 = dble(j1 - k1)/dble(n)

                             temp22 = x1 - xbar
                             ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                             if ((ttemp .LE. ra) .AND.  & 
                                  (temp22 .LT. 0D0)) then 

                                ttemp2 = abs(temp22)
                                aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                     ker1(ttemp2/ttemp1)
                                fhat(i, j) = fhat(i, j) + z2(i1, j1) * aa
                                temp11 = temp11 + aa

                             end if

                          end do
                       end do

                       if (abs(temp11) > 0D0) then

                          fhat(i, j) = fhat(i, j)/temp11

                       else

                          fhat(i, j) = z2(i,j)

                       end if

                    end if

                 end if

                 !! End of deblurring around roof edges.

              end if

           end if

           ! End of estimating f.

        end do
     end do

     !! Repeat deblurring procedure with bootstrap data.
     !! Set seed for random number generator.

     call rndstart()

     do i3 = 1, nboot

        !! Obtain bootstrap sample.

        do i = 0, n
           do j = 0, n

              !call random_number(u_temp)
              !call random_number(v_temp)
              u_temp = myrunif(0D0, 1D0)
              v_temp = myrunif(0D0, 1D0)
              u = int(u_temp * dble(n+1))
              v = int(v_temp * dble(n+1))
              !u = int(dble(rand(0)) * dble(n + 1))
              !v = int(dble(rand(0)) * dble(n + 1)) 
              zb(i, j) = fbhat(i, j) + resid(u, v)

           end do
        end do

        !! Extend to avoid boundary problems. Bootstrap sample.

        call extend(n, k1, zb, zb2(0:(n+2*k1), 0:(n+2*k1)))

        !! Start estimating surface. Bootstrap sample

        do i = k1, n + k1
           do j = k1, n + k1

              x = dble(i - k1)/dble(n)
              y = dble(j - k1)/dble(n)

              fboot(i, j)=0D0
              temp11 = 0D0

              !! nstep/nroof denotes the total number of step/roof edge pixels in the neighborhood
              !! If both nstep and nroof <= (k1-1)/2, then we estimate the surface by conventional 
              !! local linear kernel estimator. Otherwise, we do deblurring. Bootstrap sample.

              nstep = 0
              nroof = 0

              do i1 =i - k1, i + k1
                 do j1 = j - k1, j + k1

                    x1 = dble(i1 - k1)/dble(n)
                    y1 = dble(j1 - k1)/dble(n)

                    if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                         (step1(i1, j1) .EQ. 1)) then

                       nstep = nstep + 1

                    end if

                    if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                         (roof1(i1, j1) .EQ. 1)) then

                       nroof = nroof + 1

                    end if
                 end do
              end do

              if (nstep .LE. int(k1/2) .and. nroof <= int(k1/2)) then

                 do i1 = i - k1, i + k1
                    do j1 = j - k1, j + k1

                       x1 = dble(i1 - k1)/dble(n)
                       y1 = dble(j1 - k1)/dble(n)

                       ttemp1 = x1 - x
                       ttemp2 = y1 - y
                       temp = ker(ttemp1/ra, ttemp2/ra)

                       fboot(i, j) = fboot(i, j) + zb2(i1, j1) * temp * &
                            (temp1 + temp2 * ttemp1 + temp3 * ttemp2)

                    end do
                 end do

                 if (abs(det) > 0D0) then

                    fboot(i, j) = fboot(i, j)/det

                 else

                    fboot(i, j) = zb2(i, j)

                 end if

                 !! End for conventional smoothing. Bootstrap sample.

              else

                 if (nstep > int(k1/2)) then

                    !! Fit a principal component line through the detected step edge pixels
                    !! and ignore roof edges. Bootstrap sample.

                    sigmaxx = 0D0
                    sigmayy = 0D0
                    sigmaxy = 0D0
                    xbar = 0D0
                    ybar = 0D0

                    do i1 = i - k1, i + k1
                       do j1 = j - k1, j + k1

                          x1 = dble(i1 - k1)/dble(n)
                          y1 = dble(j1 - k1)/dble(n)

                          if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                               (step1(i1, j1) .EQ. 1)) then

                             xbar = xbar + x1
                             ybar = ybar + y1
                             sigmaxx = sigmaxx + x1 * x1
                             sigmayy = sigmayy + y1 * y1
                             sigmaxy = sigmaxy + x1 * y1

                          end if

                       end do
                    end do

                    xbar=xbar/dble(nstep)
                    ybar=ybar/dble(nstep)
                    sigmaxx = sigmaxx/dble(nstep) - xbar**2
                    sigmayy = sigmayy/dble(nstep) - ybar**2
                    sigmaxy = sigmaxy/dble(nstep) - xbar * ybar

                    del = (sigmaxx - sigmayy)**2 + 4D0 * sigmaxy**2
                    lambda1 = (sigmayy - sigmaxx - sqrt(del))/2D0
                    
                    !Two cases for principal component line.

                    if (sigmaxx > 0D0) then

                       prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)
                       ttemp1 = abs(prin)/sqrt(sigmaxy**2 + lambda1**2) + ra

                       if (prin .GE. 0D0) then

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  &
                                     (temp22 .GE. 0D0)) then 

                                   ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i, j)

                          end if

                       else

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  & 
                                     (temp22 .LT. 0D0)) then 

                                   ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa


                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i,j)

                          end if

                       end if

                    else !In this case, principal component line x=xbar

                       prin = x - xbar
                       ttemp1 = abs(prin) + ra

                       if (prin .GE. 0D0) then

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = x1 - xbar
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  &
                                     (temp22 .GE. 0D0)) then 

                                   ttemp2 = abs(temp22)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i, j)

                          end if

                       else

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = x1 - xbar
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  & 
                                     (temp22 .LT. 0D0)) then 

                                   ttemp2 = abs(temp22)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i,j)

                          end if

                       end if

                    end if

                    !! End for deblurring around step edges. Bootstrap sample     

                 else

                    !! Fit a principal component line through detected roof edge pixels. 
                    !! Bootstrap sample.

                    sigmaxx = 0D0
                    sigmayy = 0D0
                    sigmaxy = 0D0
                    xbar = 0D0
                    ybar = 0D0

                    do i1 = i - k1, i + k1
                       do j1 = j - k1, j + k1

                          x1 = dble(i1 - k1)/dble(n)
                          y1 = dble(j1 - k1)/dble(n)

                          if (((x1 - x)**2 + (y1 - y)**2 .LE. ra**2) .AND. &
                               (roof1(i1, j1) .EQ. 1)) then

                             xbar = xbar + x1
                             ybar = ybar + y1
                             sigmaxx = sigmaxx + x1 * x1
                             sigmayy = sigmayy + y1 * y1
                             sigmaxy = sigmaxy + x1 * y1

                          end if

                       end do
                    end do

                    xbar=xbar/dble(nroof)
                    ybar=ybar/dble(nroof)
                    sigmaxx = sigmaxx/dble(nroof) - xbar**2
                    sigmayy = sigmayy/dble(nroof) - ybar**2
                    sigmaxy = sigmaxy/dble(nroof) - xbar * ybar

                    del = (sigmaxx - sigmayy)**2 + 4D0 * sigmaxy**2
                    lambda1 = (sigmayy - sigmaxx - sqrt(del))/2D0

                    !!Two cases for principal component line.

                    if (sigmaxx > 0D0) then

                       prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)
                       ttemp1 = abs(prin)/sqrt(sigmaxy**2 + lambda1**2) + ra

                       if (prin .GE. 0D0) then

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  &
                                     (temp22 .GE. 0D0)) then 

                                   ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i, j)

                          end if

                       else

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = sigmaxy * (x1 - xbar) + lambda1 * (y1 - ybar)
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  & 
                                     (temp22 .LT. 0D0)) then 

                                   ttemp2 = abs(temp22)/sqrt(sigmaxy**2 + lambda1**2)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i,j)

                          end if

                       end if

                    else !! In this case, the PC line through roof edge pixels x=xbar

                       prin = x - xbar
                       ttemp1 = abs(prin) + ra

                       if (prin .GE. 0D0) then

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = x1 - xbar
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  &
                                     (temp22 .GE. 0D0)) then 

                                   ttemp2 = abs(temp22)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i, j)

                          end if

                       else

                          do i1 = i - k1, i + k1
                             do j1 = j - k1, j + k1

                                x1 = dble(i1 - k1)/dble(n)
                                y1 = dble(j1 - k1)/dble(n)

                                temp22 = x1 - xbar
                                ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                                if ((ttemp .LE. ra) .AND.  & 
                                     (temp22 .LT. 0D0)) then 

                                   ttemp2 = abs(temp22)
                                   aa = ker(dble(i1 - i)/dble(k1), dble(j1 - j)/dble(k1)) * &
                                        ker1(ttemp2/ttemp1)
                                   fboot(i, j) = fboot(i, j) + zb2(i1, j1) * aa
                                   temp11 = temp11 + aa

                                end if

                             end do
                          end do

                          if (abs(temp11) > 0D0) then

                             fboot(i, j) = fboot(i, j)/temp11

                          else

                             fboot(i, j) = zb2(i,j)

                          end if

                       end if

                    end if

                    !! End of deblurring around roof edges.
                    !! Booststrap sample.

                 end if

              end if

              !! End of estimating fboot.

           end do
        end do

        !! Calculate msecv.

        do i = k1, k1 + n
           do j = k1, k1 + n

              if (step1(i, j) + roof1(i, j) > 0) then

                 mse(i2) = mse(i2) + (fboot(i, j) - fhat(i, j))**2

              else

                 cv(i2) = cv(i2) + (fboot(i, j) - zb2(i, j))**2

              end if

           end do
        end do

     end do !! End boostrap loop.

     call rndend()

     mse(i2) = mse(i2)/dble(nedge)/dble(nboot)
     cv(i2) = cv(i2)/dble((n + 1)**2 - nedge)/dble(nboot)
     msecv(i2) = mse(i2) + cv(i2)


  end do !! End of searching bandwidth loop.


end subroutine deblur_3stage_bandwidth
