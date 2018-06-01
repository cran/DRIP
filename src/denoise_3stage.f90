!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                            !
! This is a Fortran souce file to recover jump surfaces from !
! noisy data using the method proposed in the paper          !
! Kang, Y. and Qiu, P, 'Blind Image Deblurring Using Jump    !
!                       Regression Analysis'.                !
! Creator: Yicheng Kang                                      !
! Date: May 7, 2013                                          !
!                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine denoise_3stage(n, obsImg, bandwidth, edge1, edge2, &
     estImg)

  implicit none

  integer :: n, bandwidth

  integer :: i, j, i1, j1, edge1(0:n, 0:n), k, nroof, nstep, &
       step1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), edge2(0:n, 0:n), &
       roof1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), step(0:n, 0:n), roof(0:n, 0:n)

  double precision :: z(0:n, 0:n), z1(0:(n+2*bandwidth), 0:(n+2*bandwidth)), x, y, &
       temp, x1, y1, ra, w00, ttemp1, ttemp2, ker, fhat, del, &
       temp11, sigmaxx, sigmayy, prin, ttemp, temp22, aa, xbar, &
       ybar, lambda1, sigmaxy, obsImg(0:n, 0:n), estImg(0:n, 0:n)

  external :: extend, extend1, ker

  ! Assign values to parameters.

  k = bandwidth
  ra = dble(k)/dble(n)

  ! Read in data and the detected edges.

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)
        step(i, j) = edge1(i, j)
        roof(i, j) = edge2(i, j)

     end do
  end do

  ! Extend to avoid boundary problems.

  call extend(n, k, z, z1)
  call extend1(n, k, step, step1)
  call extend1(n, k, roof, roof1)

  ! In a nbhd with bdwdth k of a given point, if the number
  ! of detected edge pixels is smaller than (k-1)/2, then the
  ! surface is estimated as usual. Otherwise, fit a PC line  
  ! through the detected edge pixels and then estimate surface 
  ! locally.


  ! These parameters are used in conventional surface estimation

  w00 = 0D0

  do i = -k, k
     do j = -k, k

        if (i**2 + j**2 <= k**2) then

           ttemp1 = dble(i)/dble(n)
           ttemp2 = dble(j)/dble(n)
           temp = ker(ttemp1/ra, ttemp2/ra)
           w00 = w00 + temp

        end if

     end do
  end do

  !! Start estimating surface

  do i = k, n + k
     do j = k, n + k

        fhat = 0D0
        temp11 = 0D0
        x = dble(i - k)/dble(n)
        y = dble(j - k)/dble(n)

        ! nstep/nroof denotes the total number of step/roof edge pixels 
        ! in the neighborhood. If both nstep and nroof <= (k-1)/2, then 
        ! we estimate the surface by conventional local linear kernel 
        ! estimator. Otherwise, we preserve edge when denoising.

        nstep = 0
        nroof = 0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 if (step1(i1, j1) == 1) then

                    nstep = nstep + 1

                 end if

                 if (roof1(i1, j1) == 1) then

                    nroof = nroof + 1

                 end if

              end if

           end do
        end do

        if (nstep <= int(k/2) .and. nroof <= int(k/2)) then

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                    ttemp1 = dble(i1 - i)/dble(n)
                    ttemp2 = dble(j1 - j)/dble(n)
                    temp = ker(ttemp1/ra, ttemp2/ra)
                    fhat = fhat + z1(i1, j1) * temp

                 end if

              end do
           end do

           if (abs(w00) > 0D0) then

              fhat = fhat/w00

           else

              fhat = z1(i, j)

           end if

           ! End for conventional smoothing.

        else

           if (nstep > int(k/2)) then

              ! Fit a principal component line through the 
              ! detected step edge pixels and ignore roof edges.

              sigmaxx = 0D0
              sigmayy = 0D0
              sigmaxy = 0D0
              xbar = 0D0
              ybar = 0D0

              do i1 = i - k, i + k
                 do j1 = j - k, j + k

                    if (((i1 - i)**2 + (j1 - j)**2 <= k**2) .AND. &
                         (step1(i1, j1) == 1)) then

                       x1 = dble(i1 - k)/dble(n)
                       y1 = dble(j1 - k)/dble(n)
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

              ! Two cases for principal component line.

              if (sigmaxx > 0D0) then

                 prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)

                 if (prin >= 0D0) then

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = sigmaxy * (x1 - xbar) + lambda1 &
                               * (y1 - ybar)
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 >= 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i, j)

                    end if

                 else

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = sigmaxy * (x1 - xbar) + lambda1 &
                               * (y1 - ybar)
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 < 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i,j)

                    end if

                 end if

              else !In this case, principal component line x=xbar

                 prin = x - xbar

                 if (prin >= 0D0) then

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = x1 - xbar
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 >= 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i, j)

                    end if

                 else

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = x1 - xbar
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 < 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k)) 
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i,j)

                    end if

                 end if
              end if

              ! End for deblurring around step edges.       

           else

              ! Fit a principal component line through detected
              ! roof edge pixels. 

              sigmaxx = 0D0
              sigmayy = 0D0
              sigmaxy = 0D0
              xbar = 0D0
              ybar = 0D0

              do i1 = i - k, i + k
                 do j1 = j - k, j + k

                    if (((i1 - i)**2 + (j1 - j)**2 <= k**2) .AND. &
                         (roof1(i1, j1) == 1)) then

                       x1 = dble(i1 - k)/dble(n)
                       y1 = dble(j1 - k)/dble(n)
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

              ! Two cases for principal component line.

              if (sigmaxx > 0D0) then

                 prin = sigmaxy * (x - xbar) + lambda1 * (y - ybar)

                 if (prin >= 0D0) then

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = sigmaxy * (x1 - xbar) + lambda1 &
                               * (y1 - ybar)
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 >= 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i, j)

                    end if

                 else

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = sigmaxy * (x1 - xbar) + lambda1 &
                               * (y1 - ybar)
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 < 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i,j)

                    end if

                 end if

              else 

                 !In this case, the PC line through roof edge pixels x=xbar

                 prin = x - xbar

                 if (prin >= 0D0) then

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = x1 - xbar
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 >= 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i, j)

                    end if

                 else

                    do i1 = i - k, i + k
                       do j1 = j - k, j + k

                          x1 = dble(i1 - k)/dble(n)
                          y1 = dble(j1 - k)/dble(n)

                          temp22 = x1 - xbar
                          ttemp = sqrt((x1 - x)**2 + (y1 - y)**2)

                          if ((ttemp <= ra) .AND. (temp22 < 0D0)) then 

                             aa = ker(dble(i1 - i)/dble(k), &
                                  dble(j1 - j)/dble(k))
                             fhat = fhat + z1(i1, j1) * aa
                             temp11 = temp11 + aa

                          end if

                       end do
                    end do

                    if (abs(temp11) > 0D0) then

                       fhat = fhat/temp11

                    else

                       fhat = z1(i,j)

                    end if

                 end if
              end if

              ! End of deblurring around roof edges.

           end if

        end if

        ! End of estimating f.

        estImg(i - k, j - k) = fhat

     end do
  end do

end subroutine denoise_3stage

