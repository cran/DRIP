!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                            !
! This is a Fortran subroutine to DEBLUR an image using local!
! pixel clustering and nonparametric regression. This method !
! is proposed in the paper                                   !
! Kang, Y. and Qiu, P. 'Efficient Blind Image Deblurring     !
! Using Nonparametric Regression and Local Pixel Clustering'.!
! Creator: Yicheng Kang                                      !
! Date: April 30, 2013                                       !
!                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cluster_deblur(n, obsImg, bandwidth, mmt2, mmt4, &
  zq, estImg)

  implicit none

  double precision :: z(0:600, 0:600), mmt2, small, varnum, &
       z1(0:600, 0:600), ra, w00, w20, ttemp1, ttemp2, ker, temp, &
       ahat, bhat, chat, obsImg(0:n, 0:n), e, zmin, zmax, zbar, &
       z1bar, z2bar, varden, sd, z1ss, z2ss, bwr(1:1000), &
       bwrLocalMax, cutoff, ttemp, ker1, fhat, zascend(1:1000), &
       mmt4, estImg(0:n, 0:n), zq, tolerance

  integer :: n, i, j, k, i1, j1, s, n1clust, n2clust, sl, su, &
       nitem, l1, l2, locsm, nnghb, bandwidth

  external :: extend, ker, ker1

  !! Assign value to parameters.

  tolerance = 1D-8

  !! Read in data.

  do i = 0, n
     do j = 0, n

        z(i, j) = obsImg(i, j)

     end do
  end do

  !! Deblur bandwidth

  k = bandwidth
  ra = dble(k)/dble(n)

  !! Extend to avoid boundary problems.

  call extend(n, k, z, z1)

  !! The following quantities are used in the conventional
  !! local linear CIRCULARLY SYMMETRIC kernel smoothing in
  !! a CIRCULAR neighborhood.

  w00 = 0D0
  w20 = 0D0
  varnum = 0D0
  varden = 0D0

  do i1 = -k, k
     do j1 = -k, k

        if (i1**2 + j1**2 <= k**2) then

           ttemp1 = dble(i1)/dble(n)
           ttemp2 = dble(j1)/dble(n)
           temp = ker(ttemp1/ra, ttemp2/ra)
           w00 = w00 + temp
           w20 = w20 + temp * ttemp1**2
           varnum = varnum + temp**2
           varden = varden + temp

        end if

     end do
  end do

  varnum = sqrt(varnum)
  sd = varnum/varden * sqrt(mmt4 - mmt2**2)

  !! Start estimating surface

  do i = k, n + k
     do j = k, n + k

        fhat = 0D0


!!!! Fit local plane using 0%-trimmed observations.

        ahat = 0D0
        bhat = 0D0
        chat = 0D0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 temp = ker(ttemp1/ra, ttemp2/ra)
                 ahat = ahat + z1(i1, j1) * temp 
                 bhat = bhat + z1(i1, j1) * ttemp1 * temp
                 chat = chat + z1(i1, j1) * ttemp2 * temp

              end if

           end do
        end do

        ahat = ahat/w00
        bhat = bhat/w20
        chat = chat/w20

!!!! Compute the local WRMS error for design point (i, j).
!!!! The local WRMS error is an estimate for sigma in the
!!!! neighborhood. We will use it to decide whether there
!!!! are more than one cluster in the neighborhood.

        e = 0D0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                 ttemp1 = dble(i1 - i)/dble(n)
                 ttemp2 = dble(j1 - j)/dble(n)
                 temp = ker(ttemp1/ra, ttemp2/ra)
                 e = e + temp * & 
                      (z1(i1, j1) - ahat - bhat * ttemp1 &
                      - chat * ttemp2)**2

              end if

           end do
        end do


        e = e/w00

!!!! If (e(i, j)- sigma^2) < threshold, then we estimate the surface by
!!!! the conventional LLK smoothing. Else, we do clustering.

        if (abs(e - mmt2) <= zq * sd) then

           fhat = ahat

        else

!!!! Order the observations in the neighborhood.

           nitem = 0

           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                    nitem = nitem + 1
                    zascend(nitem) = z1(i1, j1)

                 end if

              end do
           end do

           do l1 = 1, nitem - 1

              small = zascend(l1)
              locsm = l1

              do l2 = l1 + 1, nitem

                 if (zascend(l2) < small) then

                    small = zascend(l2)
                    locsm = l2

                 end if

              end do

              zascend(locsm) = zascend(l1)
              zascend(l1) = small

           end do

           zmin = zascend(1)
           zmax = zascend(nitem)


!!!! Search for the best cut-off constant.

           sl = int(dble(nitem) * 0.25D0)
           su = int(dble(nitem) * 0.75D0)

           zbar = 0D0
           z1bar = 0D0
           z1ss = 0D0
           z2bar = 0D0
           z2ss = 0D0

!!!! Initial cut.

           do  l1 = 1, nitem

              zbar = zbar + zascend(l1)

              if (l1 < sl) then

                 z1bar = z1bar + zascend(l1)
                 z1ss = z1ss + zascend(l1)**2

              else

                 z2bar = z2bar + zascend(l1)
                 z2ss = z2ss + zascend(l1)**2

              end if

           end do

           zbar = zbar/dble(nitem)
           n1clust = sl - 1
           z1bar = z1bar/dble(n1clust)
           n2clust = nitem - sl + 1
           z2bar = z2bar/dble(n2clust)

!!!! Keep cutting.

           do s = sl, su

              z1bar = dble(n1clust) * z1bar + zascend(s)
              z1ss = z1ss + zascend(s)**2
              z2bar = dble(n2clust) * z2bar - zascend(s)
              z2ss = z2ss - zascend(s)**2
              n1clust = n1clust + 1
              n2clust = n2clust - 1
              z1bar = z1bar/dble(n1clust)
              z2bar = z2bar/dble(n2clust)


!!!! Comput between-group-within-group ratio(bwr).

              bwr(s) = (dble(n1clust) * (z1bar - zbar)**2 + &
                   dble(n2clust) * (z2bar - zbar)**2) / &
                   (z1ss - dble(n1clust) * z1bar**2 + z2ss - &
                   dble(n2clust) * z2bar**2)

           end do

!!!! Maximize bwr.

           bwrLocalMax = maxval(bwr(sl:su))


!!!! Find the optimal cut-off constant.

           cutoff = 0D0

           do s = sl, su

              if (abs(bwr(s) - bwrLocalMax) < tolerance) then

                 cutoff = zascend(s)

              end if

           end do

           !! Start local smoothing within the cluster.

           temp = 0D0

!!!! If z1(i, j) <= cutoff.

           if (z1(i, j) <= cutoff) then

!!!! Inspect whether neighboring pixels are also <= cutoff.

              nnghb = 0

              do i1 = i - 1, i + 1
                 do j1 = j - 1, j + 1

                    if ((z1(i1, j1) <= cutoff) .and. &
                         (abs(i1 - i) + abs(j1 - j) > 0)) then

                       nnghb = nnghb + 1

                    end if

                 end do
              end do

!!!! (i, j) is in cluster 1.

              if (nnghb >= 4) then

                 do i1 = i - k, i + k
                    do j1 = j - k, j + k

                       if (z1(i1, j1) <= cutoff) then

                          ttemp1 = dble(i1 - i)/dble(k)
                          ttemp2 = dble(j1 - j)/dble(k)
                          ttemp =  ker(ttemp1, ttemp2) * & 
                               ker1(abs(z1(i1, j1) - cutoff)/&
                               abs(zmin - cutoff))
                          fhat = fhat + z1(i1, j1) * ttemp
                          temp = temp + ttemp

                       end if

                    end do
                 end do

                 fhat = fhat/temp

              else

!!!! (i, j) is in cluster 2.

                 do i1 = i - k, i + k
                    do j1 = j - k, j + k

                       if (z1(i1, j1) > cutoff) then

                          ttemp1 = dble(i1 - i)/dble(k)
                          ttemp2 = dble(j1 - j)/dble(k)
                          ttemp =  ker(ttemp1, ttemp2) * & 
                               ker1(abs(z1(i1, j1) - cutoff)/&
                               abs(zmax - cutoff))
                          fhat = fhat + z1(i1, j1) * ttemp
                          temp = temp + ttemp

                       end if

                    end do
                 end do

                 fhat = fhat/temp

              end if


!!!! If z1(i, j) > cutoff.

           else

!!!! Inspect whether neighboring pixels are also > cutoff.

              nnghb = 0

              do i1 = i - 1, i + 1
                 do j1 = j - 1, j + 1

                    if ((z1(i1, j1) > cutoff) .and. &
                         (abs(i1 - i) + abs(j1 - j) > 0)) then

                       nnghb = nnghb + 1

                    end if

                 end do
              end do

!!!! (i, j) is in cluster 2.

              if (nnghb >= 4) then

                 do i1 = i - k, i + k
                    do j1 = j - k, j + k

                       if (z1(i1, j1) > cutoff) then

                          ttemp1 = dble(i1 - i)/dble(k)
                          ttemp2 = dble(j1 - j)/dble(k)
                          ttemp =  ker(ttemp1, ttemp2) * & 
                               ker1(abs(z1(i1, j1) - cutoff)/&
                               abs(zmax - cutoff))
                          fhat = fhat + z1(i1, j1) * ttemp
                          temp = temp + ttemp

                       end if

                    end do
                 end do

                 fhat = fhat/temp

              else

!!!! (i, j) is in cluster 1.

                 do i1 = i - k, i + k
                    do j1 = j - k, j + k

                       if (z1(i1, j1) <= cutoff) then

                          ttemp1 = dble(i1 - i)/dble(k)
                          ttemp2 = dble(j1 - j)/dble(k)
                          ttemp =  ker(ttemp1, ttemp2) * & 
                               ker1(abs(z1(i1, j1) - cutoff)/&
                               abs(zmin - cutoff))
                          fhat = fhat + z1(i1, j1) * ttemp
                          temp = temp + ttemp

                       end if

                    end do
                 end do

                 fhat = fhat/temp

              end if

           end if

        end if

        estImg(i - k, j - k) = fhat

     end do
  end do



end subroutine cluster_deblur


