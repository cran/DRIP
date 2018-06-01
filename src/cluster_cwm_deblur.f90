! This is a Fortran 90 subroutine to deblur an image using local clustering of center weighted median
! filter output. Bilateral piece-wise constant method is used in the local clustering. Bandwidth is given 
! by users.

  Subroutine cluster_cwm_deblur(n, obsImg, k, zq, sigma, phi0, mean_std_abs, cw, estImg)

    implicit none

    integer :: n, i, j, k, i1, j1, s, n1clust, n2clust, l1, sl, su, nitem, &
         cw, ind(1:((2*k+1)**2))

    double precision :: obsImg(0:n, 0:n), sigma, z1(0:(n+2*k), 0:(n+2*k)), &
         ra, w00, ahat,  fhat(0:(n+2*k), 0:(n+2*k)), ttemp1, ttemp2, ker, temp, &
         e, zmin, zmax, zbar,z1bar, z2bar, zascend(1:((2*k+1)**2)), z1ss, z2ss, &
         bwr(1:((2*k+1)**2)), cutoff, ttemp, zq, ker1, estImg(0:n, 0:n), phi0, &
         asy_thresh, ftilde(0:n, 0:n), pi, ftilde1(0:(n+2*k), 0:(n+2*k)), mean_std_abs


    external :: extend, ker, ker1, LocalMedianFilter, qsortd


    !! Assign values to parameters.

    pi = 3.14159265D0

    !! Caculate bandwidth and asymptotic threshold.

    ra = dble(k)/dble(n)
    !asy_thresh = zq * sigma / dble(k) * sqrt(0.5D0 - 2D0/(3D0*pi))
    asy_thresh = zq / dble(k) * sigma * sqrt(4D0/(3D0*pi) + 1D0/(4D0*pi*phi0**2) - mean_std_abs/(pi*phi0))

    !! Apply local median filtering.

    call LocalMedianFilter(n, k, cw, obsImg(0:n, 0:n), ftilde(0:n, 0:n))
    call extend(n, k, ftilde(0:n, 0:n), ftilde1(0:(n+2*k), 0:(n+2*k)))


    !! Extend to avoid boundary problems.

    call extend(n, k, obsImg(0:n, 0:n), z1(0:(n+2*k), 0:(n+2*k)))

    !! The following quantities are used in the conventional
    !! local linear CIRCULARLY SYMMETRIC kernel smoothing in
    !! a CIRCULAR neighborhood.

    w00 = 0D0

    do i1 = -k, k
       do j1 = -k, k

          if (i1**2 + j1**2 <= k**2) then

             ttemp1 = dble(i1)/dble(n)
             ttemp2 = dble(j1)/dble(n)
             temp = ker(ttemp1/ra, ttemp2/ra)
             w00 = w00 + temp

          end if

       end do
    end do


    !! Start estimating surface

      do i = k, n + k
         do j = k, n + k

         fhat(i, j) = 0D0


         !! Fit local plane using 0%-trimmed observations.

         ahat = 0D0

         do i1 = i - k, i + k
            do j1 = j - k, j + k


                if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                   ttemp1 = dble(i1 - i)/dble(n)
                   ttemp2 = dble(j1 - j)/dble(n)
                   temp = ker(ttemp1/ra, ttemp2/ra)
                   ahat = ahat + z1(i1, j1) * temp 

                end if

             end do
          end do

          ahat = ahat/w00

          !! Compare LMF output with LLK fit. We will use it to decide whether (i, j) is an edge point.

          e = abs(ahat - ftilde1(i, j))

          if (e >= asy_thresh) then

             !! Order the LMF output in the neighborhood.

             nitem = 0

             do i1 = i - k, i + k
                do j1 = j - k, j + k

                   if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then
                     
                      nitem = nitem + 1
                      zascend(nitem) = ftilde1(i1, j1)

                   end if

                end do
             end do

             call qsortd(zascend(1:nitem), ind(1:nitem), nitem)

             zascend(1:nitem) = zascend(ind(1:nitem))
             zmin = zascend(1)
             zmax = zascend(nitem)

             !! Search for the best cut-off constant.

             sl = int(dble(nitem) * 0.1D0)
             su = int(dble(nitem) * 0.9D0)

             zbar = 0D0
             z1bar = 0D0
             z1ss = 0D0
             z2bar = 0D0
             z2ss = 0D0
             bwr(1:nitem) = 0D0

             !! Initial cut.

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

             !! Keep cutting.

             do s = sl, su

                z1bar = dble(n1clust) * z1bar + zascend(s)
                z1ss = z1ss + zascend(s)**2
                z2bar = dble(n2clust) * z2bar - zascend(s)
                z2ss = z2ss - zascend(s)**2
                n1clust = n1clust + 1
                n2clust = n2clust - 1
                z1bar = z1bar/dble(n1clust)
                z2bar = z2bar/dble(n2clust)


                !! Comput between-group-within-group ratio(bwr).

                bwr(s) = (dble(n1clust) * (z1bar - zbar)**2 + dble(n2clust) * (z2bar - zbar)**2) / &
                    (z1ss - dble(n1clust) * z1bar**2 + z2ss - dble(n2clust) * z2bar**2)

             end do

             !! Maximize bwr.

             cutoff = zascend(sl - 1 + maxloc(bwr(sl:su), dim=1))

             !! Start local smoothing within the cluster.

             temp = 0D0

             !! If ftilde1(i, j) <= cutoff.

             if (ftilde1(i, j) <= cutoff) then

               
                do i1 = i - k, i + k
                   do j1 = j - k, j + k

                      if (ftilde1(i1, j1) <= cutoff) then

                         ttemp1 = dble(i1 - i)/dble(k)
                         ttemp2 = dble(j1 - j)/dble(k)
                         ttemp =  ker(ttemp1, ttemp2) * ker1(abs(ftilde1(i1, j1) - cutoff)/abs(zmin - cutoff))
                         fhat(i, j) = fhat(i, j) + z1(i1, j1) * ttemp
                         temp = temp + ttemp

                      end if

                   end do
                end do

                fhat(i, j) = fhat(i, j)/temp

             else

                !! (i, j) is in cluster 2.

                do i1 = i - k, i + k
                   do j1 = j - k, j + k

                      if (ftilde1(i1, j1) > cutoff) then

                         ttemp1 = dble(i1 - i)/dble(k)
                         ttemp2 = dble(j1 - j)/dble(k)
                         ttemp =  ker(ttemp1, ttemp2) * ker1(abs(ftilde1(i1, j1) - cutoff)/abs(zmax - cutoff))
                         fhat(i, j) = fhat(i, j) + z1(i1, j1) * ttemp
                         temp = temp + ttemp

                      end if

                   end do
                end do

                fhat(i, j) = fhat(i, j)/temp

             end if

          else

             fhat(i, j) = ahat

          end if

       end do
    end do

    !! Save the deblurred image

    do i = k, n + k
       do j = k, n + k

          estImg(i-k, j-k) = fhat(i, j)

       end do
    end do


  end Subroutine cluster_cwm_deblur




 
