  ! This is a Fortran 90 subroutine for local median filter (LMF) in image processing.

subroutine LocalMedianFilter(n, k, cw, z, ftilde)

  implicit none

  integer :: n, i, j, k, i1, j1, cw, nitem, l

  double precision :: z(0:n, 0:n), z1(0:(n+2*k), 0:(n+2*k)), ftilde(0:n, 0:n), zascend(1:((n+1)**2)), &
       zmed

  external :: extend, median

  
  !! Extend to avoid boundary problems.

  call extend(n, k, z, z1)

  !! Start local filtering.

  do i = k, n + k
     do j = k, n + k

        ftilde(i-k, j-k) = 0D0

        nitem = 0

        do i1 = i - k, i + k
           do j1 = j - k, j + k

              if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then
                     
                 nitem = nitem + 1
                 zascend(nitem) = z1(i1, j1)

              end if

           end do
        end do

        !! Replicate the center observation.

        if (cw >= 2) then
           
           do l = 1, cw - 1

              nitem = nitem + 1
              zascend(nitem) = z1(i, j)

           end do

        end if

        !! Order the observations in the neighborhood.

        ! do l1 = 1, nitem - 1

        !    small = zascend(l1)
        !    locsm = l1

        !    do l2 = l1 + 1, nitem

        !       if (zascend(l2) < small) then

        !          small = zascend(l2)
        !          locsm = l2

        !       end if

        !    end do

        !    zascend(locsm) = zascend(l1)
        !    zascend(l1) = small

        ! end do
         
        ! if (MOD(nitem, 2) == 0) then

        !    ftilde(i-k, j-k) = (zascend(nitem/2) + zascend(nitem/2 + 1)) / 2D0

        ! else

        !    ftilde(i-k, j-k) = zascend((nitem+1)/2)

        ! end if

        call median(zascend(1:nitem), nitem, zmed)

        ftilde(i-k, j-k) = zmed

                
     end do
  end do

  
end subroutine LocalMedianFilter
