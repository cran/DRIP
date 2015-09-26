!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
! This is Fortran subroutine to compute the edge detector   !
! performance measure d_KQ proposed in the paper            !
! Kang, Y. and Qiu, P., 'Jump Detection in Blurred          !
!                        Regression Surfaces'.              !
! Creator: Yicheng Kang                                     !
! Date: April 24, 2013                                      !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine d_KQ(n, edge1, edge2, dKQ)

  implicit none

  integer :: n, edge1(0:n, 0:n), edge2(0:n, 0:n), nedge1, &
       nedge2, i, j, i1, j1

  double precision :: x, y, x1, y1, dKQ, dist1, dist2, d1, &
       d2

  ! Calculate performance measure                       

  dist1 = 0D0
  nedge1 = 0
  dist2 = 0D0
  nedge2 = 0

  do i = 0, n
     do j = 0, n

        if (edge1(i,j) == 1) then

           nedge1 = nedge1 + 1
           x = dble(i)/dble(n)
           y = dble(j)/dble(n)
           d1 = 2D0

           do i1 = 0, n
              do j1 = 0, n

                 if (edge2(i1, j1) == 1) then

                    x1 = dble(i1)/dble(n)
                    y1 = dble(j1)/dble(n)
                    d1 = min(d1, sqrt((x - x1)**2 + (y - y1)**2))

                 end if

              end do
           end do

           dist1 = dist1 + d1

        end if

     end do
  end do

  dist1 = dist1/dble(nedge1)

  do i = 0, n
     do j = 0, n

        if (edge2(i,j) == 1) then

           nedge2 = nedge2 + 1
           x = dble(i)/dble(n)
           y = dble(j)/dble(n)
           d2 = 2D0

           do i1 = 0, n
              do j1 = 0, n

                 if (edge1(i1, j1) == 1) then

                    x1 = dble(i1)/dble(n)
                    y1 = dble(j1)/dble(n)
                    d2 = min(d2, sqrt((x - x1)**2 + (y - y1)**2))

                 end if

              end do
           end do

           dist2 = dist2 + d2

        end if

     end do
  end do

  dist2 = dist2/dble(nedge2)

  dKQ = dist1 + dist2

end subroutine d_KQ
