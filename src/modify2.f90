!  THIS PROGRAM IS TO DELETE SOME JUMP CANDIDATES THAT ARE        
!  SCATTERED IN THE DESIGN SPACE. WE USE THE CITERION THAT IF THE 
!  NUMBER OF JUMP CANDIDATES IN A NEIGHBORHOOD OF A GIVEN JUMP    !
!  CANDIDATE IS NOT BIGGER THAN K/2 THEN WE DELETE THAT JUMP CAN- !
!  DIDATE.                                                        !
!                                                                 


subroutine modify2(n, k, bound, edge)

  integer :: i, j, i1, j1, n, k, k1, sum, bound, edge(0:n + 2 * bound, 0:n + 2 * bound)

  k1 = (k + 1)/2

  !  ********************** modify four corners ************************

  !  *******************  up-left  **********************

  sum = 0

  do i = bound + 1, bound + k
     do j = bound + 1, bound + k

        sum = sum + edge(i, j)

     end do
  end do

  if (sum <= k1) then

     do  i = bound + 1, bound + k1 - 1
        do  j=bound + 1, bound + k1 - 1

           edge(i, j) = 0

        end do
     end do

  end if

  !  *******************  up-right  *********************

  sum = 0

  do  i = n + bound - k + 1, n + bound
     do j = bound + 1, bound + k

        sum = sum + edge(i, j)

     end do
  end do

  if (sum <= k1) then

     do  i = n + bound - k1 + 2, n + bound
        do  j = bound + 1, bound + k1 - 1

           edge(i, j) = 0

        end do
     end do

  end if

  !  ******************  down-left  *********************

  sum = 0

  do i = bound + 1, bound + k
     do j = bound + n- k + 1, bound + n

        sum = sum + edge(i, j)

     end do
  end do

  if (sum <= k1) then

     do i = bound + 1, bound + k1 - 1
        do j = n + bound - k1 + 2, n + bound

           edge(i, j) = 0

        end do
     end do

  end if

  !  ******************  down-right  *********************

  sum = 0

  do i = bound + n - k + 1, bound + n
     do j = bound + n - k + 1, bound + n

        sum = sum + edge(i, j)

     end do
  end do

  if (sum <= k1) then

     do i=n+bound-k1+2, n+bound
        do j=n+bound-k1+2, n+bound

           edge(i,j)=0

        end do
     end do

  end if


  !  ***************** Modify four boundaries *********************

  !  ******************  upper  ********************

  do i=bound+k1, bound+n-k1+1
     do j=bound+1, bound+k1-1

        sum=0

        if (edge(i,j) .EQ. 1) then

           do i1=i-k1+1, i+k1-1
              do j1=bound+1, bound+k

                 sum=sum+edge(i1,j1)

              end do
           end do

           if (sum .LE. k1) then

              edge(i,j)=0

           end if

        end if

     end do
  end do


  !  ******************* down  ********************

  do i=bound+k1, bound+n-k1+1
     do j=n+bound-k1+2, bound+n

        sum=0

        if (edge(i,j) .EQ. 1) then

           do i1=i-k1+1, i+k1-1
              do j1=bound+n-k+1, bound+n

                 sum=sum+edge(i1,j1)

              end do
           end do

           if (sum .LE. k1) then

              edge(i,j)=0

           end if

        end if

     end do
  end do


  !  *******************  left  ********************

  do j=bound+k1, bound+n-k1+1
     do i=bound+1, bound+k1-1

        sum=0

        if (edge(i,j) .EQ. 1) then

           do j1=j-k1+1, j+k1-1
              do i1=bound+1, bound+k

                 sum=sum+edge(i1,j1)

              end do
           end do

           if (sum .LE. k1) then

              edge(i,j)=0

           end if

        end if

     end do
  end do


  !  *******************  right  ********************

  do j=bound+k1, bound+n-k1+1
     do i=bound+n-k1+2,bound+n

        sum=0

        if (edge(i,j) .EQ. 1) then

           do j1=j-k1+1, j+k1-1
              do i1=bound+n-k+1, bound+n

                 sum=sum+edge(i1,j1)

              end do
           end do

           if (sum .LE. k1) then

              edge(i,j)=0

           end if

        end if

     end do
  end do


  !  *******************  center  ********************

  do i=bound+k1,bound+n-k1+1
     do j=bound+k1,bound+n-k1+1

        sum=0

        if (edge(i,j) .EQ. 1) then

           do i1=i-(k-1)/2,i+(k-1)/2
              do j1=j-(k-1)/2,j+(k-1)/2

                 sum=sum+edge(i1,j1)

              end do
           end do

           if (sum .LE. (k-1)/2) then

              edge(i,j)=0

           end if

        end if

     end do
  end do


end subroutine modify2
