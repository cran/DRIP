!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
!  This subroutine extends the data matrix so that there is no       !
!  boundary problem when we detect edges, etc., in regular region.   !
!                                                                    !
!     Creator: Peihua Qiu                                            !
!     Date:    November 24, 2008                                     !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine extend(n,k,z,z1)

  implicit none

  INTEGER :: i, j, n,k
  DOUBLE PRECISION :: z(0:n,0:n),z1(0:(n+2*k),0:(n+2*k))

  do i=k,k+n
     do j=k,k+n
        z1(i,j)=z(i-k,j-k)
     end do
  end do

  do i=0,k-1
     do j=0,k-1
        z1(i,j)=z(k-1-j,k-1-i)
     end do
  end do

  do i=k,k+n
     do j=0,k-1
        z1(i,j)=z(i-k,k-1-j)
     end do
  end do

  do i=k+n+1,2*k+n
     do j=0,k-1
        z1(i,j)=z(j+n-k+1,i-n-k-1)
     end do
  end do

  do i=k+n+1,2*k+n
     do j=k,k+n
        z1(i,j)=z(k+2*n+1-i,j-k)
     end do
  end do

  do i=k+n+1,2*k+n
     do j=k+n+1,2*k+n
        z1(i,j)=z(k+2*n+1-j,k+2*n+1-i)
     end do
  end do

  do i=k,k+n
     do j=k+n+1,2*k+n
        z1(i,j)=z(i-k,k+2*n+1-j)
     end do
  end do

  do i=0,k-1
     do j=k+n+1,2*k+n
        z1(i,j)=z(j-k-n-1,n-k+1-i)
     end do
  end do

  do i=0,k-1
     do j=k,k+n
        z1(i,j)=z(k-1-i,j-k)
     end do
  end do

end Subroutine extend
