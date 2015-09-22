!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
!  Define 1-D kernel function                                       !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function ker1(x)
  double precision :: x

  ker1 = 0D0

  if (x >= 0D0 .AND. x <= 1D0) then

     ! Truncated Gaussian kernel

     !          ker1 = exp(x**2/2D0)/1.194958D0

     ! Epanechnikov kernel

     ker1 = 3D0 * x**2

  end if

end function ker1
