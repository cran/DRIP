!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define 1-D kernel function                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DOUBLE PRECISION function ker1(x)
      DOUBLE PRECISION :: x, pi

      pi = 3.141592653589793D0

      ker1 = 0D0

      if (x >= 0D0 .AND. x <= 1D0) then

          ker1 = exp(x**2/2D0)/1.194958D0
!         ker1 = x**2/3D0

      end if

      end function ker1
