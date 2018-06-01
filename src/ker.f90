!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								    !
!  Define 2-D kernel function                                       !
!				                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DOUBLE PRECISION function ker(x, y)
      DOUBLE PRECISION :: x, y, pi

      pi = 3.141592653589793D0

      ker = 0D0

      if (x**2 + y**2 <= 1D0) then
 !        ker = (exp(-(x**2 + y**2)/2D0) - exp(-0.5D0))/ &
 !            (2D0 * pi - 3D0 * pi * exp(-0.5D0))
         ker = (1 - x**2 - y**2)/(pi/2D0)
      end if

      end function ker
