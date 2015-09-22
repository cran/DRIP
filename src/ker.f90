!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!								    !
!  Define 2-D Epanechnikov kernel function                          !
!				                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function ker(x, y)
  double precision :: x, y, pi

  pi = 3.141592653589793D0

  ker = 0D0

  if (x**2 + y**2 <= 1D0) then

     ker = (1 - x**2 - y**2) * 2D0/pi

  end if

end function ker
