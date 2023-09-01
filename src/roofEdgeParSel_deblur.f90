!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                                !
  ! This is a Fortran 90 source code to select roof/valley edge detection          !
  ! parameters (bandwidth and threshold). The selection procedure is proposed in   !
  ! the paper Qiu and Kang "Blind Image Deblurring Using Jump Regression Analysis",!
  ! Statistica Sinica, Volume 25,  Number 3, July 2015                             !
  !                                                                                !
  ! Creator: Yicheng Kang                                                          !
  ! Date: Sep 13 2015                                                              !
  !                                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
subroutine roofEdgeParSel_deblur(n, obsImg, nband, bandwidth, nthresh, thresh, nboot, &
     edge1, dKQ)

  implicit none

  integer :: n, nband, bandwidth(1:nband), nthresh

  double precision :: thresh(1:nthresh)

  integer :: nboot, bandw, k, iband, u, v, &
       ithresh, edge_orig(0:n, 0:n), edge_boot(0:n, 0:n), iboot, edge1(0:n, 0:n), &
       i, j, edge1_ext(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), &
       edge1_ext1(0:(n+2*maxval(bandwidth)), 0:(n+2*maxval(bandwidth))), n_edge1(0:n, 0:n), &
       i1, j1
  
  double precision :: obsImg(0:n, 0:n), llkbw(1:7), sigma, optb1, &
       fbhat(0:n, 0:n), resid(0:n, 0:n), dKQ(1:nband, 1:nthresh), diff_orig(0:n, 0:n), &
       u_temp, v_temp, bootImg(0:n, 0:n), diff_boot(0:n, 0:n), h, dist

  external :: surfest, roofDiff_deblur, d_KQ, extend1

  external :: rndstart, rndend
  double precision :: myrunif
  external :: myrunif


  !! Read in data. Initialize fbhat, residuals and dKQ distances.

  do i = 0, n
     do j = 0, n

        edge1_ext(i, j) = edge1(i, j)
        fbhat(i, j) = 0D0
        resid(i, j) = 0D0

     end do
  end do

  dKQ(1:nband, 1:nthresh) = 0D0
  
  !! Use a preliminary local linear kernel smoothing to obtain residuals.

  sigma = 0D0
  optb1 = 0D0
  bandw = 7 ! The length of llkbw

  do i = 1, 7

     llkbw(i) = dble(i)/dble(n) !search bdwdth for conventnl LLK.

  end do

  call surfest(n, obsImg, bandw, llkbw, sigma, fbhat, resid, optb1)

  !! Iterate through each bandwidth.

  do iband = 1, nband

     k = bandwidth(iband)

     !! Flag the neighborhood if there are step edges.
     
     call extend1(n, k, edge1_ext(0:n, 0:n), edge1_ext1(0:(n+2*k), 0:(n+2*k)))

     do i = k, n + k
        do j = k, n + k

           n_edge1(i-k, j-k) = 0
           
           do i1 = i - k, i + k
              do j1 = j - k, j + k

                 if ((i1 - i)**2 + (j1 - j)**2 <= k**2) then

                    n_edge1(i-k, j-k) = n_edge1(i-k, j-k) + edge1_ext1(i1, j1)

                 end if

              end do
           end do

        end do
     end do
              
     
     do i = 0, n
        do j = 0, n

           diff_orig(i, j) = 0D0

        end do
     end do

     !! Roof/Valley edge detection on the original sample
     
     call roofDiff_deblur(n, obsImg, k, diff_orig)

     !! Roof/Valley edge detection on the bootstrap sample

     call rndstart()

     do iboot = 1, nboot

        do i = 0, n
           do j = 0, n

              !call random_number(u_temp)
              !call random_number(v_temp)
              u_temp = myrunif(0D0, 1D0)
              v_temp = myrunif(0D0, 1D0)
              u = int(u_temp * dble(n+1))
              v = int(v_temp * dble(n+1))
              bootImg(i, j) = fbhat(i, j) + resid(u, v)
              diff_boot(i, j) = 0D0

           end do
        end do

        call roofDiff_deblur(n, bootImg, k, diff_boot)
     
        !! Iterate through each threshold.

        do ithresh = 1, nthresh

           h = thresh(ithresh)

           do i = 0, n
              do j = 0, n

                 if (diff_orig(i, j) > h .and. n_edge1(i, j) == 0) then

                    edge_orig(i, j) = 1

                 else

                    edge_orig(i, j) = 0

                 end if

                 if (diff_boot(i, j) > h .and. n_edge1(i, j) == 0) then

                    edge_boot(i, j) = 1

                 else

                    edge_boot(i, j) = 0

                 end if

              end do
           end do
        
           !! Calculate d_KQ distance.

           dist = 0D0
        
           call d_KQ(n, edge_orig, edge_boot, dist)

           dKQ(iband, ithresh) = dKQ(iband, ithresh) + dist

        end do

     end do

     call rndend()

  end do

  dKQ(1:nband, 1:nthresh) = dKQ(1:nband, 1:nthresh)/dble(nboot)

end subroutine roofEdgeParSel_deblur
