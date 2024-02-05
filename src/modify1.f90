! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     C                                                                CC
!     C  THIS SUBROUTINE MODIFIES THE DETECTED JUMP CANDIDATE SET SUCH CC
!     C  THAT SOME DECEPTIVE JUMP CANDIDATES ARE DELETED AND THE       CC
!     C  DETECTED JUMP LOCATIONS ARE MADE THIN.                        CC
!     C                                                                CC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CCCCCCCCCCCCDetect the edge candidate points CCCCCCCCCCCCCCCCCC


    SUBROUTINE modify1(n,k,bound,z,edge)

    integer :: n, k, k1, i, j, bound, loc1, loc2
    integer :: edge(0:n + 2 * bound, 0:n + 2 * bound)
    double precision :: z(0:n + 2 * bound, 0:n + 2 * bound)
    double precision :: beta0(0:n+2*bound,0:n+2*bound)
    double precision :: beta1(0:n+2*bound,0:n+2*bound)
    double precision :: beta2(0:n+2*bound,0:n+2*bound)
    double precision :: zdot(0:n+2*bound),ssquare
    double precision :: zidot, zjdot, temp1, temp2, temp3
    double precision :: temp4, temp5, temp6, cond


    k1=(k+1)/2

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     C
!     To calculate S-square.             C
!     C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    ssquare=0D0

    do 10 i=1,k
        ssquare=ssquare+dble(i-k1)**2D0
    10 END DO

    ssquare=ssquare/(dble(n)**2D0)


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     C
!     construct the least square estimates of the slopes. Use updateC
!     schedule. To calculate the slopes at (k1,k1) first. Then      C
!     update to calculate slopes at (k1,j) for j=k1+1,n-k1+1. Then  C
!     to calculate slopes at (i,j), i=k1+1,n-k1+1, j=k1,n-k1+1.     C
!     C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! CCCCCCCCCCTo calculate slopes at position (k1,k1)  CCCCCCCCCCC


    beta1(k1,k1)=0D0
    beta2(k1,k1)=0D0
    beta0(k1,k1)=0D0

!     **************************************************************
!     *****  zdot(j) is z.., summation of the observations in  *****
!     *****  the neighborhood. beta1(.,.) and beta2(.,.) are   *****
!     *****  slopes wrt x and y directions. We use the symmetry*****
!     *****  between beta1 and beta2 to calculate them.        *****
!     **************************************************************

    zdot(k1)=0D0

    do 20 i=1,k

        zidot=0D0
        zjdot=0D0

        do 30 j=1,k
            zidot=zidot+z(i,j)
            zjdot=zjdot+z(j,i)
            zdot(k1)=zdot(k1)+z(i,j)
        30 END DO

        beta1(k1,k1)=beta1(k1,k1)+dble(i-k1)/dble(n)*zidot
        beta2(k1,k1)=beta2(k1,k1)+dble(i-k1)/dble(n)*zjdot

    20 END DO

    beta0(k1,k1)=zdot(k1)/(dble(k)*dble(k))
    beta1(k1,k1)=beta1(k1,k1)/(dble(k)*ssquare)
    beta2(k1,k1)=beta2(k1,k1)/(dble(k)*ssquare)


! CCCCCCTo calculate slopes at positions (k1,j), any j. CCCCCCCC


    do 40 j=k1+1,n+2*bound-k1+1

        temp1=0D0
        temp2=0D0

        do 50 i=1,k
            temp1=temp1+dble(i)*(z(i,j+(k1-1))-z(i,j-k1))
            temp2=temp2+(z(i,j+(k1-1))-z(i,j-k1))
        50 END DO

        beta1(k1,j)=beta1(k1,j-1)+(temp1-dble(k1)*temp2)/ &
        (dble(n)*dble(k)*ssquare)

        zdot(j)=zdot(j-1)+temp2
        beta0(k1,j)=zdot(j)/(dble(k)*dble(k))

        temp3=0D0

        do 60 i=1,k
            temp3=temp3+dble(j+(k1-1))*z(i,j+(k1-1)) &
            -dble(j-k1)*z(i,j-k1)
        60 END DO

        beta2(k1,j)=beta2(k1,j-1)+(temp3-zdot(j)-dble(j-1)*temp2)/ &
        (dble(n)*dble(k)*ssquare)

    40 END DO


! CCCCTo calculate the slopes at position (i,j), any i,j. CCCCCC

          
    do 70 j=k1,n+2*bound-k1+1

        do 80 i=k1+1,n+2*bound-k1+1

            temp4=0D0
            temp5=0D0

            do 90 j1=j-k1+1, j+k1-1
                temp4=temp4+j1*(z(i+k1-1,j1)-z(i-k1,j1))
                temp5=temp5+(z(i+k1-1,j1)-z(i-k1,j1))
            90 END DO

            beta2(i,j)=beta2(i-1,j)+(temp4-dble(j)*temp5)/ &
            (dble(n)*dble(k)*ssquare)

            zdot(j)=zdot(j)+temp5
            beta0(i,j)=zdot(j)/(dble(k)*dble(k))

            temp6=0D0

            do 100 j1=j-k1+1, j+k1-1
                temp6=temp6+dble(i+k1-1)*z(i+k1-1,j1)- &
                dble(i-k1)*z(i-k1,j1)

            100 END DO

            beta1(i,j)=beta1(i-1,j)+(temp6-zdot(j)-dble(i-1)*temp5) &
            /(dble(n)*dble(k)*ssquare)

        80 END DO
    70 END DO


!     *******************************************************************
!     *** To make the edges thinner, we do the following. (1) Along x ***
!     *** direction, we modify the candidates with gradient angle     ***
!     *** -45<=theta<45 or 135<=theta<225. The modification procedure ***
!     *** is from Qiu (1994). (2) Along y direction, we deal with the ***
!     *** candidates with 45<=theta<135 or 225<=theta<315.            ***
!     *******************************************************************

    n=n+2*bound
    loc2=0

    do 210 i=bound+1, n-bound

        loc1=n-bound+2

        do 220 j=bound+1, n-bound
            cond=dabs(beta2(i,j))/dsqrt(beta1(i,j)**2D0+ &
            beta2(i,j)**2D0)
            if (edge(i,j) == 1 .AND. cond > dsqrt(2D0)/2D0) then
                loc1=j
                loc2=j
                goto 230
            endif
        220 END DO

        230 do 240 j=loc1+1, n-bound

            cond=dabs(beta2(i,j))/dsqrt(beta1(i,j)**2D0+ &
            beta2(i,j)**2D0)
                        
            if (edge(i,j) == 1 .AND. j > loc2+2*bound-1 .AND. &
            cond > dsqrt(2D0)/2D0) then
                edge(i,loc1)=0
                edge(i,(loc1+loc2)/2)=1
                loc1=j
                loc2=j
                goto 230
            endif
            	    
            if (edge(i,j) == 1 .AND. j <= loc2+2*bound-1 .AND. &
            cond > dsqrt(2D0)/2D0) then
                loc2=j
                edge(i,j)=0
                goto 240
            endif

            if (j == n-bound) then
                edge(i,loc1)=0
                edge(i,(loc1+loc2)/2)=1
            endif

        240 END DO
    210 END DO

    do 170 j=bound+1, n-bound

        loc1=n-bound+2

        do 180 i=bound+1, n-bound
            cond=dabs(beta2(i,j))/dsqrt(beta1(i,j)**2D0+ &
            beta2(i,j)**2D0)
            if (edge(i,j) == 1 .AND. cond <= dsqrt(2D0)/2D0) then
                loc1=i
                loc2=i
                goto 190
            endif

        180 END DO

        190 do 200 i=loc1+1, n-bound

            cond=dabs(beta2(i,j))/dsqrt(beta1(i,j)**2D0+ &
            beta2(i,j)**2D0)

            if (edge(i,j) == 1 .AND. i > loc2+2*bound-1 .AND. &
            cond <= dsqrt(2D0)/2D0) then
                edge(loc1,j)=0
                edge((loc1+loc2)/2,j)=1
                loc1=i
                loc2=i
                goto 190
            endif
            	    
            if (edge(i,j) == 1 .AND. i <= loc2+2*bound-1 .AND. &
            cond <= dsqrt(2D0)/2D0) then
                loc2=i
                edge(i,j)=0
                goto 200
            endif

            if (i == n-bound) then
                edge(loc1,j)=0
                edge((loc1+loc2)/2,j)=1
            endif

        200 END DO
    170 END DO

    n=n-2*bound

    end SUBROUTINE modify1
