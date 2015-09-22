CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                 C
C  THIS PROGRAM IS TO DELETE SOME JUMP CANDIDATES THAT ARE        C
C  SCATTERED IN THE DESIGN SPACE. WE USE THE CITERION THAT IF THE C
C  NUMBER OF JUMP CANDIDATES IN A NEIGHBORHOOD OF A GIVEN JUMP    C
C  CANDIDATE IS NOT BIGGER THAN K/2 THEN WE DELETE THAT JUMP CAN- C
C  DIDATE.                                                        C
C                                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE modify2(n,k,bound,edge)
      
      INTEGER i,j,i1,j1,n,k,k1,sum,bound,edge(0:300,0:300)

      k1=(k+1)/2

C  ********************** modify four corners ************************

C  *******************  up-left  **********************

      sum=0

      do 270 i=bound+1, bound+k
         do 280 j=bound+1, bound+k

            sum=sum+edge(i,j)

 280        continue
 270     continue

      if (sum .LE. k1) then

         do 300 i=bound+1, bound+k1-1
            do 310 j=bound+1, bound+k1-1

               edge(i,j)=0

 310           continue
 300        continue

      endif

C  *******************  up-right  *********************

      sum=0

      do 320 i=n+bound-k+1, n+bound
         do 330 j=bound+1, bound+k

            sum=sum+edge(i,j)

 330        continue
 320     continue

      if (sum .LE. k1) then

         do 340 i=n+bound-k1+2, n+bound
            do 350 j=bound+1, bound+k1-1

               edge(i,j)=0

 350           continue
 340        continue

      endif

C  ******************  down-left  *********************

      sum=0

      do 360 i=bound+1, bound+k
         do 370 j=bound+n-k+1, bound+n

            sum=sum+edge(i,j)

 370        continue
 360     continue

      if (sum .LE. k1) then

         do 380 i=bound+1, bound+k1-1
            do 390 j=n+bound-k1+2, n+bound

               edge(i,j)=0

 390           continue
 380        continue

      endif

C  ******************  down-right  *********************

      sum=0

      do 400 i=bound+n-k+1, bound+n
         do 410 j=bound+n-k+1, bound+n

            sum=sum+edge(i,j)

 410        continue
 400     continue

      if (sum .LE. k1) then

         do 420 i=n+bound-k1+2, n+bound
            do 430 j=n+bound-k1+2, n+bound

               edge(i,j)=0

 430           continue
 420        continue

      endif


C  ***************** Modify four boundaries *********************

C  ******************  upper  ********************

      do 440 i=bound+k1, bound+n-k1+1
         do 450 j=bound+1, bound+k1-1

            sum=0

            if (edge(i,j) .EQ. 1) then

               do 460 i1=i-k1+1, i+k1-1
                  do 470 j1=bound+1, bound+k
                     sum=sum+edge(i1,j1)
 470                 continue
 460              continue

               if (sum .LE. k1) then
                  edge(i,j)=0
                  endif

               endif

 450        continue
 440     continue


C  ******************* down  ********************

      do 480 i=bound+k1, bound+n-k1+1
         do 490 j=n+bound-k1+2, bound+n

            sum=0

            if (edge(i,j) .EQ. 1) then

               do 500 i1=i-k1+1, i+k1-1
                  do 510 j1=bound+n-k+1, bound+n
                     sum=sum+edge(i1,j1)
 510                 continue
 500              continue

               if (sum .LE. k1) then
                  edge(i,j)=0
                  endif

               endif

 490        continue
 480     continue


C  *******************  left  ********************

      do 520 j=bound+k1, bound+n-k1+1
         do 530 i=bound+1, bound+k1-1

            sum=0

            if (edge(i,j) .EQ. 1) then

               do 540 j1=j-k1+1, j+k1-1
                  do 550 i1=bound+1, bound+k
                     sum=sum+edge(i1,j1)
 550                 continue
 540              continue

               if (sum .LE. k1) then
                  edge(i,j)=0
                  endif

               endif

 530        continue
 520     continue


C  *******************  right  ********************

      do 560 j=bound+k1, bound+n-k1+1
         do 570 i=bound+n-k1+2,bound+n

            sum=0

            if (edge(i,j) .EQ. 1) then

               do 580 j1=j-k1+1, j+k1-1
                  do 590 i1=bound+n-k+1, bound+n
                     sum=sum+edge(i1,j1)
 590                 continue
 580              continue

               if (sum .LE. k1) then
                  edge(i,j)=0
                  endif

               endif

 570        continue
 560     continue


C  *******************  center  ********************

      do 10 i=bound+k1,bound+n-k1+1
         do 20 j=bound+k1,bound+n-k1+1

            sum=0

            if (edge(i,j) .EQ. 1) then

               do 30 i1=i-(k-1)/2,i+(k-1)/2
                  do 40 j1=j-(k-1)/2,j+(k-1)/2
                     sum=sum+edge(i1,j1)
 40                 continue
 30              continue

               if (sum .LE. (k-1)/2) then
                  edge(i,j)=0
                  endif

               endif

 20        continue
 10     continue     


        end
