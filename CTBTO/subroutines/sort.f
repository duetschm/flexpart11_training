      subroutine sort(x,index,num,opt)
*
*************************************************************
*     DESCRIPTION: This program sorts the real array x      *
*     OPT=1: In numerical order from small to big numbers   *
*     OPT=-1: In numerical order from big to small numbers  *
*     index: integer array, which stores the indices:       *
*                                                           *
*     INDEX(J): Index of X(J) before execution of SORT      *
*                                                           *
*     NUM: Dimension of X and INDEX                         *
*************************************************************
*
      implicit none
      integer num
      real x(num),b
      integer index(num),opt,i,j,l
*
      if((opt.ne.1).and.(opt.ne.-1)) then
      write(*,*) ' ERROR: FAULT IN DETERMINATION OF SORT DIRECTION '
      write(*,*) ' SORT DIRECTION HAS TO BE SET 1 OR -1            '
      write(*,*) ' PLEASE CHECK !!!                                '
      return
      endif

      do 5 i=1,num
         index(i)=i
5     continue
*
      if(opt.eq.1) then
*
*     X(I+1) >= X(I)
*
      do 10 i=1,num-1
         do 10 j=1,num-i

            if(x(j).gt.x(j+1)) then
               b=x(j)
               l=index(j)
               x(j)=x(j+1)
               index(j)=index(j+1)
               x(j+1)=b
               index(j+1)=l
            endif

10    continue
*
      endif
*
      if(opt.eq.-1) then
*
*     X(I+1) <= X(I)
*
      do 20 i=1,num-1
         do 20 j=1,num-i

            if(x(j).lt.x(j+1)) then
               b=x(j)
               l=index(j)
               x(j)=x(j+1)
               index(j)=index(j+1)
               x(j+1)=b
               index(j+1)=l
            endif

20    continue
*
      endif
*
      return
      end
