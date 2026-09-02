      real function mean(x,num)
*
* calculation of the mean value from variable x
*
      implicit none

      integer i,num
      real x(num)
      real sum

      sum=0.
      do 10 i=1,num
         sum=sum+x(i)
10    continue
      mean=sum/float(num)

      return
      end
