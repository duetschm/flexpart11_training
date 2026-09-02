      real function varx(x,xmean,num)
*
* calculation of the empirical variance of variable x
*
      implicit none

      integer i,num
      real x(num),xmean
      real sum

      sum=0.

      do 10 i=1,num

         sum=sum+(x(i)-xmean)**2

10    continue

      varx=sum/float(num)

      return
      end
