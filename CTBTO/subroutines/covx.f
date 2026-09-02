      real function covx(x,xmean,y,ymean,num)
*
* calculation of the empirical covariance of variable x and y
*
      implicit none

      integer i,num
      real x(num),y(num),xmean,ymean
      real sum

      sum=0.
      do 10 i=1,num
         sum=sum+(x(i)-xmean)*(y(i)-ymean)
10    continue
      covx=sum/float(num)

      return
      end
