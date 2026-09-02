      real function corr(x,xmean,y,ymean,num)
*
* calculation of the correlation coefficient between variable x and y
*
* INPUT:
*
* x(num)    variable 1
* xmean     mean of variable 1
* y(num)    variable 2
* ymean     mean of variable 2
* num       number of values
*
      implicit none

      integer num
      real x(num),y(num),xmean,ymean
      real covx,varx

      corr =
     &       covx(x,xmean,y,ymean,num)/
     &       (sqrt(varx(x,xmean,num)*varx(y,ymean,num)))

      return
      end
