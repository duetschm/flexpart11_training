      subroutine regression(x,xmean,y,ymean,num,a,b,r)
*
* calculation of the parameters a and b for linear correlation of
* the variables x and y
*
* Y = A + B*X
* r...Korrelationskoeffizient
*
      implicit none

      integer num
      real x(num),y(num),xmean,ymean,a,b,r
      real covx,varx,corr

      b=covx(x,xmean,y,ymean,num)/varx(x,xmean,num)
      a=ymean-b*xmean
      r=corr(x,xmean,y,ymean,num)

      return
      end
