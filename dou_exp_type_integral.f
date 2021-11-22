      program double
      implicit real*8 (a-h,o-z)
      parameter(nn=10000)
      dimension fac(0:nn),fp(0:nn),xp(0:nn),xq(0:nn)
      data ne,rl, a, b/100, 10.d0, -1.d0, 1.d0/
      pi=dacos(-1.d0)
      h=2.d0*rl/dble(ne)
      do 70 j=0, ne
      fac(j)=1.d0
      if(j.eq.0.or.j.eq.ne)then
        fac(j)=0.5d0
      endif
      xp(j)=-rl+h*dble(j)
      xpp=xp(j)
      yy=dtanh(pi/2.d0*dsin(xpp))
      xq(j)=(b-a)/2.d0*yy+(a+b)/2.d0
      if(xpp.ge.0.d0)then
        fp(j)=2.d0*(1.d0+dexp(-2.d0*xpp))*dexp(-pi*dsinh(xpp))
     &         /(dexp(-xpp)+2.d0*dexp(-xpp-pi*dsinh(xpp))
     &         +dexp(-xpp-2.d0*pi*dsinh(xpp)))
      else
        fp(j)=2.d0*(1.d0+dexp(2.d0*xpp))*dexp(pi*dsinh(xpp))
     &         /(dexp(xpp)+2.d0*dexp(xpp+pi*dsinh(xpp))
     &         +dexp(xpp+2.d0*pi*dsinh(xpp)))
      endif
  70  continue
      sss=0.d0
      do 10 j=0, ne
        xqq=xq(j)
        call func(xqq,ff)
        sss=sss+fac(j)*ff*fp(j)
  10  continue
      sss=sss*h*pi/4.d0*(b-a)
      write(6,100) ne,sss
 100  format(2x,' ne = ', i5, ' sss = ', 1pe16.9)
      end program double

      subroutine func(x,f)
      implicit real*8 (a-h,o-z)
      pi=dacos(-1.d0)
      f=dsin(pi/2.d0*(x+1.d0))
      return
      end subroutine func

