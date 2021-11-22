      program legtst
      implicit real*8 (a-h,o-z)
      dimension xp(0:200),xq(0:200)
        m=40
        call leg(m,xp,xq)
        sss=0.d0
        do 10 k=1, m-1
          sss=sss+xq(k)*dsqrt(1.d0-xp(k)**2)
   10   continue
        write(6,*) '2*sss = ',2*sss
      stop
      end program legtst

      subroutine leg(m,xp,xq)
      implicit real*8 (a-h,o-z)
      dimension xp(0:m),xq(0:m),p(3),q(3),r(3)
        n=m-1
        pi=4*datan(1.d0)
        eps=1.d-14
        do 50 k=1,n
          kp1=k+1
          km1=k-1
          if(n/2*2.eq.n)then
            mm=n/2
            if(k.lt.mm)then
              xn=dble(n+1-2*k)/dble(2*n+1)
              xnp1=dble(n+1-2*kp1)/dble(2*n+1)
            else if(k.eq.mm)then
              xn=dble(n+1-2*k)/dble(2*n+1)
              xnp1=0.d0
            else if(k.eq.mm+1)then
              xn=0.d0
              xnp1=dble(n+1-2*k)/dble(2*n+1)
            else
              xn=dble(n+1-2*km1)/dble(2*n+1)
              xnp1=dble(n+1-2*k)/dble(2*n+1)
            endif
          else
            mm=(n+1)/2
            if(k.lt.mm-1)then
              xn=dble(n+1-2*k)/dble(2*n+1)
              xnp1=dble(n+1-2*kp1)/dble(2*n+1)
            else if(k.eq.mm-1)then
              xn=dble(n+1-2*k)/dble(2*n+1)
              xnp1=-eps
            else if(k.eq.mm)then
              xn=0.d0
              xnp1=0.d0
            else if(k.eq.mm+1)then
              xn=eps
              xnp1=dble(n+1-2*k)/dble(2*n+1)
            else
              xn=dble(n+1-2*km1)/dble(2*n+1)
              xnp1=dble(n+1-2*k)/dble(2*n+1)
            endif
          endif
          x=dsin(xn*pi)
          xp1=dsin(xnp1*pi)
          write(6,*) ' k = ', k, ' x = ', x, ' xp1 = ', xp1
   15     xmd=(x+xp1)/2.d0
          if(dabs(x-xp1).lt.eps) go to 30
   10       p(1)=1.d0
            p(2)=x
            q(1)=1.d0
            q(2)=xp1
            r(1)=1.d0
            r(2)=xmd
            do 20 i=2,n
              p(3)=x*p(2)+dble(i-1)/dble(i)*(x*p(2)-p(1))
              p(1)=p(2)
              p(2)=p(3)
              q(3)=xp1*q(2)+dble(i-1)/dble(i)*(xp1*q(2)-q(1))
              q(1)=q(2)
              q(2)=q(3)
              r(3)=xmd*r(2)+dble(i-1)/dble(i)*(xmd*r(2)-r(1))
              r(1)=r(2)
              r(2)=r(3)
   20       continue
            fac1=p(2)*r(2)
            fac2=q(2)*r(2)
            if(fac1.lt.0.d0)then
              xp1=xmd
              go to 15
            endif
            if(fac2.lt.0.d0)then
              x=xmd
              go to 15
            endif
            if(fac2.lt.0.d0)then
              x=xmd
              go to 15
            endif
   30     continue
   40     w=2.d0*(1.d0-x**2)/(dble(n)*p(1))**2
          y=-x
          xp(k)=y
          xq(k)=w
   50   continue
        write(6,650)(k,xp(k),xq(k),k=1,m-1)
  650   format(2x,'k=',i3,6x,'y=',d14.7,6x,'w=',d14.7)
      return
      end subroutine leg
