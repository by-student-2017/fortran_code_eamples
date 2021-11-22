      program epsilon
      real*8 a,b,c,eps
         a = 4.0d0/3.0d0
         b = a - 1.0d0
         c = b + b + b
         eps = dabs(c - 1.d0)
         write(*,*) eps
         stop
      end program epsilon
