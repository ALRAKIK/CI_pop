module Tools

      implicit none

      contains

      function erfcx(t) result(sum)

        implicit none 
        double precision, intent(in) :: t
        double precision             :: sum 

        integer                      :: i 
        double precision,parameter   :: pi = 3.14159265358979323846D00

        sum = 0 

        if (t > 10.0d0) then 
         do i = 1 , 20
          sum = sum + (-1.d0)**(i-1) * factorial2((2*i-3)) / ( sqrt(pi) * 2.d0**(i-1) * t**(2*i-1) )
         end do 
        else
          sum = erfc(t) * dexp(t*t)
        end if 

      end function erfcx 


      pure function factorial(x)  result(fac)

      implicit none 

      integer ,intent(in) :: x
      integer             :: fac 

      integer             :: i 

      if (x <=1) then 
        fac = 1
        return 
      end if 

      fac = 1
      do i = x , 2 , -1 
        fac = fac * i 
      end do 

      end function factorial

      pure function factorial2(x)  result(fac)

      implicit none 

      integer ,intent(in) :: x
      integer             :: fac 

      integer             :: i 

      if (x <= 1 ) then 
        fac = 1 
        return 
      end if 

      fac = 1.d0 
      do i = x , 2 , -2
        fac = fac * i 
      end do 

      end function factorial2












end module 