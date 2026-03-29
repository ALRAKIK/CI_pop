module bessel_derivatives
      
      implicit none

      contains

      double complex function der_I_A(s, c, n, I_A,A)

      implicit none
      integer,          intent(in) :: s, c, n
      double precision, intent(in) :: I_A(0:)
      double precision, intent(in) :: A
      double complex, parameter    :: I_dp = (0.0D0, 1.0D0)
      integer                      :: sc

      sc = s*10 + c

      select case (sc)
      
      case(10)
        
        der_I_A = ( n * I_A(n) / A ) * I_dp

        if (A == 0.d0) then
          if  (n == 1) then 
            der_I_A = 0.5d0 * I_dp
          else 
            der_I_A = 0.d0
          end if 
        end if 

      case(01)
        der_I_A =           0.50d0 * ( I_A(abs(n-1)) + I_A(n+1) )

      case(11)
        
        der_I_A =           0.50d0 * ( (n+1)*I_A(n+1) + (n-1)*I_A(abs(n-1)) ) / A * I_dp

        if (A == 0.d0) then
          if  (n == 2) then 
            der_I_A = 0.25d0 * I_dp
          else 
            der_I_A = 0.d0
          end if 
        end if 

      case(20)
         !der_I_A = I_A(n) - 0.25d0 * ( I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(n+2) )  ! care sin^2(x) = 1 - cos^2(x)
         
         der_I_A = 0.50d0 * ( (n+1)*I_A(n+1) - (n-1)*I_A(abs(n-1)) ) / A

        if (A == 0.d0) then
          if  (n == 2) then 
            der_I_A = 0.25d0
          else 
            der_I_A = 0.d0
          end if 
        end if 

      case(02)
        der_I_A =           0.25d0 * ( I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(n+2) )
      
      case default
        der_I_A = 0.d0
      end select

      end function der_I_A


      double complex function der_I_B(s, c, n, I_B,B)

        implicit none

      integer,          intent(in) :: s, c, n
      double precision, intent(in) :: I_B(0:)
      double precision, intent(in) :: B
      double complex, parameter    :: I_dp = (0.0D0, 1.0D0)
      integer                      :: sc
      
      sc = s*10 + c

      select case (sc)
        
        case(10)
          
          der_I_B = - ( n * I_B(n) / B ) * I_dp

           if (B == 0.d0) then 
             if (n == 1) then 
               der_I_B = - 0.5d0 * I_dp
             else
               der_I_B = 0.d0
             end if 
           end if 

        case(01)
          der_I_B =          0.50d0 * ( I_B(abs(n-1)) + I_B(n+1) )
        case(11)
          der_I_B =        - 0.50d0 * ( (n+1)*I_B(n+1) + (n-1)*I_B(abs(n-1)) ) / B * I_dp

        if (B == 0.d0) then
          if  (n == 2) then 
            der_I_B = 0.25d0 * I_dp
          else 
            der_I_B = 0.d0
          end if 
        end if 

        case(20)
          !der_I_B = I_B(n) - 0.25d0 * ( I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(n+2) ) ! care sin^2(x) = 1 - cos^2(x)
          der_I_B = 0.50d0 * ( (n+1)*I_B(n+1) - (n-1)*I_B(abs(n-1)) ) / B

        if (B == 0.d0) then
          if  (n == 2) then 
            der_I_B = 0.25d0
          else 
            der_I_B = 0.d0
          end if 
        end if



        case(02)
          der_I_B =          0.25d0 * ( I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(n+2) )
        
        case default
          der_I_B = 0.d0
      end select

      end function der_I_B


end module bessel_derivatives