module functions

    contains 

double precision function Ireal(n,A)

      use torus_init
      use bessel_functions
      use tools 
      implicit none 

      integer, intent(in)          :: n
      double precision, intent(in) :: A
      double precision,parameter   :: pi = 3.14159265358979323846D00


      ! n mean the power of the term

      if (mod(n,2) == 0) then
        Ireal = dble(factorial2(n-1)) / (2.d0*A)**(n/2) * dsqrt(pi/A)
      else
        Ireal = 0.d0
      end if

end function Ireal

double precision function Icliff(n,m,A)

      use torus_init
      use bessel_functions
      implicit none 

      integer, intent(in)          :: n , m
      double precision, intent(in) :: A
      integer                      :: nm

      ! Combine n and m into a two-digit number
      ! For n=1, m=0 -> nm=10  
      ! For n=0, m=1 -> nm=01
      ! n mean sin and m mean cos

      nm = n*10 + m

      select case (nm)
        case (00)  
            Icliff = Lx * iv_scaled(0, A) 
        case (01)  
            Icliff = Lx * iv_scaled(1, A)
        case (02)  
            Icliff = Lx * ( iv_scaled(1, A) / A  +  iv_scaled(2, A) ) 
            if (A == 0.d0) then 
              Icliff = Lx * 0.5d0  
            end if 
        case (03)  
            Icliff = Lx * ( 3.d0 * iv_scaled(2, A) / A + iv_scaled(3, A) )
            if (A == 0.d0) then 
              Icliff = 0.d0   
            end if 
        case (04)
            Icliff = Lx * (3.d0 + A * A) * iv_scaled(2, A) / (A * A)
            if (A == 0.d0) then 
              Icliff = Lx * 0.375d0
            end if 
             
            
        case (10)  
            Icliff = 0.d0 
        case (11)  
            Icliff = 0.d0 
        case (12)  
            Icliff = 0.d0
        case (13)  
            Icliff = 0.d0
        case (14)  
            Icliff = 0.d0

        case (20)  
            Icliff = Lx *  iv_scaled(1, A) / A 
            if (A == 0.d0) then 
              Icliff = Lx * 0.500000d0
            end if 
        case (21)  
            Icliff = Lx * iv_scaled(2, A) / A 
            if (A == 0.d0) then 
              Icliff = 0.d0
            end if 
        case (22)  
            Icliff = Lx * (A * iv_scaled(1, A) - 3.d0 * iv_scaled(2, A)) / ( A * A )
            if (A == 0.d0) then 
              Icliff = Lx * 0.125d0
            end if 
        case (23)  
            Icliff = Lx * (A * iv_scaled(2, A) - 3.d0 * iv_scaled(3, A)) / ( A * A )
            if (A == 0.d0) then 
              Icliff = 0.d0 
            end if 
        case (24)  
            Icliff = Lx * (  ( 15.d0 + A * A ) * iv_scaled(3, A) - 2.d0 * A * iv_scaled(2, A)  ) / ( A * A * A )
            if (A == 0.d0) then 
              Icliff = Lx * 0.0625d0 
            end if 


        case (30)  
            Icliff = 0.d0 
        case (31)  
            Icliff = 0.d0
        case (32)  
            Icliff = 0.d0
        case (33)
            Icliff = 0.d0
        case (34)  
            Icliff = 0.d0
          
          
        case (40)  
            Icliff = Lx * 3.d0 * iv_scaled(2, A) / ( A * A )
            if (A == 0.d0) then 
              Icliff = Lx * 0.375d0
            end if 
        case (41)  
            Icliff = Lx * 3.d0 * iv_scaled(3, A) / ( A * A )
            if (A == 0.d0) then 
              Icliff = 0.d0
            end if
        case (42)  
            Icliff = Lx * 3.d0 * ( (20.d0 + A * A ) * iv_scaled(2, A) - 5.d0 * A * iv_scaled(1, A) ) / ( A * A * A * A )
            if (A == 0.d0) then 
              Icliff = Lx * 0.0625d0
            end if
        case (43)  
            Icliff = Lx * 3.d0 * ( (30.d0 + A * A ) * iv_scaled(3, A) - 5.d0 * A * iv_scaled(2, A) ) / ( A * A * A * A )
            if (A == 0.d0) then 
              Icliff = 0.d0
            end if
        case (44)  
            Icliff = Lx * 3.d0 * ( A * (35.d0 + A * A ) * iv_scaled(2, A) - 10.d0 * ( 21.d0 + A * A ) * iv_scaled(3, A) ) / ( A * A * A * A * A )
            if (A == 0.d0) then 
              Icliff = Lx * 0.0234375d0 
            end if

        case default
            Icliff = 0.d0
        end select

end function Icliff

end module functions