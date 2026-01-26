subroutine ERI_integral_4_function_toroidal(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI
      use HeavisideModule

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 

      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision,parameter      :: pi     = 3.14159265358979323846D00
      
      double precision                :: xa  , ya , za 
      double precision                :: xb  , yb , zb
      double precision                :: xc  , yc , zc
      double precision                :: xd  , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: zAB , zCD
      double precision                :: xp  , xq 
      double precision                :: yp  , yq
      double precision                :: zp  , zq
      double precision                :: c1  , c2 , c3 , c4 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: const  
      double precision                :: value_s 
      double precision                :: mu_x, mu_y , mu_z
      double precision                :: nu_x, nu_y , nu_z
      double precision                :: mu  , nu 
      double precision                :: xpA , xpB , xqC , xqD , phi 
      double precision                :: test 
      integer                         :: pattern_id, encode_orbital_pattern



      !-----------------------------------------------------------------!

      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital

          mu_x  = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB))))
          mu_y  = alpha + beta
          mu_z  = alpha + beta

          call bary_center_toroidal(alpha,beta,xa,xb,xp)
          yp    = 0.d0
          zp    = 0.d0

          mu = alpha+beta 
          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital

              nu_x  = dsqrt(dabs(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD))))
              nu_y  =  gamma+delta
              nu_z  =  gamma+delta

              call bary_center_toroidal(gamma,delta,xc,xd,xq)
              yq     = 0.d0 
              zq     = 0.d0

              nu     = gamma + delta

              const  = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi)*Lx*Lx

              !test = dexp(-(alpha+beta-mu_x)*(Lx**2)/(2.d0*pi**2)) * dexp(-(gamma+delta-nu_x)*(Lx**2)/(2.d0*pi**2))

              xpA     = ax*(xp - xa)
              xpB     = ax*(xp - xb)
              xqC     = ax*(xq - xc)
              xqD     = ax*(xq - xd)
              phi     = ax*(xp - xq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)

              call integrate_ERI_sum(pattern_id,mu,nu,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,value_s)
              !call integrate_ERI_integral(pattern_id,px_count,mu,nu,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,value_s)
              !call integrate_ERI_integral_mod(pattern_id,mu,nu,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,value_s)

              value  = value    + const * value_s  
              
              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal