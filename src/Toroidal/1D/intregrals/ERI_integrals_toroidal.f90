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
      double precision                :: kc1 , kc2 
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

          kc1   = dexp(-(alpha+beta)*(Lx**2)/(2.d0*pi**2))

          mu_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB)))
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

              kc2   = dexp(-(gamma+delta)*(Lx**2)/(2.d0*pi**2))

              nu_x  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD)))
              nu_y  =  gamma+delta
              nu_z  =  gamma+delta

              call bary_center_toroidal(gamma,delta,xc,xd,xq)
              yq    = 0.d0 
              zq    = 0.d0 
               
               nu    = gamma + delta

              const   = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi)*Lx*Lx 

              test = dexp(-(alpha+beta-mu_x)*(Lx**2)/(2.d0*pi**2)) * dexp(-(gamma+delta-nu_x)*(Lx**2)/(2.d0*pi**2))

              if (test < 1e-15) cycle

              xpA     = ax*(xp - xa)
              xpB     = ax*(xp - xb) 
              xqC     = ax*(xq - xc)
              xqD     = ax*(xq - xd)
              phi     = ax*(xp - xq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)

              call integrate_ERI_mod_mod(pattern_id,mu,nu,mu_x,nu_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,value_s)

              value  = value    + const * value_s

              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal

subroutine integrate_ERI_mod_mod(pattern_id,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,result)
      
      use quadpack , only : dqagi
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , p
      double precision, intent(in)       :: q_x , q
      double precision, intent(in)       :: phi
      double precision, intent(in)       :: xpA , xpB 
      double precision, intent(in)       :: xqC , xqD
      double precision, intent(in)       :: xa , xb , xc ,xd ,xp ,xq 
      integer         , intent(in)       :: pattern_id
      

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      integer,parameter                  :: Nmax = 60
      integer,parameter                  :: Nmax_t = 500
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi
      
      
      call dqagi(f_decay, bound, inf, epsabs, epsrel, result, abserr,   &
      &          neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)
        double precision, intent(in) :: t
        double precision             :: ft

         ft =  S1(t)

      end function f_decay


      double precision function S(t) result(sum)

      use gsl_bessel_mod

      implicit none
      double precision, intent(in)         :: t
      double precision                     :: A,  B, C 
      double precision                     :: D, D2
      double precision                     :: const 
      double precision                     :: tol  = 1D-12
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n , small_count
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      double precision                     :: termC
      double precision                     :: sum1 
      COMPLEX(KIND=KIND(1.0D0))            :: term1  , term2 



      A   = 2.d0*p_x/(ax*ax)
      B   = 2.d0*q_x/(ax*ax)
      C   = 2.d0*t*t/(ax*ax)

      D   = 1.d0/(dsqrt(p*q+(p+q)*t*t))
      D2  = D * D 
      
      select case(pattern_id)

      case (0000) ! | s   s   s   s    ( 1 ) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       = const * bessi_scaled(n, C) * termAn * termBn
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do

      case (0001) ! | s   s   s   px   ( 2 ) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
      case (0010) ! | s   s   px  s    ( 5 ) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0011) ! | s   s   px  px   ( 6 ) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0022) ! | s   s   py  py   ( 11) 
        n         = 0
        const     = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (0033) ! | s   s   pz  pz   ( 16) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0100) ! | s   px  s   s    ( 17) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0101) ! | s   px  s   px   ( 18) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (0110) ! | s   px  px  s    ( 21) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0111) ! | s   px  px  px   ( 22) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (0122) ! | s   px  py  py   ( 27) 
        n         = 0
        const     = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0133) ! | s   px  pz  pz   ( 32) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (0202) ! | s   py  s   py   ( 35) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0212) ! | s   py  px  py   ( 39) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0220) ! | s   py  py  s    ( 41) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0221) ! | s   py  py  px   ( 42) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (0303) ! | s   pz  s   pz   ( 52) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0313) ! | s   pz  px  pz   ( 56) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0330) ! | s   pz  pz  s    ( 61) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (0331) ! | s   pz  pz  px   ( 62) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1000) ! | px  s   s   s    ( 65) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1001) ! | px  s   s   px   ( 66) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1010) ! | px  s   px  s    ( 69) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1011) ! | px  s   px  px   ( 70) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
          
        case (1022) ! | px  s   py  py   ( 75) 
        n         = 0
        const     = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1033) ! | px  s   pz  pz   ( 80) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1100) ! | px  px  s   s    ( 81) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1101) ! | px  px  s   px   ( 82) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1110) ! | px  px  px  s    ( 85) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1111) ! | px  px  px  px   ( 86) 
        n         = 0
        const     =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1122) ! | px  px  py  py   ( 91) 
        n         = 0
        const     = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1133) ! | px  px  pz  pz   ( 96) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1202) ! | px  py  s   py   ( 99) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1212) ! | px  py  px  py   ( 103) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1220) ! | px  py  py  s    ( 105) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1221) ! | px  py  py  px   ( 106) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1303) ! | px  pz  s   pz   ( 116) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1313) ! | px  pz  px  pz   ( 120) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1330) ! | px  pz  pz  s    ( 125) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (1331) ! | px  pz  pz  px   ( 126) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (2002) ! | py  s   s   py   ( 131) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2012) ! | py  s   px  py   ( 135) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        
        case (2020) ! | py  s   py  s    ( 137) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2021) ! | py  s   py  px   ( 138) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2102) ! | py  px  s   py   ( 147) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2112) ! | py  px  px  py   ( 151) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2120) ! | py  px  py  s    ( 153) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2121) ! | py  px  py  px   ( 154) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2200) ! | py  py  s   s    ( 161) 
        n         = 0
        const     = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2201) ! | py  py  s   px   ( 162) 
        n         = 0
        const     = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2210) ! | py  py  px  s    ( 165) 
        n         = 0
        const     = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2211) ! | py  py  px  px   ( 166) 
        n         = 0
        const     = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2222) ! | py  py  py  py   ( 171) 
        n         = 0
        const     =   (0.25d0 * ( 1.d0 + 3.d0 * t*t*t*t*D2 ) * D2 * pi * D )  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2233) ! | py  py  pz  pz   ( 176) 
        n         = 0
        const     = ( 0.5d0  * (q+t*t) * D2 * pi * D ) * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (2323) ! | py  pz  py  pz   ( 188) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (2332) ! | py  pz  pz  py   ( 191) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3003) ! | pz  s   s   pz   ( 196) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3013) ! | pz  s   px  pz   ( 200) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3030) ! | pz  s   pz  s    ( 205) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3031) ! | pz  s   pz  px   ( 206) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3103) ! | pz  px  s   pz   ( 212) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3113) ! | pz  px  px  pz   ( 216) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
          
        case (3130) ! | pz  px  pz  s    ( 221) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3131) ! | pz  px  pz  px   ( 222) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3223) ! | pz  py  py  pz   ( 236) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
                
        case (3232) ! | pz  py  pz  py   ( 239) 
        n         = 0
        const     = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3300) ! | pz  pz  s   s    ( 241) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3301) ! | pz  pz  s   px   ( 242) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3310) ! | pz  pz  px  s    ( 245) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3311) ! | pz  pz  px  px   ( 246) 
        n         = 0
        const     =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3322) ! | pz  pz  py  py   ( 251) 
        n         = 0
        const     = ( 0.5d0 * (p+t*t) * D2 * pi * D ) * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
        case (3333) ! | pz  pz  pz  pz   ( 256) 
        n         = 0
        const     =  (pi * D)  *   (0.25d0 * ( 1.d0 + 3.d0 * t*t*t*t*D2 ) * D2 * pi * D )   * exp(A+B-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, A)
        termBn    = bessi_scaled(n, B)
        sum       =  const * bessi_scaled(n, C) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, A)
          termBn  = bessi_scaled(n, B)
          termc   = bessi_scaled(n, C) 
          term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum = sum + real(term1+term2) * const
        end do
        
      case default
        sum = 0.0d0
      end select

end function S


      double precision function S1(t) result(sum)

      use gsl_bessel_mod

      implicit none
      double precision, intent(in)         :: t
      double precision                     :: A,  B, C 
      double precision                     :: D, D2
      double precision                     :: const 
      double precision                     :: tol  = 1D-12
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n , small_count
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      double precision                     :: termC
      double precision                     :: sum1 
      COMPLEX(KIND=KIND(1.0D0))            :: term1  , term2 



      A   = 2.d0*p_x/(ax*ax)
      B   = 2.d0*q_x/(ax*ax)
      C   = 2.d0*t*t/(ax*ax)

      D   = 1.d0/(dsqrt(p*q+(p+q)*t*t))
      D2  = D * D 
      
      select case(pattern_id)

        case (0000) ! | s   s   s   s    ( 1 ) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0001) ! | s   s   s   px   ( 2 ) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0002) ! | s   s   s   py   ( 3 ) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0003) ! | s   s   s   pz   ( 4 ) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0010) ! | s   s   px  s    ( 5 ) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0011) ! | s   s   px  px   ( 6 ) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0012) ! | s   s   px  py   ( 7 ) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0013) ! | s   s   px  pz   ( 8 ) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0020) ! | s   s   py  s    ( 9 ) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0021) ! | s   s   py  px   ( 10) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0022) ! | s   s   py  py   ( 11) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0023) ! | s   s   py  pz   ( 12) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0030) ! | s   s   pz  s    ( 13) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0031) ! | s   s   pz  px   ( 14) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0032) ! | s   s   pz  py   ( 15) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0033) ! | s   s   pz  pz   ( 16) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0100) ! | s   px  s   s    ( 17) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0101) ! | s   px  s   px   ( 18) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0102) ! | s   px  s   py   ( 19) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0103) ! | s   px  s   pz   ( 20) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0110) ! | s   px  px  s    ( 21) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0111) ! | s   px  px  px   ( 22) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0112) ! | s   px  px  py   ( 23) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0113) ! | s   px  px  pz   ( 24) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0120) ! | s   px  py  s    ( 25) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0121) ! | s   px  py  px   ( 26) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0122) ! | s   px  py  py   ( 27) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0123) ! | s   px  py  pz   ( 28) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0130) ! | s   px  pz  s    ( 29) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0131) ! | s   px  pz  px   ( 30) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0132) ! | s   px  pz  py   ( 31) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0133) ! | s   px  pz  pz   ( 32) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0200) ! | s   py  s   s    ( 33) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0201) ! | s   py  s   px   ( 34) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0202) ! | s   py  s   py   ( 35) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0203) ! | s   py  s   pz   ( 36) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0210) ! | s   py  px  s    ( 37) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0211) ! | s   py  px  px   ( 38) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0212) ! | s   py  px  py   ( 39) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0213) ! | s   py  px  pz   ( 40) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0220) ! | s   py  py  s    ( 41) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0221) ! | s   py  py  px   ( 42) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0222) ! | s   py  py  py   ( 43) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0223) ! | s   py  py  pz   ( 44) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0230) ! | s   py  pz  s    ( 45) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0231) ! | s   py  pz  px   ( 46) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0232) ! | s   py  pz  py   ( 47) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0233) ! | s   py  pz  pz   ( 48) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0300) ! | s   pz  s   s    ( 49) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0301) ! | s   pz  s   px   ( 50) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0302) ! | s   pz  s   py   ( 51) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0303) ! | s   pz  s   pz   ( 52) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0310) ! | s   pz  px  s    ( 53) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0311) ! | s   pz  px  px   ( 54) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0312) ! | s   pz  px  py   ( 55) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0313) ! | s   pz  px  pz   ( 56) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0320) ! | s   pz  py  s    ( 57) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0321) ! | s   pz  py  px   ( 58) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0322) ! | s   pz  py  py   ( 59) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0323) ! | s   pz  py  pz   ( 60) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0330) ! | s   pz  pz  s    ( 61) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0331) ! | s   pz  pz  px   ( 62) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0332) ! | s   pz  pz  py   ( 63) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (0333) ! | s   pz  pz  pz   ( 64) 
n           = 0
const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1000) ! | px  s   s   s    ( 65) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1001) ! | px  s   s   px   ( 66) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1002) ! | px  s   s   py   ( 67) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1003) ! | px  s   s   pz   ( 68) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1010) ! | px  s   px  s    ( 69) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1011) ! | px  s   px  px   ( 70) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1012) ! | px  s   px  py   ( 71) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1013) ! | px  s   px  pz   ( 72) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1020) ! | px  s   py  s    ( 73) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1021) ! | px  s   py  px   ( 74) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1022) ! | px  s   py  py   ( 75) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1023) ! | px  s   py  pz   ( 76) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1030) ! | px  s   pz  s    ( 77) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1031) ! | px  s   pz  px   ( 78) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1032) ! | px  s   pz  py   ( 79) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1033) ! | px  s   pz  pz   ( 80) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1100) ! | px  px  s   s    ( 81) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1101) ! | px  px  s   px   ( 82) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1102) ! | px  px  s   py   ( 83) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1103) ! | px  px  s   pz   ( 84) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1110) ! | px  px  px  s    ( 85) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1111) ! | px  px  px  px   ( 86) 
n           = 0
const       =  (pi * D)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1112) ! | px  px  px  py   ( 87) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1113) ! | px  px  px  pz   ( 88) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1120) ! | px  px  py  s    ( 89) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1121) ! | px  px  py  px   ( 90) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1122) ! | px  px  py  py   ( 91) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1123) ! | px  px  py  pz   ( 92) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1130) ! | px  px  pz  s    ( 93) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1131) ! | px  px  pz  px   ( 94) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1132) ! | px  px  pz  py   ( 95) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1133) ! | px  px  pz  pz   ( 96) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (cos(xpA)*cos(xpB)*bessi_scaled(n,A)-cos(ax*(2.d0*xp-xA-xB))*(0.25d0*(bessi_scaled(n-2,A)+2.d0*bessi_scaled(n,A)+bessi_scaled(n+2,A)))+I_dp/A*n*sin(ax*(2.d0*xp-xA-xB))*(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))))/ax/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1200) ! | px  py  s   s    ( 97) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1201) ! | px  py  s   px   ( 98) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1202) ! | px  py  s   py   ( 99) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1203) ! | px  py  s   pz   ( 100) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1210) ! | px  py  px  s    ( 101) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1211) ! | px  py  px  px   ( 102) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1212) ! | px  py  px  py   ( 103) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1213) ! | px  py  px  pz   ( 104) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1220) ! | px  py  py  s    ( 105) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1221) ! | px  py  py  px   ( 106) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1222) ! | px  py  py  py   ( 107) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1223) ! | px  py  py  pz   ( 108) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1230) ! | px  py  pz  s    ( 109) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1231) ! | px  py  pz  px   ( 110) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1232) ! | px  py  pz  py   ( 111) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1233) ! | px  py  pz  pz   ( 112) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1300) ! | px  pz  s   s    ( 113) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1301) ! | px  pz  s   px   ( 114) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1302) ! | px  pz  s   py   ( 115) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1303) ! | px  pz  s   pz   ( 116) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1310) ! | px  pz  px  s    ( 117) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1311) ! | px  pz  px  px   ( 118) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1312) ! | px  pz  px  py   ( 119) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1313) ! | px  pz  px  pz   ( 120) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1320) ! | px  pz  py  s    ( 121) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1321) ! | px  pz  py  px   ( 122) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1322) ! | px  pz  py  py   ( 123) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1323) ! | px  pz  py  pz   ( 124) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1330) ! | px  pz  pz  s    ( 125) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1331) ! | px  pz  pz  px   ( 126) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1332) ! | px  pz  pz  py   ( 127) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (1333) ! | px  pz  pz  pz   ( 128) 
n           = 0
const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpA)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpA))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2000) ! | py  s   s   s    ( 129) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2001) ! | py  s   s   px   ( 130) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2002) ! | py  s   s   py   ( 131) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2003) ! | py  s   s   pz   ( 132) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2010) ! | py  s   px  s    ( 133) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2011) ! | py  s   px  px   ( 134) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2012) ! | py  s   px  py   ( 135) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2013) ! | py  s   px  pz   ( 136) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2020) ! | py  s   py  s    ( 137) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2021) ! | py  s   py  px   ( 138) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2022) ! | py  s   py  py   ( 139) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2023) ! | py  s   py  pz   ( 140) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2030) ! | py  s   pz  s    ( 141) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2031) ! | py  s   pz  px   ( 142) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2032) ! | py  s   pz  py   ( 143) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2033) ! | py  s   pz  pz   ( 144) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2100) ! | py  px  s   s    ( 145) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2101) ! | py  px  s   px   ( 146) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2102) ! | py  px  s   py   ( 147) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2103) ! | py  px  s   pz   ( 148) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2110) ! | py  px  px  s    ( 149) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2111) ! | py  px  px  px   ( 150) 
n           = 0
const       =  (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2112) ! | py  px  px  py   ( 151) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2113) ! | py  px  px  pz   ( 152) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2120) ! | py  px  py  s    ( 153) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2121) ! | py  px  py  px   ( 154) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2122) ! | py  px  py  py   ( 155) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2123) ! | py  px  py  pz   ( 156) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2130) ! | py  px  pz  s    ( 157) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2131) ! | py  px  pz  px   ( 158) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2132) ! | py  px  pz  py   ( 159) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2133) ! | py  px  pz  pz   ( 160) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2200) ! | py  py  s   s    ( 161) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2201) ! | py  py  s   px   ( 162) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2202) ! | py  py  s   py   ( 163) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2203) ! | py  py  s   pz   ( 164) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2210) ! | py  py  px  s    ( 165) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2211) ! | py  py  px  px   ( 166) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2212) ! | py  py  px  py   ( 167) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2213) ! | py  py  px  pz   ( 168) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2220) ! | py  py  py  s    ( 169) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2221) ! | py  py  py  px   ( 170) 
n           = 0
const       =   (0.d0)  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2222) ! | py  py  py  py   ( 171) 
n           = 0
const       =   (0.25d0 * ( 1.d0 + 3.d0 * t*t*t*t*D2 ) * D2 * pi * D )  *  (pi * D)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2223) ! | py  py  py  pz   ( 172) 
n           = 0
const       =   (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2230) ! | py  py  pz  s    ( 173) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2231) ! | py  py  pz  px   ( 174) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2232) ! | py  py  pz  py   ( 175) 
n           = 0
const       =   (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2233) ! | py  py  pz  pz   ( 176) 
n           = 0
const       = ( 0.5d0  * (q+t*t) * D2 * pi * D ) * ( 0.5d0 * (p+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2300) ! | py  pz  s   s    ( 177) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2301) ! | py  pz  s   px   ( 178) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2302) ! | py  pz  s   py   ( 179) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2303) ! | py  pz  s   pz   ( 180) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2310) ! | py  pz  px  s    ( 181) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2311) ! | py  pz  px  px   ( 182) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2312) ! | py  pz  px  py   ( 183) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2313) ! | py  pz  px  pz   ( 184) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2320) ! | py  pz  py  s    ( 185) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2321) ! | py  pz  py  px   ( 186) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2322) ! | py  pz  py  py   ( 187) 
n           = 0
const       =   (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2323) ! | py  pz  py  pz   ( 188) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2330) ! | py  pz  pz  s    ( 189) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2331) ! | py  pz  pz  px   ( 190) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2332) ! | py  pz  pz  py   ( 191) 
n           = 0
const       = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (2333) ! | py  pz  pz  pz   ( 192) 
n           = 0
const       =  (0.d0)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3000) ! | pz  s   s   s    ( 193) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3001) ! | pz  s   s   px   ( 194) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3002) ! | pz  s   s   py   ( 195) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3003) ! | pz  s   s   pz   ( 196) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3010) ! | pz  s   px  s    ( 197) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3011) ! | pz  s   px  px   ( 198) 
n           = 0
const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3012) ! | pz  s   px  py   ( 199) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3013) ! | pz  s   px  pz   ( 200) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3020) ! | pz  s   py  s    ( 201) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3021) ! | pz  s   py  px   ( 202) 
n           = 0
const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3022) ! | pz  s   py  py   ( 203) 
n           = 0
const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3023) ! | pz  s   py  pz   ( 204) 
n           = 0
const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

case (3030) ! | pz  s   pz  s    ( 205) 
n           = 0
const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
termAn      = bessi_scaled(n, A)
termBn      = bessi_scaled(n, B)
sum         = const * bessi_scaled(n, C) * termAn * termBn 
n           = 1
small_count = 0
do
  termAn  = bessi_scaled(n, A)
  termBn  = bessi_scaled(n, B)
  termc   = bessi_scaled(n, C) 
  term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
  term2   = conjg(term1)
  sum     = sum + real(term1+term2) * const
  if (abs(real(term1+term2) * const) < tol) then
    small_count = small_count + 1
    if (small_count >= 2) exit
  else
    small_count = 0
  end if
  n = n + 1
  if (n > Nmax_t) exit
end do

      case (3031) ! | pz  s   pz  px   ( 206) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3032) ! | pz  s   pz  py   ( 207) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3033) ! | pz  s   pz  pz   ( 208) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3100) ! | pz  px  s   s    ( 209) 
      n           = 0
      const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3101) ! | pz  px  s   px   ( 210) 
      n           = 0
      const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3102) ! | pz  px  s   py   ( 211) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3103) ! | pz  px  s   pz   ( 212) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3110) ! | pz  px  px  s    ( 213) 
      n           = 0
      const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3111) ! | pz  px  px  px   ( 214) 
      n           = 0
      const       =  (pi * D)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3112) ! | pz  px  px  py   ( 215) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3113) ! | pz  px  px  pz   ( 216) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3120) ! | pz  px  py  s    ( 217) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3121) ! | pz  px  py  px   ( 218) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3122) ! | pz  px  py  py   ( 219) 
      n           = 0
      const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3123) ! | pz  px  py  pz   ( 220) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3130) ! | pz  px  pz  s    ( 221) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3131) ! | pz  px  pz  px   ( 222) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3132) ! | pz  px  pz  py   ( 223) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3133) ! | pz  px  pz  pz   ( 224) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = (n*I_dp/A*cos(xpB)*bessi_scaled(n,A)+(0.5d0*(bessi_scaled(n-1,A)+bessi_scaled(n+1,A))) * sin(xpB))/ax
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3200) ! | pz  py  s   s    ( 225) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3201) ! | pz  py  s   px   ( 226) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3202) ! | pz  py  s   py   ( 227) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3203) ! | pz  py  s   pz   ( 228) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3210) ! | pz  py  px  s    ( 229) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3211) ! | pz  py  px  px   ( 230) 
      n           = 0
      const       =  (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3212) ! | pz  py  px  py   ( 231) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3213) ! | pz  py  px  pz   ( 232) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3220) ! | pz  py  py  s    ( 233) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3221) ! | pz  py  py  px   ( 234) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3222) ! | pz  py  py  py   ( 235) 
      n           = 0
      const       =   (0.d0)  *  (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3223) ! | pz  py  py  pz   ( 236) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3230) ! | pz  py  pz  s    ( 237) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3231) ! | pz  py  pz  px   ( 238) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3232) ! | pz  py  pz  py   ( 239) 
      n           = 0
      const       = ( 0.5d0 *  t*t * D2 * pi * D ) * ( 0.5d0 *  t*t * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3233) ! | pz  py  pz  pz   ( 240) 
      n           = 0
      const       =  (0.d0)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3300) ! | pz  pz  s   s    ( 241) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3301) ! | pz  pz  s   px   ( 242) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3302) ! | pz  pz  s   py   ( 243) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3303) ! | pz  pz  s   pz   ( 244) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3310) ! | pz  pz  px  s    ( 245) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3311) ! | pz  pz  px  px   ( 246) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (cos(xqC)*cos(xqD)*bessi_scaled(n,B)-cos(ax*(2.d0*xq-xC-xD))*(0.25d0*(bessi_scaled(n-2,B)+2.d0*bessi_scaled(n,B)+bessi_scaled(n+2,B)))-I_dp/B*n*sin(ax*(2.d0*xq-xC-xD))*(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B))))/ax/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3312) ! | pz  pz  px  py   ( 247) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3313) ! | pz  pz  px  pz   ( 248) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqC)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqC))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3320) ! | pz  pz  py  s    ( 249) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3321) ! | pz  pz  py  px   ( 250) 
      n           = 0
      const       =  (0.d0)  * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3322) ! | pz  pz  py  py   ( 251) 
      n           = 0
      const       = ( 0.5d0 * (p+t*t) * D2 * pi * D ) * ( 0.5d0  * (q+t*t) * D2 * pi * D )  * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3323) ! | pz  pz  py  pz   ( 252) 
      n           = 0
      const       =  (0.d0)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3330) ! | pz  pz  pz  s    ( 253) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3331) ! | pz  pz  pz  px   ( 254) 
      n           = 0
      const       =  (pi * D)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = (-n*I_dp/B*cos(xqD)*bessi_scaled(n,B)+(0.5d0*(bessi_scaled(n-1,B)+bessi_scaled(n+1,B)))*sin(xqD))/ax
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3332) ! | pz  pz  pz  py   ( 255) 
      n           = 0
      const       =  (0.d0)  *   (0.d0)   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case (3333) ! | pz  pz  pz  pz   ( 256) 
      n           = 0
      const       =  (pi * D)  *   (0.25d0 * ( 1.d0 + 3.d0 * t*t*t*t*D2 ) * D2 * pi * D )   * exp(A+B-2.d0*(p+q)/(ax*ax))
      termAn      = bessi_scaled(n, A)
      termBn      = bessi_scaled(n, B)
      sum         = const * bessi_scaled(n, C) * termAn * termBn 
      n           = 1
      small_count = 0
      do
        termAn  = bessi_scaled(n, A)
        termBn  = bessi_scaled(n, B)
        termc   = bessi_scaled(n, C) 
        term1   = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        term2   = conjg(term1)
        sum     = sum + real(term1+term2) * const
        if (abs(real(term1+term2) * const) < tol) then
          small_count = small_count + 1
          if (small_count >= 2) exit
        else
          small_count = 0
        end if
        n = n + 1
        if (n > Nmax_t) exit
      end do

      case default
        sum = 0.0d0
      end select

      end function S1



end subroutine integrate_ERI_mod_mod