subroutine ERI_integral_4_function_toroidal_3D(one,two,three,four,value)
      
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
      double precision                :: der 
      double precision                :: mu_x, mu_y , mu_z
      double precision                :: nu_x, nu_y , nu_z
      double precision                :: mu  , nu 
      double precision                :: xpA , xpB , xqC , xqD , phix
      double precision                :: ypA , ypB , yqC , yqD , phiy
      double precision                :: zpA , zpB , zqC , zqD , phiz 
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

          mu_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB)))
          mu_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(YAB)))
          mu_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(ZAB)))

          call bary_center_toroidal(alpha,beta,xa,xb,xp)
          call bary_center_toroidal(alpha,beta,ya,yb,yp)
          call bary_center_toroidal(alpha,beta,za,zb,zp)

          mu = alpha+beta 
          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital


              nu_x  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD)))
              nu_y  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ay*(YCD)))
              nu_z  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(az*(ZCD)))

              call bary_center_toroidal(gamma,delta,xc,xd,xq)
              call bary_center_toroidal(gamma,delta,yc,yd,yq)
              call bary_center_toroidal(gamma,delta,zc,zd,zq)
               
              nu    = gamma + delta

              const   = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi) * (Lx*Lx) * (Ly*Ly) * (Lz*Lz)  

              !test = dexp(-(alpha+beta-mu_x)*(Lx**2)/(2.d0*pi**2)) * dexp(-(gamma+delta-nu_x)*(Lx**2)/(2.d0*pi**2))

              !if (test < 1e-12) cycle

              xpA     = ax*(xp - xa) ; ypA     = ay*(yp - ya) ; zpA     = az*(zp - za)
              xpB     = ax*(xp - xb) ; ypB     = ay*(yp - yb) ; zpB     = az*(zp - zb) 
              xqC     = ax*(xq - xc) ; yqC     = ay*(yq - yc) ; zqC     = az*(zq - zc)
              xqD     = ax*(xq - xd) ; yqD     = ay*(yq - yd) ; zqD     = az*(zq - zd)
              phix    = ax*(xp - xq) ; phiy    = ay*(yp - yq) ; phiz    = az*(zp - zq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)
 
              call integrate_ERI_3D(pattern_id,mu,nu,mu_x,nu_x,phix,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                                     mu_y,nu_y,phiy,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                                     mu_z,nu_z,phiz,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,value_s)

              value  = value    + const * value_s

            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_3D

subroutine integrate_ERI_3D(pattern_id,p,q,p_x,q_x,phix,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                           p_y,q_y,phiy,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                           p_z,q_z,phiz,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,result)
      
      use quadpack , only : dqagi
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , q_x 
      double precision, intent(in)       :: p_y , q_y
      double precision, intent(in)       :: p_z , q_z 
      double precision, intent(in)       :: p   , q 
      double precision, intent(in)       :: phix , phiy , phiz
      double precision, intent(in)       :: xpA , xpB , xqC , xqD
      double precision, intent(in)       :: ypA , ypB , yqC , yqD
      double precision, intent(in)       :: zpA , zpB , zqC , zqD
      double precision, intent(in)       :: xa , xb , xc ,xd ,xp ,xq 
      double precision, intent(in)       :: ya , yb , yc ,yd ,yp ,yq
      double precision, intent(in)       :: za , zb , zc ,zd ,zp ,zq  
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
      integer,parameter                  :: Nmax = 150
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi
      
      
      call dqagi(f_decay , bound, inf, epsabs, epsrel, result, abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)
        double precision, intent(in) :: t
        double precision             :: ft

         ft =  S(t)

      end function f_decay

      double precision function S(t) result(sum)

      use gsl_bessel_mod

      implicit none
      double precision, intent(in)         :: t
      double precision                     :: AAx,  BBx, CCx
      double precision                     :: AAy,  BBy, CCy
      double precision                     :: AAz,  BBz, CCz
      double precision                     :: D, D2
      double precision                     :: const 
      double precision                     :: tol  = 1D-12
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n  
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      double precision                     :: termC 
      double precision                     :: sum1, sum2, sum3
      COMPLEX(KIND=KIND(1.0D0))            :: term1  , term2 


      AAx   = 2.d0*p_x/(ax*ax)
      BBx   = 2.d0*q_x/(ax*ax)
      CCx   = 2.d0*t*t/(ax*ax)

      AAy   = 2.d0*p_y/(ay*ay)
      BBy   = 2.d0*q_y/(ay*ay)
      CCy   = 2.d0*t*t/(ay*ay)

      AAz   = 2.d0*p_z/(az*az)
      BBz   = 2.d0*q_z/(az*az)
      CCz   = 2.d0*t*t/(az*az)
      
      select case(pattern_id)
        
      case (0000) ! | s   s   s   s    ( 1 ) 

        sum1      = 0.d0
        sum2      = 0.d0
        sum3      = 0.d0 
        sum       = 0.d0

        n         = 0
        const     = exp(AAx+BBx-2.d0*(p+q)/(ax*ax))
        termAn    = bessi_scaled(n, AAx)
        termBn    = bessi_scaled(n, BBx)
        sum1      = const * bessi_scaled(n, CCx) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, AAx)
          termBn  = bessi_scaled(n, BBx)
          termc   = bessi_scaled(n, CCx) 
          term1   = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum1 = sum1 + real(term1+term2) * const
        end do
        
        n         = 0
        const     = exp(AAy+BBy-2.d0*(p+q)/(ay*ay))
        termAn    = bessi_scaled(n, AAy)
        termBn    = bessi_scaled(n, BBy)
        sum2      = const * bessi_scaled(n, CCy) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, AAy)
          termBn  = bessi_scaled(n, BBy)
          termc   = bessi_scaled(n, CCy) 
          term1   = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum2 = sum2 + real(term1+term2) * const
        end do

        n         = 0
        const     = exp(AAz+BBz-2.d0*(p+q)/(az*az))
        termAn    = bessi_scaled(n, AAz)
        termBn    = bessi_scaled(n, BBz)
        sum3      = const * bessi_scaled(n, CCz) * termAn * termBn 
        do n      = 1 , Nmax
          termAn  = bessi_scaled(n, AAz)
          termBn  = bessi_scaled(n, BBz)
          termc   = bessi_scaled(n, CCz) 
          term1   = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
          term2   = conjg(term1)
          if (abs(term1) < tol) exit
          if (abs(term2) < tol) exit
          sum3 = sum3 + real(term1+term2) * const
        end do

        sum = sum1 * sum2 * sum3 

      case default
        sum = 0.0d0
      end select

      end function S

end subroutine integrate_ERI_3D