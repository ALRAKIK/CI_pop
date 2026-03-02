subroutine integrate_ERI_sum(pattern_id,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq &
     &                                                                 ,ya,yb,yc,yd,yp,yq &
     &                                                                 ,za,zb,zc,zd,zp,zq,result)
      
      use quadpack , only : dqagi , dqags
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod
      use bessel_functions
      use files
      use table_lookup_module

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , p
      double precision, intent(in)       :: q_x , q
      double precision, intent(in)       :: phi
      double precision, intent(in)       :: xpA , xpB 
      double precision, intent(in)       :: xqC , xqD
      double precision, intent(in)       :: xa, xb, xc, xd, xp, xq 
      double precision, intent(in)       :: ya, yb, yc, yd, yp, yq
      double precision, intent(in)       :: za, zb, zc, zd, zp, zq
      integer         , intent(in)       :: pattern_id
      

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0d-10 , epsrel = 1.0d-8
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      
      integer                            :: Nmax
      !double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi

      double precision                   :: inv_ax 
      double precision                   :: inv_ax2

      double precision                   :: cxpa   , sxpa
      double precision                   :: cxpb   , sxpb
      double precision                   :: cxqc   , sxqc
      double precision                   :: cxqd   , sxqd
      double precision                   :: c2xpab , s2xpab
      double precision                   :: c2xqcd , s2xqcd
      double precision                   :: ypq    , zpq 
      double precision                   :: ypa    , ypb 
      double precision                   :: yqc    , yqd 
      double precision                   :: zpa    , zpb 
      double precision                   :: zqc    , zqd 

      double precision                   :: A , B , C
      !integer                            :: Peak

      inv_ax  = 1.d0/ax 
      inv_ax2 = inv_ax * inv_ax 

      cxpa  = dcos(xpA)
      sxpa  = dsin(xpA)

      cxpb  = dcos(xpB)
      sxpb  = dsin(xpB)
      
      cxqc  = dcos(xqC)
      sxqc  = dsin(xqC)

      cxqd  = dcos(xqD)
      sxqd  = dsin(xqD)

      c2xpab = dcos(ax*(2.d0*xp-xA-xB))
      s2xpab = dsin(ax*(2.d0*xp-xA-xB))

      c2xqcd = dcos(ax*(2.d0*xq-xC-xD))
      s2xqcd = dsin(ax*(2.d0*xq-xC-xD))


      ! /////////////////////////////////////////////////////////////// !
      ! Real Gaussian ! 

      ypq    = yp - yq 

      ypa    = yp - ya 
      ypb    = yp - yb 

      yqc    = yq - yc
      yqd    = yq - yd

      zpq    = zp - zq

      zpa    = zp - za 
      zpb    = zp - zb 

      zqc    = zq - zc
      zqd    = zq - zd

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      
      A   = 2.d0  * p_x   * inv_ax2
      B   = 2.d0  * q_x   * inv_ax2


      !call dqagi(f_decay, bound, inf, epsabs, epsrel, result, abserr,   &
      !&          neval, ier,Limit,Lenw,Last,Iwork,Work)

      call dqags(transformed_integrand, 0.0d0, 1.0d0, epsabs, epsrel, &
                 result, abserr, neval, ier, limit, lenw, last, iwork, work)


      if (ier > 2) then
        write(*,'(A,I8,A)') 'Error code = ', ier
        write(outfile,'(A,I8,A)') 'you have an integral did not converge '
        !stop
      end if

      contains

      function f_decay(t) result(ft)
        double precision, intent(in) :: t
        double precision             :: ft
         ft =  S(t)
      end function f_decay

      ! Transformed integrand for u in [0,1)

      function transformed_integrand(u) result(y)
        double precision, intent(in) :: u
        double precision :: y, t
        
        if (u >= 1.0d0) then
          y = 0.0d0
          return
        endif
        
        ! Transformation: t = u/(1-u)
        t = u / (1.0d0 - u)
        
        ! Jacobian: dt/du = 1/(1-u)²
        y = S(t) / ((1.0d0 - u) * (1.0d0 - u))
        
        ! Check for numerical issues
        if (.not. ieee_is_finite(y)) then
          y = 0.0d0
        endif
      end function transformed_integrand


      double precision function S(t) result(sum)

      use gsl_bessel_mod

      implicit none
      double precision, intent(in)         :: t
      double precision                     :: D, D2
      double precision                     :: const 
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n 
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      double precision                     :: termC
      COMPLEX(KIND=KIND(1.0D0))            :: term
      double precision                     :: t2 , t3 , t4 
      double precision,parameter           :: eps = 2.22d-16
      COMPLEX(KIND=KIND(1.0D0))            :: current_term 
      COMPLEX(KIND=KIND(1.0D0))            :: expo_term


      double precision                     :: eta_t   , eta_p  , eta_q 
      double precision                     :: eta_t2  , eta_p2 , eta_q2 
      double precision                     :: eta_p_p , eta_q_p

      ! /////////////////////////////////////////////////////////////// !
      double precision                     :: y1D , y1D2
      double precision                     :: z1D , z1D2 
      double precision                     :: y2S , y2S2 , y2S3 , y2S4
      double precision                     :: z2S , z2S2 , z2S3 , z2S4 
      ! /////////////////////////////////////////////////////////////// !
      double precision                     :: ypq2 , ypq3 , ypq4
      double precision                     :: zpq2 , zpq3 , zpq4
      ! /////////////////////////////////////////////////////////////// !


      ! --------------------------------------------------------------- !
      ! backward ! 
      double precision                     :: I_A(0:10000)
      double precision                     :: I_B(0:10000)
      double precision                     :: I_C(0:10000)
      ! forward !
      !double precision                     :: I_A_prev, I_A_curr, I_A_next
      !double precision                     :: I_B_prev, I_B_curr, I_B_next
      !double precision                     :: I_C_prev, I_C_curr, I_C_next
      ! --------------------------------------------------------------- !

      t2 = t  * t 
      t3 = t2 * t 
      t4 = t3 * t 

      C   = 2.d0  * t * t * inv_ax2

      D   = 1.d0/(dsqrt(p*q+(p+q)*t2))

      D2  = D * D

      ! /////////////////////////////////////////////////////////////// !

      eta_t   = t2 / (p+t2)    ; eta_t2 = eta_t * eta_t  
      eta_p   =  p*t2  * D2    ; eta_p2 = eta_p * eta_p 
      eta_q   =  q*t2  * D2    ; eta_q2 = eta_q * eta_q
      eta_p_p = (p+t2) * D2
      eta_q_p = (q+t2) * D2

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      ypq2 = ypq  * ypq    ; zpq2 = zpq  * zpq 
      ypq3 = ypq2 * ypq    ; zpq3 = zpq2 * zpq
      ypq4 = ypq3 * ypq    ; zpq4 = zpq3 * zpq 
      
      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !

      y1D  = 0.d0 ; y1D2 = 0.5d0 / (p+t2)
      z1D  = 0.d0 ; z1D2 = 0.5d0 / (p+t2)
      
      y2S  = 0.d0 
      y2S2 = 0.5d0  * eta_p_p
      y2S3 = 0.d0
      y2S4 = 0.75d0 * eta_p_p * eta_p_p

      z2S  = 0.d0
      z2S2 = 0.5d0  * eta_p_p
      z2S3 = 0.d0
      z2S4 = 0.75d0 * eta_p_p * eta_p_p 

      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ !


      termAn = 0.d0 
      termBn = 0.d0 
      termC  = 0.d0

      !               the peak of the sum              !
      !Peak         = ceiling(min(A,B,C))
      !         the maximum of terms in the sum        !
      !Nmax        = Peak+10
      Nmax         = get_Nmax(A, B, C, Phi) + 10
      !                    Phase term                  !
      expo_term    = exp(I_dp*phi)
      !                   Initial term                 !
      current_term = expo_term

      

      call bessel_I_scaled_backward(Nmax, A, I_A)
      call bessel_I_scaled_backward(Nmax, B, I_B)
      call bessel_I_scaled_backward(Nmax, C, I_C)
      
      select case(pattern_id)
      
      ! /////////////////////////////////////////////////////////////// !

      case (0000) ! | s   s   s   s    ( 1 ) 

      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0001) ! | s   s   s   px   ( 2 ) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do


      case (0002) ! | s   s   s   py   ( 3 ) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0003) ! | s   s   s   pz   ( 4 ) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0010) ! | s   s   px  s    ( 5 ) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0011) ! | s   s   px  px   ( 6 ) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0012) ! | s   s   px  py   ( 7 ) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0013) ! | s   s   px  pz   ( 8 ) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0020) ! | s   s   py  s    ( 9 ) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0021) ! | s   s   py  px   ( 10) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0022) ! | s   s   py  py   ( 11) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq2 + 2*eta_p*y2S*ypq + eta_p*ypq*yqc + eta_p*ypq*yqd + y2S2 + y2S*yqc + y2S*yqd + yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0023) ! | s   s   py  pz   ( 12) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqd + eta_p*yqc*zpq + y2S*z2S + y2S*zqd + yqc*z2S + yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0030) ! | s   s   pz  s    ( 13) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0031) ! | s   s   pz  px   ( 14) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0032) ! | s   s   pz  py   ( 15) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqc + eta_p*yqd*zpq + y2S*z2S + y2S*zqc + yqd*z2S + yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0033) ! | s   s   pz  pz   ( 16) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*zpq2 + 2*eta_p*z2S*zpq + eta_p*zpq*zqc + eta_p*zpq*zqd + z2S2 + z2S*zqc + z2S*zqd + zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0100) ! | s   px  s   s    ( 17) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0101) ! | s   px  s   px   ( 18) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0102) ! | s   px  s   py   ( 19) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0103) ! | s   px  s   pz   ( 20) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0110) ! | s   px  px  s    ( 21) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0111) ! | s   px  px  px   ( 22) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0112) ! | s   px  px  py   ( 23) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0113) ! | s   px  px  pz   ( 24) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0120) ! | s   px  py  s    ( 25) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0121) ! | s   px  py  px   ( 26) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0122) ! | s   px  py  py   ( 27) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq2 + 2*eta_p*y2S*ypq + eta_p*ypq*yqc + eta_p*ypq*yqd + y2S2 + y2S*yqc + y2S*yqd + yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0123) ! | s   px  py  pz   ( 28) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqd + eta_p*yqc*zpq + y2S*z2S + y2S*zqd + yqc*z2S + yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0130) ! | s   px  pz  s    ( 29) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0131) ! | s   px  pz  px   ( 30) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0132) ! | s   px  pz  py   ( 31) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqc + eta_p*yqd*zpq + y2S*z2S + y2S*zqc + yqd*z2S + yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0133) ! | s   px  pz  pz   ( 32) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*zpq2 + 2*eta_p*z2S*zpq + eta_p*zpq*zqc + eta_p*zpq*zqd + z2S2 + z2S*zqc + z2S*zqd + zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0200) ! | s   py  s   s    ( 33) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0201) ! | s   py  s   px   ( 34) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0202) ! | s   py  s   py   ( 35) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypb + ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0203) ! | s   py  s   pz   ( 36) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypb*z2S + ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0210) ! | s   py  px  s    ( 37) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0211) ! | s   py  px  px   ( 38) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0212) ! | s   py  px  py   ( 39) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypb + ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0213) ! | s   py  px  pz   ( 40) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypb*z2S + ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0220) ! | s   py  py  s    ( 41) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypb + ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0221) ! | s   py  py  px   ( 42) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypb + ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0222) ! | s   py  py  py   ( 43) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq3 + eta_p2*eta_t*y2S*ypq2 + eta_p2*y1D*ypq2 + eta_p2*ypb*ypq2 - 2*eta_p*eta_q*y2S*ypq2 - eta_p*eta_q*ypq2*yqc - eta_p*eta_q*ypq2*yqd + 2*eta_p*eta_t*y2S2*ypq + eta_p*eta_t*y2S*ypq*yqc + eta_p*eta_t*y2S*ypq*yqd + 2*eta_p*y1D*y2S*ypq + eta_p*y1D*ypq*yqc + eta_p*y1D*ypq*yqd + 2*eta_p*y2S*ypb*ypq + eta_p*ypb*ypq*yqc + eta_p*ypb*ypq*yqd - eta_q*y2S2*ypq - eta_q*y2S*ypq*yqc - eta_q*y2S*ypq*yqd - eta_q*ypq*yqc*yqd + eta_t*y2S3 + eta_t*y2S2*yqc + eta_t*y2S2*yqd + eta_t*y2S*yqc*yqd + y1D*y2S2 + y1D*y2S*yqc + y1D*y2S*yqd + y1D*yqc*yqd + y2S2*ypb + y2S*ypb*yqc + y2S*ypb*yqd + ypb*yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0223) ! | s   py  py  pz   ( 44) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypb*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqd - eta_p*eta_q*ypq*yqc*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqd + eta_p*eta_t*y2S*yqc*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqd + eta_p*y1D*yqc*zpq + eta_p*y2S*ypb*zpq + eta_p*ypb*ypq*z2S + eta_p*ypb*ypq*zqd + eta_p*ypb*yqc*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqd - eta_q*ypq*yqc*z2S - eta_q*ypq*yqc*zqd + eta_t*y2S2*z2S + eta_t*y2S2*zqd + eta_t*y2S*yqc*z2S + eta_t*y2S*yqc*zqd + y1D*y2S*z2S + y1D*y2S*zqd + y1D*yqc*z2S + y1D*yqc*zqd + y2S*ypb*z2S + y2S*ypb*zqd + ypb*yqc*z2S + ypb*yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0230) ! | s   py  pz  s    ( 45) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypb*z2S + ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0231) ! | s   py  pz  px   ( 46) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypb*z2S + ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0232) ! | s   py  pz  py   ( 47) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypb*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqc - eta_p*eta_q*ypq*yqd*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqc + eta_p*eta_t*y2S*yqd*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqc + eta_p*y1D*yqd*zpq + eta_p*y2S*ypb*zpq + eta_p*ypb*ypq*z2S + eta_p*ypb*ypq*zqc + eta_p*ypb*yqd*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqc - eta_q*ypq*yqd*z2S - eta_q*ypq*yqd*zqc + eta_t*y2S2*z2S + eta_t*y2S2*zqc + eta_t*y2S*yqd*z2S + eta_t*y2S*yqd*zqc + y1D*y2S*z2S + y1D*y2S*zqc + y1D*yqd*z2S + y1D*yqd*zqc + y2S*ypb*z2S + y2S*ypb*zqc + ypb*yqd*z2S + ypb*yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0233) ! | s   py  pz  pz   ( 48) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*y2S*zpq2 + eta_p2*y1D*zpq2 + eta_p2*ypb*zpq2 - 2*eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*ypq*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*y2S*zpq*zqc + eta_p*eta_t*y2S*zpq*zqd + 2*eta_p*y1D*z2S*zpq + eta_p*y1D*zpq*zqc + eta_p*y1D*zpq*zqd + 2*eta_p*ypb*z2S*zpq + eta_p*ypb*zpq*zqc + eta_p*ypb*zpq*zqd - eta_q*ypq*z2S2 - eta_q*ypq*z2S*zqc - eta_q*ypq*z2S*zqd - eta_q*ypq*zqc*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*y2S*z2S*zqd + eta_t*y2S*zqc*zqd + y1D*z2S2 + y1D*z2S*zqc + y1D*z2S*zqd + y1D*zqc*zqd + ypb*z2S2 + ypb*z2S*zqc + ypb*z2S*zqd + ypb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0300) ! | s   pz  s   s    ( 49) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0301) ! | s   pz  s   px   ( 50) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0302) ! | s   pz  s   py   ( 51) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpb + yqd*z1D + yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0303) ! | s   pz  s   pz   ( 52) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpb + zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0310) ! | s   pz  px  s    ( 53) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0311) ! | s   pz  px  px   ( 54) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0312) ! | s   pz  px  py   ( 55) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpb + yqd*z1D + yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0313) ! | s   pz  px  pz   ( 56) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpb + zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0320) ! | s   pz  py  s    ( 57) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpb + yqc*z1D + yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0321) ! | s   pz  py  px   ( 58) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpb + yqc*z1D + yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0322) ! | s   pz  py  py   ( 59) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*ypq2*z2S + eta_p2*ypq2*z1D + eta_p2*ypq2*zpb - 2*eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq*yqc*zpq - eta_p*eta_q*ypq*yqd*zpq + 2*eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*ypq*yqc*z2S + eta_p*eta_t*ypq*yqd*z2S + 2*eta_p*y2S*ypq*z1D + 2*eta_p*y2S*ypq*zpb + eta_p*ypq*yqc*z1D + eta_p*ypq*yqc*zpb + eta_p*ypq*yqd*z1D + eta_p*ypq*yqd*zpb - eta_q*y2S2*zpq - eta_q*y2S*yqc*zpq - eta_q*y2S*yqd*zpq - eta_q*yqc*yqd*zpq + eta_t*y2S2*z2S + eta_t*y2S*yqc*z2S + eta_t*y2S*yqd*z2S + eta_t*yqc*yqd*z2S + y2S2*z1D + y2S2*zpb + y2S*yqc*z1D + y2S*yqc*zpb + y2S*yqd*z1D + y2S*yqd*zpb + yqc*yqd*z1D + yqc*yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0323) ! | s   pz  py  pz   ( 60) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpb*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqd - eta_p*eta_q*yqc*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqd + eta_p*eta_t*yqc*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpb*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqd + eta_p*ypq*z2S*zpb + eta_p*ypq*zpb*zqd + eta_p*yqc*z1D*zpq + eta_p*yqc*zpb*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqd - eta_q*yqc*z2S*zpq - eta_q*yqc*zpq*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqd + eta_t*yqc*z2S2 + eta_t*yqc*z2S*zqd + y2S*z1D*z2S + y2S*z1D*zqd + y2S*z2S*zpb + y2S*zpb*zqd + yqc*z1D*z2S + yqc*z1D*zqd + yqc*z2S*zpb + yqc*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0330) ! | s   pz  pz  s    ( 61) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpb + zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0331) ! | s   pz  pz  px   ( 62) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpb + zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0332) ! | s   pz  pz  py   ( 63) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpb*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*yqd*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqc + eta_p*eta_t*yqd*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpb*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqc + eta_p*ypq*z2S*zpb + eta_p*ypq*zpb*zqc + eta_p*yqd*z1D*zpq + eta_p*yqd*zpb*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqc - eta_q*yqd*z2S*zpq - eta_q*yqd*zpq*zqc + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*yqd*z2S2 + eta_t*yqd*z2S*zqc + y2S*z1D*z2S + y2S*z1D*zqc + y2S*z2S*zpb + y2S*zpb*zqc + yqd*z1D*z2S + yqd*z1D*zqc + yqd*z2S*zpb + yqd*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (0333) ! | s   pz  pz  pz   ( 64) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*zpq3 + eta_p2*eta_t*z2S*zpq2 + eta_p2*z1D*zpq2 + eta_p2*zpb*zpq2 - 2*eta_p*eta_q*z2S*zpq2 - eta_p*eta_q*zpq2*zqc - eta_p*eta_q*zpq2*zqd + 2*eta_p*eta_t*z2S2*zpq + eta_p*eta_t*z2S*zpq*zqc + eta_p*eta_t*z2S*zpq*zqd + 2*eta_p*z1D*z2S*zpq + eta_p*z1D*zpq*zqc + eta_p*z1D*zpq*zqd + 2*eta_p*z2S*zpb*zpq + eta_p*zpb*zpq*zqc + eta_p*zpb*zpq*zqd - eta_q*z2S2*zpq - eta_q*z2S*zpq*zqc - eta_q*z2S*zpq*zqd - eta_q*zpq*zqc*zqd + eta_t*z2S3 + eta_t*z2S2*zqc + eta_t*z2S2*zqd + eta_t*z2S*zqc*zqd + z1D*z2S2 + z1D*z2S*zqc + z1D*z2S*zqd + z1D*zqc*zqd + z2S2*zpb + z2S*zpb*zqc + z2S*zpb*zqd + zpb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1000) ! | px  s   s   s    ( 65) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1001) ! | px  s   s   px   ( 66) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1002) ! | px  s   s   py   ( 67) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1003) ! | px  s   s   pz   ( 68) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1010) ! | px  s   px  s    ( 69) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1011) ! | px  s   px  px   ( 70) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1012) ! | px  s   px  py   ( 71) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1013) ! | px  s   px  pz   ( 72) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1020) ! | px  s   py  s    ( 73) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1021) ! | px  s   py  px   ( 74) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1022) ! | px  s   py  py   ( 75) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq2 + 2*eta_p*y2S*ypq + eta_p*ypq*yqc + eta_p*ypq*yqd + y2S2 + y2S*yqc + y2S*yqd + yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1023) ! | px  s   py  pz   ( 76) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqd + eta_p*yqc*zpq + y2S*z2S + y2S*zqd + yqc*z2S + yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1030) ! | px  s   pz  s    ( 77) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1031) ! | px  s   pz  px   ( 78) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1032) ! | px  s   pz  py   ( 79) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqc + eta_p*yqd*zpq + y2S*z2S + y2S*zqc + yqd*z2S + yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1033) ! | px  s   pz  pz   ( 80) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*zpq2 + 2*eta_p*z2S*zpq + eta_p*zpq*zqc + eta_p*zpq*zqd + z2S2 + z2S*zqc + z2S*zqd + zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1100) ! | px  px  s   s    ( 81) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1101) ! | px  px  s   px   ( 82) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1102) ! | px  px  s   py   ( 83) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1103) ! | px  px  s   pz   ( 84) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1110) ! | px  px  px  s    ( 85) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1111) ! | px  px  px  px   ( 86) 
      
      n         = 0
      const     = pi2 * D2 * 1.d0 * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1112) ! | px  px  px  py   ( 87) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqd + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1113) ! | px  px  px  pz   ( 88) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqd + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1120) ! | px  px  py  s    ( 89) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1121) ! | px  px  py  px   ( 90) 
      
      n         = 0
      const     = pi2 * D2 * (y2S + yqc + eta_p * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1122) ! | px  px  py  py   ( 91) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq2 + 2*eta_p*y2S*ypq + eta_p*ypq*yqc + eta_p*ypq*yqd + y2S2 + y2S*yqc + y2S*yqd + yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1123) ! | px  px  py  pz   ( 92) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqd + eta_p*yqc*zpq + y2S*z2S + y2S*zqd + yqc*z2S + yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1130) ! | px  px  pz  s    ( 93) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1131) ! | px  px  pz  px   ( 94) 
      
      n         = 0
      const     = pi2 * D2 * (z2S + zqc + eta_p * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1132) ! | px  px  pz  py   ( 95) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*ypq*zpq + eta_p*y2S*zpq + eta_p*ypq*z2S + eta_p*ypq*zqc + eta_p*yqd*zpq + y2S*z2S + y2S*zqc + yqd*z2S + yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1133) ! | px  px  pz  pz   ( 96) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*zpq2 + 2*eta_p*z2S*zpq + eta_p*zpq*zqc + eta_p*zpq*zqd + z2S2 + z2S*zqc + z2S*zqd + zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax2 * (cxpa * cxpb * I_A(0) - c2xpab * (0.25d0 * (I_A(abs(0-2)) + 2.d0 * I_A(0) + I_A(abs(0+2)) ) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        if (dabs(A) < 1.d-300) then
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) )
        else
          termAn  = inv_ax2 * (cxpa * cxpb * I_A(n) - c2xpab * (0.25d0 * (I_A(abs(n-2)) + 2.d0 * I_A(n) + I_A(abs(n+2)) ) ) + I_dp/A * n * s2xpab * (0.5d0 * ( I_A(abs(n-1)) + I_A(abs(n+1)) ) ) )
        end if
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1200) ! | px  py  s   s    ( 97) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1201) ! | px  py  s   px   ( 98) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1202) ! | px  py  s   py   ( 99) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypb + ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1203) ! | px  py  s   pz   ( 100) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypb*z2S + ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1210) ! | px  py  px  s    ( 101) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1211) ! | px  py  px  px   ( 102) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypb + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1212) ! | px  py  px  py   ( 103) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypb + ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1213) ! | px  py  px  pz   ( 104) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypb*z2S + ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1220) ! | px  py  py  s    ( 105) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypb + ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1221) ! | px  py  py  px   ( 106) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypb*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypb + ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1222) ! | px  py  py  py   ( 107) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq3 + eta_p2*eta_t*y2S*ypq2 + eta_p2*y1D*ypq2 + eta_p2*ypb*ypq2 - 2*eta_p*eta_q*y2S*ypq2 - eta_p*eta_q*ypq2*yqc - eta_p*eta_q*ypq2*yqd + 2*eta_p*eta_t*y2S2*ypq + eta_p*eta_t*y2S*ypq*yqc + eta_p*eta_t*y2S*ypq*yqd + 2*eta_p*y1D*y2S*ypq + eta_p*y1D*ypq*yqc + eta_p*y1D*ypq*yqd + 2*eta_p*y2S*ypb*ypq + eta_p*ypb*ypq*yqc + eta_p*ypb*ypq*yqd - eta_q*y2S2*ypq - eta_q*y2S*ypq*yqc - eta_q*y2S*ypq*yqd - eta_q*ypq*yqc*yqd + eta_t*y2S3 + eta_t*y2S2*yqc + eta_t*y2S2*yqd + eta_t*y2S*yqc*yqd + y1D*y2S2 + y1D*y2S*yqc + y1D*y2S*yqd + y1D*yqc*yqd + y2S2*ypb + y2S*ypb*yqc + y2S*ypb*yqd + ypb*yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1223) ! | px  py  py  pz   ( 108) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypb*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqd - eta_p*eta_q*ypq*yqc*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqd + eta_p*eta_t*y2S*yqc*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqd + eta_p*y1D*yqc*zpq + eta_p*y2S*ypb*zpq + eta_p*ypb*ypq*z2S + eta_p*ypb*ypq*zqd + eta_p*ypb*yqc*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqd - eta_q*ypq*yqc*z2S - eta_q*ypq*yqc*zqd + eta_t*y2S2*z2S + eta_t*y2S2*zqd + eta_t*y2S*yqc*z2S + eta_t*y2S*yqc*zqd + y1D*y2S*z2S + y1D*y2S*zqd + y1D*yqc*z2S + y1D*yqc*zqd + y2S*ypb*z2S + y2S*ypb*zqd + ypb*yqc*z2S + ypb*yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1230) ! | px  py  pz  s    ( 109) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypb*z2S + ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1231) ! | px  py  pz  px   ( 110) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypb*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypb*z2S + ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1232) ! | px  py  pz  py   ( 111) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypb*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqc - eta_p*eta_q*ypq*yqd*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqc + eta_p*eta_t*y2S*yqd*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqc + eta_p*y1D*yqd*zpq + eta_p*y2S*ypb*zpq + eta_p*ypb*ypq*z2S + eta_p*ypb*ypq*zqc + eta_p*ypb*yqd*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqc - eta_q*ypq*yqd*z2S - eta_q*ypq*yqd*zqc + eta_t*y2S2*z2S + eta_t*y2S2*zqc + eta_t*y2S*yqd*z2S + eta_t*y2S*yqd*zqc + y1D*y2S*z2S + y1D*y2S*zqc + y1D*yqd*z2S + y1D*yqd*zqc + y2S*ypb*z2S + y2S*ypb*zqc + ypb*yqd*z2S + ypb*yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1233) ! | px  py  pz  pz   ( 112) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*y2S*zpq2 + eta_p2*y1D*zpq2 + eta_p2*ypb*zpq2 - 2*eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*ypq*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*y2S*zpq*zqc + eta_p*eta_t*y2S*zpq*zqd + 2*eta_p*y1D*z2S*zpq + eta_p*y1D*zpq*zqc + eta_p*y1D*zpq*zqd + 2*eta_p*ypb*z2S*zpq + eta_p*ypb*zpq*zqc + eta_p*ypb*zpq*zqd - eta_q*ypq*z2S2 - eta_q*ypq*z2S*zqc - eta_q*ypq*z2S*zqd - eta_q*ypq*zqc*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*y2S*z2S*zqd + eta_t*y2S*zqc*zqd + y1D*z2S2 + y1D*z2S*zqc + y1D*z2S*zqd + y1D*zqc*zqd + ypb*z2S2 + ypb*z2S*zqc + ypb*z2S*zqd + ypb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1300) ! | px  pz  s   s    ( 113) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1301) ! | px  pz  s   px   ( 114) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1302) ! | px  pz  s   py   ( 115) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpb + yqd*z1D + yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1303) ! | px  pz  s   pz   ( 116) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpb + zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1310) ! | px  pz  px  s    ( 117) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1311) ! | px  pz  px  px   ( 118) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpb + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1312) ! | px  pz  px  py   ( 119) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpb + yqd*z1D + yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1313) ! | px  pz  px  pz   ( 120) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpb + zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1320) ! | px  pz  py  s    ( 121) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpb + yqc*z1D + yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1321) ! | px  pz  py  px   ( 122) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpb - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpb + yqc*z1D + yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1322) ! | px  pz  py  py   ( 123) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*ypq2*z2S + eta_p2*ypq2*z1D + eta_p2*ypq2*zpb - 2*eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq*yqc*zpq - eta_p*eta_q*ypq*yqd*zpq + 2*eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*ypq*yqc*z2S + eta_p*eta_t*ypq*yqd*z2S + 2*eta_p*y2S*ypq*z1D + 2*eta_p*y2S*ypq*zpb + eta_p*ypq*yqc*z1D + eta_p*ypq*yqc*zpb + eta_p*ypq*yqd*z1D + eta_p*ypq*yqd*zpb - eta_q*y2S2*zpq - eta_q*y2S*yqc*zpq - eta_q*y2S*yqd*zpq - eta_q*yqc*yqd*zpq + eta_t*y2S2*z2S + eta_t*y2S*yqc*z2S + eta_t*y2S*yqd*z2S + eta_t*yqc*yqd*z2S + y2S2*z1D + y2S2*zpb + y2S*yqc*z1D + y2S*yqc*zpb + y2S*yqd*z1D + y2S*yqd*zpb + yqc*yqd*z1D + yqc*yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1323) ! | px  pz  py  pz   ( 124) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpb*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqd - eta_p*eta_q*yqc*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqd + eta_p*eta_t*yqc*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpb*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqd + eta_p*ypq*z2S*zpb + eta_p*ypq*zpb*zqd + eta_p*yqc*z1D*zpq + eta_p*yqc*zpb*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqd - eta_q*yqc*z2S*zpq - eta_q*yqc*zpq*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqd + eta_t*yqc*z2S2 + eta_t*yqc*z2S*zqd + y2S*z1D*z2S + y2S*z1D*zqd + y2S*z2S*zpb + y2S*zpb*zqd + yqc*z1D*z2S + yqc*z1D*zqd + yqc*z2S*zpb + yqc*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1330) ! | px  pz  pz  s    ( 125) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpb + zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1331) ! | px  pz  pz  px   ( 126) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpb*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpb + zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1332) ! | px  pz  pz  py   ( 127) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpb*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*yqd*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqc + eta_p*eta_t*yqd*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpb*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqc + eta_p*ypq*z2S*zpb + eta_p*ypq*zpb*zqc + eta_p*yqd*z1D*zpq + eta_p*yqd*zpb*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqc - eta_q*yqd*z2S*zpq - eta_q*yqd*zpq*zqc + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*yqd*z2S2 + eta_t*yqd*z2S*zqc + y2S*z1D*z2S + y2S*z1D*zqc + y2S*z2S*zpb + y2S*zpb*zqc + yqd*z1D*z2S + yqd*z1D*zqc + yqd*z2S*zpb + yqd*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (1333) ! | px  pz  pz  pz   ( 128) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*zpq3 + eta_p2*eta_t*z2S*zpq2 + eta_p2*z1D*zpq2 + eta_p2*zpb*zpq2 - 2*eta_p*eta_q*z2S*zpq2 - eta_p*eta_q*zpq2*zqc - eta_p*eta_q*zpq2*zqd + 2*eta_p*eta_t*z2S2*zpq + eta_p*eta_t*z2S*zpq*zqc + eta_p*eta_t*z2S*zpq*zqd + 2*eta_p*z1D*z2S*zpq + eta_p*z1D*zpq*zqc + eta_p*z1D*zpq*zqd + 2*eta_p*z2S*zpb*zpq + eta_p*zpb*zpq*zqc + eta_p*zpb*zpq*zqd - eta_q*z2S2*zpq - eta_q*z2S*zpq*zqc - eta_q*z2S*zpq*zqd - eta_q*zpq*zqc*zqd + eta_t*z2S3 + eta_t*z2S2*zqc + eta_t*z2S2*zqd + eta_t*z2S*zqc*zqd + z1D*z2S2 + z1D*z2S*zqc + z1D*z2S*zqd + z1D*zqc*zqd + z2S2*zpb + z2S*zpb*zqc + z2S*zpb*zqd + zpb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpa * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpa * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpa * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2000) ! | py  s   s   s    ( 129) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2001) ! | py  s   s   px   ( 130) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2002) ! | py  s   s   py   ( 131) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypa + ypa*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2003) ! | py  s   s   pz   ( 132) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypa*z2S + ypa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2010) ! | py  s   px  s    ( 133) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2011) ! | py  s   px  px   ( 134) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2012) ! | py  s   px  py   ( 135) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypa + ypa*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2013) ! | py  s   px  pz   ( 136) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypa*z2S + ypa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2020) ! | py  s   py  s    ( 137) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypa + ypa*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2021) ! | py  s   py  px   ( 138) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypa + ypa*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2022) ! | py  s   py  py   ( 139) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq3 + eta_p2*eta_t*y2S*ypq2 + eta_p2*y1D*ypq2 + eta_p2*ypa*ypq2 - 2*eta_p*eta_q*y2S*ypq2 - eta_p*eta_q*ypq2*yqc - eta_p*eta_q*ypq2*yqd + 2*eta_p*eta_t*y2S2*ypq + eta_p*eta_t*y2S*ypq*yqc + eta_p*eta_t*y2S*ypq*yqd + 2*eta_p*y1D*y2S*ypq + eta_p*y1D*ypq*yqc + eta_p*y1D*ypq*yqd + 2*eta_p*y2S*ypa*ypq + eta_p*ypa*ypq*yqc + eta_p*ypa*ypq*yqd - eta_q*y2S2*ypq - eta_q*y2S*ypq*yqc - eta_q*y2S*ypq*yqd - eta_q*ypq*yqc*yqd + eta_t*y2S3 + eta_t*y2S2*yqc + eta_t*y2S2*yqd + eta_t*y2S*yqc*yqd + y1D*y2S2 + y1D*y2S*yqc + y1D*y2S*yqd + y1D*yqc*yqd + y2S2*ypa + y2S*ypa*yqc + y2S*ypa*yqd + ypa*yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2023) ! | py  s   py  pz   ( 140) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypa*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqd - eta_p*eta_q*ypq*yqc*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqd + eta_p*eta_t*y2S*yqc*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqd + eta_p*y1D*yqc*zpq + eta_p*y2S*ypa*zpq + eta_p*ypa*ypq*z2S + eta_p*ypa*ypq*zqd + eta_p*ypa*yqc*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqd - eta_q*ypq*yqc*z2S - eta_q*ypq*yqc*zqd + eta_t*y2S2*z2S + eta_t*y2S2*zqd + eta_t*y2S*yqc*z2S + eta_t*y2S*yqc*zqd + y1D*y2S*z2S + y1D*y2S*zqd + y1D*yqc*z2S + y1D*yqc*zqd + y2S*ypa*z2S + y2S*ypa*zqd + ypa*yqc*z2S + ypa*yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2030) ! | py  s   pz  s    ( 141) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypa*z2S + ypa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2031) ! | py  s   pz  px   ( 142) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypa*z2S + ypa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2032) ! | py  s   pz  py   ( 143) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypa*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqc - eta_p*eta_q*ypq*yqd*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqc + eta_p*eta_t*y2S*yqd*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqc + eta_p*y1D*yqd*zpq + eta_p*y2S*ypa*zpq + eta_p*ypa*ypq*z2S + eta_p*ypa*ypq*zqc + eta_p*ypa*yqd*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqc - eta_q*ypq*yqd*z2S - eta_q*ypq*yqd*zqc + eta_t*y2S2*z2S + eta_t*y2S2*zqc + eta_t*y2S*yqd*z2S + eta_t*y2S*yqd*zqc + y1D*y2S*z2S + y1D*y2S*zqc + y1D*yqd*z2S + y1D*yqd*zqc + y2S*ypa*z2S + y2S*ypa*zqc + ypa*yqd*z2S + ypa*yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2033) ! | py  s   pz  pz   ( 144) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*y2S*zpq2 + eta_p2*y1D*zpq2 + eta_p2*ypa*zpq2 - 2*eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*ypq*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*y2S*zpq*zqc + eta_p*eta_t*y2S*zpq*zqd + 2*eta_p*y1D*z2S*zpq + eta_p*y1D*zpq*zqc + eta_p*y1D*zpq*zqd + 2*eta_p*ypa*z2S*zpq + eta_p*ypa*zpq*zqc + eta_p*ypa*zpq*zqd - eta_q*ypq*z2S2 - eta_q*ypq*z2S*zqc - eta_q*ypq*z2S*zqd - eta_q*ypq*zqc*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*y2S*z2S*zqd + eta_t*y2S*zqc*zqd + y1D*z2S2 + y1D*z2S*zqc + y1D*z2S*zqd + y1D*zqc*zqd + ypa*z2S2 + ypa*z2S*zqc + ypa*z2S*zqd + ypa*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2100) ! | py  px  s   s    ( 145) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2101) ! | py  px  s   px   ( 146) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2102) ! | py  px  s   py   ( 147) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypa + ypa*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2103) ! | py  px  s   pz   ( 148) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypa*z2S + ypa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2110) ! | py  px  px  s    ( 149) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2111) ! | py  px  px  px   ( 150) 
      
      n         = 0
      const     = pi2 * D2 * (y1D + ypa + eta_t * y2S - eta_q * ypq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2112) ! | py  px  px  py   ( 151) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqd + eta_t*y2S2 + eta_t*y2S*yqd + y1D*y2S + y1D*yqd + y2S*ypa + ypa*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2113) ! | py  px  px  pz   ( 152) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqd + eta_t*y2S*z2S + eta_t*y2S*zqd + y1D*z2S + y1D*zqd + ypa*z2S + ypa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2120) ! | py  px  py  s    ( 153) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypa + ypa*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2121) ! | py  px  py  px   ( 154) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq2 + eta_p*eta_t*y2S*ypq + eta_p*y1D*ypq + eta_p*ypa*ypq - eta_q*y2S*ypq - eta_q*ypq*yqc + eta_t*y2S2 + eta_t*y2S*yqc + y1D*y2S + y1D*yqc + y2S*ypa + ypa*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2122) ! | py  px  py  py   ( 155) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq3 + eta_p2*eta_t*y2S*ypq2 + eta_p2*y1D*ypq2 + eta_p2*ypa*ypq2 - 2*eta_p*eta_q*y2S*ypq2 - eta_p*eta_q*ypq2*yqc - eta_p*eta_q*ypq2*yqd + 2*eta_p*eta_t*y2S2*ypq + eta_p*eta_t*y2S*ypq*yqc + eta_p*eta_t*y2S*ypq*yqd + 2*eta_p*y1D*y2S*ypq + eta_p*y1D*ypq*yqc + eta_p*y1D*ypq*yqd + 2*eta_p*y2S*ypa*ypq + eta_p*ypa*ypq*yqc + eta_p*ypa*ypq*yqd - eta_q*y2S2*ypq - eta_q*y2S*ypq*yqc - eta_q*y2S*ypq*yqd - eta_q*ypq*yqc*yqd + eta_t*y2S3 + eta_t*y2S2*yqc + eta_t*y2S2*yqd + eta_t*y2S*yqc*yqd + y1D*y2S2 + y1D*y2S*yqc + y1D*y2S*yqd + y1D*yqc*yqd + y2S2*ypa + y2S*ypa*yqc + y2S*ypa*yqd + ypa*yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2123) ! | py  px  py  pz   ( 156) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypa*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqd - eta_p*eta_q*ypq*yqc*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqd + eta_p*eta_t*y2S*yqc*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqd + eta_p*y1D*yqc*zpq + eta_p*y2S*ypa*zpq + eta_p*ypa*ypq*z2S + eta_p*ypa*ypq*zqd + eta_p*ypa*yqc*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqd - eta_q*ypq*yqc*z2S - eta_q*ypq*yqc*zqd + eta_t*y2S2*z2S + eta_t*y2S2*zqd + eta_t*y2S*yqc*z2S + eta_t*y2S*yqc*zqd + y1D*y2S*z2S + y1D*y2S*zqd + y1D*yqc*z2S + y1D*yqc*zqd + y2S*ypa*z2S + y2S*ypa*zqd + ypa*yqc*z2S + ypa*yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2130) ! | py  px  pz  s    ( 157) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypa*z2S + ypa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2131) ! | py  px  pz  px   ( 158) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*y2S*zpq + eta_p*y1D*zpq + eta_p*ypa*zpq - eta_q*ypq*z2S - eta_q*ypq*zqc + eta_t*y2S*z2S + eta_t*y2S*zqc + y1D*z2S + y1D*zqc + ypa*z2S + ypa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2132) ! | py  px  pz  py   ( 159) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*y2S*ypq*zpq + eta_p2*y1D*ypq*zpq + eta_p2*ypa*ypq*zpq - eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq2*z2S - eta_p*eta_q*ypq2*zqc - eta_p*eta_q*ypq*yqd*zpq + eta_p*eta_t*y2S2*zpq + eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*y2S*ypq*zqc + eta_p*eta_t*y2S*yqd*zpq + eta_p*y1D*y2S*zpq + eta_p*y1D*ypq*z2S + eta_p*y1D*ypq*zqc + eta_p*y1D*yqd*zpq + eta_p*y2S*ypa*zpq + eta_p*ypa*ypq*z2S + eta_p*ypa*ypq*zqc + eta_p*ypa*yqd*zpq - eta_q*y2S*ypq*z2S - eta_q*y2S*ypq*zqc - eta_q*ypq*yqd*z2S - eta_q*ypq*yqd*zqc + eta_t*y2S2*z2S + eta_t*y2S2*zqc + eta_t*y2S*yqd*z2S + eta_t*y2S*yqd*zqc + y1D*y2S*z2S + y1D*y2S*zqc + y1D*yqd*z2S + y1D*yqd*zqc + y2S*ypa*z2S + y2S*ypa*zqc + ypa*yqd*z2S + ypa*yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2133) ! | py  px  pz  pz   ( 160) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*y2S*zpq2 + eta_p2*y1D*zpq2 + eta_p2*ypa*zpq2 - 2*eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*ypq*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*y2S*zpq*zqc + eta_p*eta_t*y2S*zpq*zqd + 2*eta_p*y1D*z2S*zpq + eta_p*y1D*zpq*zqc + eta_p*y1D*zpq*zqd + 2*eta_p*ypa*z2S*zpq + eta_p*ypa*zpq*zqc + eta_p*ypa*zpq*zqd - eta_q*ypq*z2S2 - eta_q*ypq*z2S*zqc - eta_q*ypq*z2S*zqd - eta_q*ypq*zqc*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*y2S*z2S*zqd + eta_t*y2S*zqc*zqd + y1D*z2S2 + y1D*z2S*zqc + y1D*z2S*zqd + y1D*zqc*zqd + ypa*z2S2 + ypa*z2S*zqc + ypa*z2S*zqd + ypa*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2200) ! | py  py  s   s    ( 161) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq2 - 2*eta_q*eta_t*y2S*ypq - 2*eta_q*y1D*ypq - eta_q*ypa*ypq - eta_q*ypb*ypq + eta_t2*y2S2 + 2*eta_t*y1D*y2S + eta_t*y2S*ypa + eta_t*y2S*ypb + y1D2 + y1D*ypa + y1D*ypb + ypa*ypb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2201) ! | py  py  s   px   ( 162) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq2 - 2*eta_q*eta_t*y2S*ypq - 2*eta_q*y1D*ypq - eta_q*ypa*ypq - eta_q*ypb*ypq + eta_t2*y2S2 + 2*eta_t*y1D*y2S + eta_t*y2S*ypa + eta_t*y2S*ypb + y1D2 + y1D*ypa + y1D*ypb + ypa*ypb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2202) ! | py  py  s   py   ( 163) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq3 - 2*eta_p*eta_q*eta_t*y2S*ypq2 - 2*eta_p*eta_q*y1D*ypq2 - eta_p*eta_q*ypa*ypq2 - eta_p*eta_q*ypb*ypq2 + eta_p*eta_t2*y2S2*ypq + 2*eta_p*eta_t*y1D*y2S*ypq + eta_p*eta_t*y2S*ypa*ypq + eta_p*eta_t*y2S*ypb*ypq + eta_p*y1D2*ypq + eta_p*y1D*ypa*ypq + eta_p*y1D*ypb*ypq + eta_p*ypa*ypb*ypq + eta_q2*y2S*ypq2 + eta_q2*ypq2*yqd - 2*eta_q*eta_t*y2S2*ypq - 2*eta_q*eta_t*y2S*ypq*yqd - 2*eta_q*y1D*y2S*ypq - 2*eta_q*y1D*ypq*yqd - eta_q*y2S*ypa*ypq - eta_q*y2S*ypb*ypq - eta_q*ypa*ypq*yqd - eta_q*ypb*ypq*yqd + eta_t2*y2S3 + eta_t2*y2S2*yqd + 2*eta_t*y1D*y2S2 + 2*eta_t*y1D*y2S*yqd + eta_t*y2S2*ypa + eta_t*y2S2*ypb + eta_t*y2S*ypa*yqd + eta_t*y2S*ypb*yqd + y1D2*y2S + y1D2*yqd + y1D*y2S*ypa + y1D*y2S*ypb + y1D*ypa*yqd + y1D*ypb*yqd + y2S*ypa*ypb + ypa*ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2203) ! | py  py  s   pz   ( 164) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq + eta_p*eta_t2*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*zpq + eta_p*eta_t*y2S*ypa*zpq + eta_p*eta_t*y2S*ypb*zpq + eta_p*y1D2*zpq + eta_p*y1D*ypa*zpq + eta_p*y1D*ypb*zpq + eta_p*ypa*ypb*zpq + eta_q2*ypq2*z2S + eta_q2*ypq2*zqd - 2*eta_q*eta_t*y2S*ypq*z2S - 2*eta_q*eta_t*y2S*ypq*zqd - 2*eta_q*y1D*ypq*z2S - 2*eta_q*y1D*ypq*zqd - eta_q*ypa*ypq*z2S - eta_q*ypa*ypq*zqd - eta_q*ypb*ypq*z2S - eta_q*ypb*ypq*zqd + eta_t2*y2S2*z2S + eta_t2*y2S2*zqd + 2*eta_t*y1D*y2S*z2S + 2*eta_t*y1D*y2S*zqd + eta_t*y2S*ypa*z2S + eta_t*y2S*ypa*zqd + eta_t*y2S*ypb*z2S + eta_t*y2S*ypb*zqd + y1D2*z2S + y1D2*zqd + y1D*ypa*z2S + y1D*ypa*zqd + y1D*ypb*z2S + y1D*ypb*zqd + ypa*ypb*z2S + ypa*ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2210) ! | py  py  px  s    ( 165) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq2 - 2*eta_q*eta_t*y2S*ypq - 2*eta_q*y1D*ypq - eta_q*ypa*ypq - eta_q*ypb*ypq + eta_t2*y2S2 + 2*eta_t*y1D*y2S + eta_t*y2S*ypa + eta_t*y2S*ypb + y1D2 + y1D*ypa + y1D*ypb + ypa*ypb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2211) ! | py  py  px  px   ( 166) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq2 - 2*eta_q*eta_t*y2S*ypq - 2*eta_q*y1D*ypq - eta_q*ypa*ypq - eta_q*ypb*ypq + eta_t2*y2S2 + 2*eta_t*y1D*y2S + eta_t*y2S*ypa + eta_t*y2S*ypb + y1D2 + y1D*ypa + y1D*ypb + ypa*ypb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2212) ! | py  py  px  py   ( 167) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq3 - 2*eta_p*eta_q*eta_t*y2S*ypq2 - 2*eta_p*eta_q*y1D*ypq2 - eta_p*eta_q*ypa*ypq2 - eta_p*eta_q*ypb*ypq2 + eta_p*eta_t2*y2S2*ypq + 2*eta_p*eta_t*y1D*y2S*ypq + eta_p*eta_t*y2S*ypa*ypq + eta_p*eta_t*y2S*ypb*ypq + eta_p*y1D2*ypq + eta_p*y1D*ypa*ypq + eta_p*y1D*ypb*ypq + eta_p*ypa*ypb*ypq + eta_q2*y2S*ypq2 + eta_q2*ypq2*yqd - 2*eta_q*eta_t*y2S2*ypq - 2*eta_q*eta_t*y2S*ypq*yqd - 2*eta_q*y1D*y2S*ypq - 2*eta_q*y1D*ypq*yqd - eta_q*y2S*ypa*ypq - eta_q*y2S*ypb*ypq - eta_q*ypa*ypq*yqd - eta_q*ypb*ypq*yqd + eta_t2*y2S3 + eta_t2*y2S2*yqd + 2*eta_t*y1D*y2S2 + 2*eta_t*y1D*y2S*yqd + eta_t*y2S2*ypa + eta_t*y2S2*ypb + eta_t*y2S*ypa*yqd + eta_t*y2S*ypb*yqd + y1D2*y2S + y1D2*yqd + y1D*y2S*ypa + y1D*y2S*ypb + y1D*ypa*yqd + y1D*ypb*yqd + y2S*ypa*ypb + ypa*ypb*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2213) ! | py  py  px  pz   ( 168) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq + eta_p*eta_t2*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*zpq + eta_p*eta_t*y2S*ypa*zpq + eta_p*eta_t*y2S*ypb*zpq + eta_p*y1D2*zpq + eta_p*y1D*ypa*zpq + eta_p*y1D*ypb*zpq + eta_p*ypa*ypb*zpq + eta_q2*ypq2*z2S + eta_q2*ypq2*zqd - 2*eta_q*eta_t*y2S*ypq*z2S - 2*eta_q*eta_t*y2S*ypq*zqd - 2*eta_q*y1D*ypq*z2S - 2*eta_q*y1D*ypq*zqd - eta_q*ypa*ypq*z2S - eta_q*ypa*ypq*zqd - eta_q*ypb*ypq*z2S - eta_q*ypb*ypq*zqd + eta_t2*y2S2*z2S + eta_t2*y2S2*zqd + 2*eta_t*y1D*y2S*z2S + 2*eta_t*y1D*y2S*zqd + eta_t*y2S*ypa*z2S + eta_t*y2S*ypa*zqd + eta_t*y2S*ypb*z2S + eta_t*y2S*ypb*zqd + y1D2*z2S + y1D2*zqd + y1D*ypa*z2S + y1D*ypa*zqd + y1D*ypb*z2S + y1D*ypb*zqd + ypa*ypb*z2S + ypa*ypb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2220) ! | py  py  py  s    ( 169) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq3 - 2*eta_p*eta_q*eta_t*y2S*ypq2 - 2*eta_p*eta_q*y1D*ypq2 - eta_p*eta_q*ypa*ypq2 - eta_p*eta_q*ypb*ypq2 + eta_p*eta_t2*y2S2*ypq + 2*eta_p*eta_t*y1D*y2S*ypq + eta_p*eta_t*y2S*ypa*ypq + eta_p*eta_t*y2S*ypb*ypq + eta_p*y1D2*ypq + eta_p*y1D*ypa*ypq + eta_p*y1D*ypb*ypq + eta_p*ypa*ypb*ypq + eta_q2*y2S*ypq2 + eta_q2*ypq2*yqc - 2*eta_q*eta_t*y2S2*ypq - 2*eta_q*eta_t*y2S*ypq*yqc - 2*eta_q*y1D*y2S*ypq - 2*eta_q*y1D*ypq*yqc - eta_q*y2S*ypa*ypq - eta_q*y2S*ypb*ypq - eta_q*ypa*ypq*yqc - eta_q*ypb*ypq*yqc + eta_t2*y2S3 + eta_t2*y2S2*yqc + 2*eta_t*y1D*y2S2 + 2*eta_t*y1D*y2S*yqc + eta_t*y2S2*ypa + eta_t*y2S2*ypb + eta_t*y2S*ypa*yqc + eta_t*y2S*ypb*yqc + y1D2*y2S + y1D2*yqc + y1D*y2S*ypa + y1D*y2S*ypb + y1D*ypa*yqc + y1D*ypb*yqc + y2S*ypa*ypb + ypa*ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2221) ! | py  py  py  px   ( 170) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq3 - 2*eta_p*eta_q*eta_t*y2S*ypq2 - 2*eta_p*eta_q*y1D*ypq2 - eta_p*eta_q*ypa*ypq2 - eta_p*eta_q*ypb*ypq2 + eta_p*eta_t2*y2S2*ypq + 2*eta_p*eta_t*y1D*y2S*ypq + eta_p*eta_t*y2S*ypa*ypq + eta_p*eta_t*y2S*ypb*ypq + eta_p*y1D2*ypq + eta_p*y1D*ypa*ypq + eta_p*y1D*ypb*ypq + eta_p*ypa*ypb*ypq + eta_q2*y2S*ypq2 + eta_q2*ypq2*yqc - 2*eta_q*eta_t*y2S2*ypq - 2*eta_q*eta_t*y2S*ypq*yqc - 2*eta_q*y1D*y2S*ypq - 2*eta_q*y1D*ypq*yqc - eta_q*y2S*ypa*ypq - eta_q*y2S*ypb*ypq - eta_q*ypa*ypq*yqc - eta_q*ypb*ypq*yqc + eta_t2*y2S3 + eta_t2*y2S2*yqc + 2*eta_t*y1D*y2S2 + 2*eta_t*y1D*y2S*yqc + eta_t*y2S2*ypa + eta_t*y2S2*ypb + eta_t*y2S*ypa*yqc + eta_t*y2S*ypb*yqc + y1D2*y2S + y1D2*yqc + y1D*y2S*ypa + y1D*y2S*ypb + y1D*ypa*yqc + y1D*ypb*yqc + y2S*ypa*ypb + ypa*ypb*yqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2222) ! | py  py  py  py   ( 171) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq4 - 2*eta_p2*eta_q*eta_t*y2S*ypq3 - 2*eta_p2*eta_q*y1D*ypq3 - eta_p2*eta_q*ypa*ypq3 - eta_p2*eta_q*ypb*ypq3 + eta_p2*eta_t2*y2S2*ypq2 + 2*eta_p2*eta_t*y1D*y2S*ypq2 + eta_p2*eta_t*y2S*ypa*ypq2 + eta_p2*eta_t*y2S*ypb*ypq2 + eta_p2*y1D2*ypq2 + eta_p2*y1D*ypa*ypq2 + eta_p2*y1D*ypb*ypq2 + eta_p2*ypa*ypb*ypq2 + 2*eta_p*eta_q2*y2S*ypq3 + eta_p*eta_q2*ypq3*yqc + eta_p*eta_q2*ypq3*yqd - 4*eta_p*eta_q*eta_t*y2S2*ypq2 - 2*eta_p*eta_q*eta_t*y2S*ypq2*yqc - 2*eta_p*eta_q*eta_t*y2S*ypq2*yqd - 4*eta_p*eta_q*y1D*y2S*ypq2 - 2*eta_p*eta_q*y1D*ypq2*yqc - 2*eta_p*eta_q*y1D*ypq2*yqd - 2*eta_p*eta_q*y2S*ypa*ypq2 - 2*eta_p*eta_q*y2S*ypb*ypq2 - eta_p*eta_q*ypa*ypq2*yqc - eta_p*eta_q*ypa*ypq2*yqd - eta_p*eta_q*ypb*ypq2*yqc - eta_p*eta_q*ypb*ypq2*yqd + 2*eta_p*eta_t2*y2S3*ypq + eta_p*eta_t2*y2S2*ypq*yqc + eta_p*eta_t2*y2S2*ypq*yqd + 4*eta_p*eta_t*y1D*y2S2*ypq + 2*eta_p*eta_t*y1D*y2S*ypq*yqc + 2*eta_p*eta_t*y1D*y2S*ypq*yqd + 2*eta_p*eta_t*y2S2*ypa*ypq + 2*eta_p*eta_t*y2S2*ypb*ypq + eta_p*eta_t*y2S*ypa*ypq*yqc + eta_p*eta_t*y2S*ypa*ypq*yqd + eta_p*eta_t*y2S*ypb*ypq*yqc + eta_p*eta_t*y2S*ypb*ypq*yqd + 2*eta_p*y1D2*y2S*ypq + eta_p*y1D2*ypq*yqc + eta_p*y1D2*ypq*yqd + 2*eta_p*y1D*y2S*ypa*ypq + 2*eta_p*y1D*y2S*ypb*ypq + eta_p*y1D*ypa*ypq*yqc + eta_p*y1D*ypa*ypq*yqd + eta_p*y1D*ypb*ypq*yqc + eta_p*y1D*ypb*ypq*yqd + 2*eta_p*y2S*ypa*ypb*ypq + eta_p*ypa*ypb*ypq*yqc + eta_p*ypa*ypb*ypq*yqd + eta_q2*y2S2*ypq2 + eta_q2*y2S*ypq2*yqc + eta_q2*y2S*ypq2*yqd + eta_q2*ypq2*yqc*yqd - 2*eta_q*eta_t*y2S3*ypq - 2*eta_q*eta_t*y2S2*ypq*yqc - 2*eta_q*eta_t*y2S2*ypq*yqd - 2*eta_q*eta_t*y2S*ypq*yqc*yqd - 2*eta_q*y1D*y2S2*ypq - 2*eta_q*y1D*y2S*ypq*yqc - 2*eta_q*y1D*y2S*ypq*yqd - 2*eta_q*y1D*ypq*yqc*yqd - eta_q*y2S2*ypa*ypq - eta_q*y2S2*ypb*ypq - eta_q*y2S*ypa*ypq*yqc - eta_q*y2S*ypa*ypq*yqd - eta_q*y2S*ypb*ypq*yqc - eta_q*y2S*ypb*ypq*yqd - eta_q*ypa*ypq*yqc*yqd - eta_q*ypb*ypq*yqc*yqd + eta_t2*y2S4 + eta_t2*y2S3*yqc + eta_t2*y2S3*yqd + eta_t2*y2S2*yqc*yqd + 2*eta_t*y1D*y2S3 + 2*eta_t*y1D*y2S2*yqc + 2*eta_t*y1D*y2S2*yqd + 2*eta_t*y1D*y2S*yqc*yqd + eta_t*y2S3*ypa + eta_t*y2S3*ypb + eta_t*y2S2*ypa*yqc + eta_t*y2S2*ypa*yqd + eta_t*y2S2*ypb*yqc + eta_t*y2S2*ypb*yqd + eta_t*y2S*ypa*yqc*yqd + eta_t*y2S*ypb*yqc*yqd + y1D2*y2S2 + y1D2*y2S*yqc + y1D2*y2S*yqd + y1D2*yqc*yqd + y1D*y2S2*ypa + y1D*y2S2*ypb + y1D*y2S*ypa*yqc + y1D*y2S*ypa*yqd + y1D*y2S*ypb*yqc + y1D*y2S*ypb*yqd + y1D*ypa*yqc*yqd + y1D*ypb*yqc*yqd + y2S2*ypa*ypb + y2S*ypa*ypb*yqc + y2S*ypa*ypb*yqd + ypa*ypb*yqc*yqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2223) ! | py  py  py  pz   ( 172) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq3*zpq - 2*eta_p2*eta_q*eta_t*y2S*ypq2*zpq - 2*eta_p2*eta_q*y1D*ypq2*zpq - eta_p2*eta_q*ypa*ypq2*zpq - eta_p2*eta_q*ypb*ypq2*zpq + eta_p2*eta_t2*y2S2*ypq*zpq + 2*eta_p2*eta_t*y1D*y2S*ypq*zpq + eta_p2*eta_t*y2S*ypa*ypq*zpq + eta_p2*eta_t*y2S*ypb*ypq*zpq + eta_p2*y1D2*ypq*zpq + eta_p2*y1D*ypa*ypq*zpq + eta_p2*y1D*ypb*ypq*zpq + eta_p2*ypa*ypb*ypq*zpq + eta_p*eta_q2*y2S*ypq2*zpq + eta_p*eta_q2*ypq3*z2S + eta_p*eta_q2*ypq3*zqd + eta_p*eta_q2*ypq2*yqc*zpq - 2*eta_p*eta_q*eta_t*y2S2*ypq*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq2*z2S - 2*eta_p*eta_q*eta_t*y2S*ypq2*zqd - 2*eta_p*eta_q*eta_t*y2S*ypq*yqc*zpq - 2*eta_p*eta_q*y1D*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq2*z2S - 2*eta_p*eta_q*y1D*ypq2*zqd - 2*eta_p*eta_q*y1D*ypq*yqc*zpq - eta_p*eta_q*y2S*ypa*ypq*zpq - eta_p*eta_q*y2S*ypb*ypq*zpq - eta_p*eta_q*ypa*ypq2*z2S - eta_p*eta_q*ypa*ypq2*zqd - eta_p*eta_q*ypa*ypq*yqc*zpq - eta_p*eta_q*ypb*ypq2*z2S - eta_p*eta_q*ypb*ypq2*zqd - eta_p*eta_q*ypb*ypq*yqc*zpq + eta_p*eta_t2*y2S3*zpq + eta_p*eta_t2*y2S2*ypq*z2S + eta_p*eta_t2*y2S2*ypq*zqd + eta_p*eta_t2*y2S2*yqc*zpq + 2*eta_p*eta_t*y1D*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*ypq*z2S + 2*eta_p*eta_t*y1D*y2S*ypq*zqd + 2*eta_p*eta_t*y1D*y2S*yqc*zpq + eta_p*eta_t*y2S2*ypa*zpq + eta_p*eta_t*y2S2*ypb*zpq + eta_p*eta_t*y2S*ypa*ypq*z2S + eta_p*eta_t*y2S*ypa*ypq*zqd + eta_p*eta_t*y2S*ypa*yqc*zpq + eta_p*eta_t*y2S*ypb*ypq*z2S + eta_p*eta_t*y2S*ypb*ypq*zqd + eta_p*eta_t*y2S*ypb*yqc*zpq + eta_p*y1D2*y2S*zpq + eta_p*y1D2*ypq*z2S + eta_p*y1D2*ypq*zqd + eta_p*y1D2*yqc*zpq + eta_p*y1D*y2S*ypa*zpq + eta_p*y1D*y2S*ypb*zpq + eta_p*y1D*ypa*ypq*z2S + eta_p*y1D*ypa*ypq*zqd + eta_p*y1D*ypa*yqc*zpq + eta_p*y1D*ypb*ypq*z2S + eta_p*y1D*ypb*ypq*zqd + eta_p*y1D*ypb*yqc*zpq + eta_p*y2S*ypa*ypb*zpq + eta_p*ypa*ypb*ypq*z2S + eta_p*ypa*ypb*ypq*zqd + eta_p*ypa*ypb*yqc*zpq + eta_q2*y2S*ypq2*z2S + eta_q2*y2S*ypq2*zqd + eta_q2*ypq2*yqc*z2S + eta_q2*ypq2*yqc*zqd - 2*eta_q*eta_t*y2S2*ypq*z2S - 2*eta_q*eta_t*y2S2*ypq*zqd - 2*eta_q*eta_t*y2S*ypq*yqc*z2S - 2*eta_q*eta_t*y2S*ypq*yqc*zqd - 2*eta_q*y1D*y2S*ypq*z2S - 2*eta_q*y1D*y2S*ypq*zqd - 2*eta_q*y1D*ypq*yqc*z2S - 2*eta_q*y1D*ypq*yqc*zqd - eta_q*y2S*ypa*ypq*z2S - eta_q*y2S*ypa*ypq*zqd - eta_q*y2S*ypb*ypq*z2S - eta_q*y2S*ypb*ypq*zqd - eta_q*ypa*ypq*yqc*z2S - eta_q*ypa*ypq*yqc*zqd - eta_q*ypb*ypq*yqc*z2S - eta_q*ypb*ypq*yqc*zqd + eta_t2*y2S3*z2S + eta_t2*y2S3*zqd + eta_t2*y2S2*yqc*z2S + eta_t2*y2S2*yqc*zqd + 2*eta_t*y1D*y2S2*z2S + 2*eta_t*y1D*y2S2*zqd + 2*eta_t*y1D*y2S*yqc*z2S + 2*eta_t*y1D*y2S*yqc*zqd + eta_t*y2S2*ypa*z2S + eta_t*y2S2*ypa*zqd + eta_t*y2S2*ypb*z2S + eta_t*y2S2*ypb*zqd + eta_t*y2S*ypa*yqc*z2S + eta_t*y2S*ypa*yqc*zqd + eta_t*y2S*ypb*yqc*z2S + eta_t*y2S*ypb*yqc*zqd + y1D2*y2S*z2S + y1D2*y2S*zqd + y1D2*yqc*z2S + y1D2*yqc*zqd + y1D*y2S*ypa*z2S + y1D*y2S*ypa*zqd + y1D*y2S*ypb*z2S + y1D*y2S*ypb*zqd + y1D*ypa*yqc*z2S + y1D*ypa*yqc*zqd + y1D*ypb*yqc*z2S + y1D*ypb*yqc*zqd + y2S*ypa*ypb*z2S + y2S*ypa*ypb*zqd + ypa*ypb*yqc*z2S + ypa*ypb*yqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2230) ! | py  py  pz  s    ( 173) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq + eta_p*eta_t2*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*zpq + eta_p*eta_t*y2S*ypa*zpq + eta_p*eta_t*y2S*ypb*zpq + eta_p*y1D2*zpq + eta_p*y1D*ypa*zpq + eta_p*y1D*ypb*zpq + eta_p*ypa*ypb*zpq + eta_q2*ypq2*z2S + eta_q2*ypq2*zqc - 2*eta_q*eta_t*y2S*ypq*z2S - 2*eta_q*eta_t*y2S*ypq*zqc - 2*eta_q*y1D*ypq*z2S - 2*eta_q*y1D*ypq*zqc - eta_q*ypa*ypq*z2S - eta_q*ypa*ypq*zqc - eta_q*ypb*ypq*z2S - eta_q*ypb*ypq*zqc + eta_t2*y2S2*z2S + eta_t2*y2S2*zqc + 2*eta_t*y1D*y2S*z2S + 2*eta_t*y1D*y2S*zqc + eta_t*y2S*ypa*z2S + eta_t*y2S*ypa*zqc + eta_t*y2S*ypb*z2S + eta_t*y2S*ypb*zqc + y1D2*z2S + y1D2*zqc + y1D*ypa*z2S + y1D*ypa*zqc + y1D*ypb*z2S + y1D*ypb*zqc + ypa*ypb*z2S + ypa*ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2231) ! | py  py  pz  px   ( 174) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq + eta_p*eta_t2*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*zpq + eta_p*eta_t*y2S*ypa*zpq + eta_p*eta_t*y2S*ypb*zpq + eta_p*y1D2*zpq + eta_p*y1D*ypa*zpq + eta_p*y1D*ypb*zpq + eta_p*ypa*ypb*zpq + eta_q2*ypq2*z2S + eta_q2*ypq2*zqc - 2*eta_q*eta_t*y2S*ypq*z2S - 2*eta_q*eta_t*y2S*ypq*zqc - 2*eta_q*y1D*ypq*z2S - 2*eta_q*y1D*ypq*zqc - eta_q*ypa*ypq*z2S - eta_q*ypa*ypq*zqc - eta_q*ypb*ypq*z2S - eta_q*ypb*ypq*zqc + eta_t2*y2S2*z2S + eta_t2*y2S2*zqc + 2*eta_t*y1D*y2S*z2S + 2*eta_t*y1D*y2S*zqc + eta_t*y2S*ypa*z2S + eta_t*y2S*ypa*zqc + eta_t*y2S*ypb*z2S + eta_t*y2S*ypb*zqc + y1D2*z2S + y1D2*zqc + y1D*ypa*z2S + y1D*ypa*zqc + y1D*ypb*z2S + y1D*ypb*zqc + ypa*ypb*z2S + ypa*ypb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2232) ! | py  py  pz  py   ( 175) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq3*zpq - 2*eta_p2*eta_q*eta_t*y2S*ypq2*zpq - 2*eta_p2*eta_q*y1D*ypq2*zpq - eta_p2*eta_q*ypa*ypq2*zpq - eta_p2*eta_q*ypb*ypq2*zpq + eta_p2*eta_t2*y2S2*ypq*zpq + 2*eta_p2*eta_t*y1D*y2S*ypq*zpq + eta_p2*eta_t*y2S*ypa*ypq*zpq + eta_p2*eta_t*y2S*ypb*ypq*zpq + eta_p2*y1D2*ypq*zpq + eta_p2*y1D*ypa*ypq*zpq + eta_p2*y1D*ypb*ypq*zpq + eta_p2*ypa*ypb*ypq*zpq + eta_p*eta_q2*y2S*ypq2*zpq + eta_p*eta_q2*ypq3*z2S + eta_p*eta_q2*ypq3*zqc + eta_p*eta_q2*ypq2*yqd*zpq - 2*eta_p*eta_q*eta_t*y2S2*ypq*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq2*z2S - 2*eta_p*eta_q*eta_t*y2S*ypq2*zqc - 2*eta_p*eta_q*eta_t*y2S*ypq*yqd*zpq - 2*eta_p*eta_q*y1D*y2S*ypq*zpq - 2*eta_p*eta_q*y1D*ypq2*z2S - 2*eta_p*eta_q*y1D*ypq2*zqc - 2*eta_p*eta_q*y1D*ypq*yqd*zpq - eta_p*eta_q*y2S*ypa*ypq*zpq - eta_p*eta_q*y2S*ypb*ypq*zpq - eta_p*eta_q*ypa*ypq2*z2S - eta_p*eta_q*ypa*ypq2*zqc - eta_p*eta_q*ypa*ypq*yqd*zpq - eta_p*eta_q*ypb*ypq2*z2S - eta_p*eta_q*ypb*ypq2*zqc - eta_p*eta_q*ypb*ypq*yqd*zpq + eta_p*eta_t2*y2S3*zpq + eta_p*eta_t2*y2S2*ypq*z2S + eta_p*eta_t2*y2S2*ypq*zqc + eta_p*eta_t2*y2S2*yqd*zpq + 2*eta_p*eta_t*y1D*y2S2*zpq + 2*eta_p*eta_t*y1D*y2S*ypq*z2S + 2*eta_p*eta_t*y1D*y2S*ypq*zqc + 2*eta_p*eta_t*y1D*y2S*yqd*zpq + eta_p*eta_t*y2S2*ypa*zpq + eta_p*eta_t*y2S2*ypb*zpq + eta_p*eta_t*y2S*ypa*ypq*z2S + eta_p*eta_t*y2S*ypa*ypq*zqc + eta_p*eta_t*y2S*ypa*yqd*zpq + eta_p*eta_t*y2S*ypb*ypq*z2S + eta_p*eta_t*y2S*ypb*ypq*zqc + eta_p*eta_t*y2S*ypb*yqd*zpq + eta_p*y1D2*y2S*zpq + eta_p*y1D2*ypq*z2S + eta_p*y1D2*ypq*zqc + eta_p*y1D2*yqd*zpq + eta_p*y1D*y2S*ypa*zpq + eta_p*y1D*y2S*ypb*zpq + eta_p*y1D*ypa*ypq*z2S + eta_p*y1D*ypa*ypq*zqc + eta_p*y1D*ypa*yqd*zpq + eta_p*y1D*ypb*ypq*z2S + eta_p*y1D*ypb*ypq*zqc + eta_p*y1D*ypb*yqd*zpq + eta_p*y2S*ypa*ypb*zpq + eta_p*ypa*ypb*ypq*z2S + eta_p*ypa*ypb*ypq*zqc + eta_p*ypa*ypb*yqd*zpq + eta_q2*y2S*ypq2*z2S + eta_q2*y2S*ypq2*zqc + eta_q2*ypq2*yqd*z2S + eta_q2*ypq2*yqd*zqc - 2*eta_q*eta_t*y2S2*ypq*z2S - 2*eta_q*eta_t*y2S2*ypq*zqc - 2*eta_q*eta_t*y2S*ypq*yqd*z2S - 2*eta_q*eta_t*y2S*ypq*yqd*zqc - 2*eta_q*y1D*y2S*ypq*z2S - 2*eta_q*y1D*y2S*ypq*zqc - 2*eta_q*y1D*ypq*yqd*z2S - 2*eta_q*y1D*ypq*yqd*zqc - eta_q*y2S*ypa*ypq*z2S - eta_q*y2S*ypa*ypq*zqc - eta_q*y2S*ypb*ypq*z2S - eta_q*y2S*ypb*ypq*zqc - eta_q*ypa*ypq*yqd*z2S - eta_q*ypa*ypq*yqd*zqc - eta_q*ypb*ypq*yqd*z2S - eta_q*ypb*ypq*yqd*zqc + eta_t2*y2S3*z2S + eta_t2*y2S3*zqc + eta_t2*y2S2*yqd*z2S + eta_t2*y2S2*yqd*zqc + 2*eta_t*y1D*y2S2*z2S + 2*eta_t*y1D*y2S2*zqc + 2*eta_t*y1D*y2S*yqd*z2S + 2*eta_t*y1D*y2S*yqd*zqc + eta_t*y2S2*ypa*z2S + eta_t*y2S2*ypa*zqc + eta_t*y2S2*ypb*z2S + eta_t*y2S2*ypb*zqc + eta_t*y2S*ypa*yqd*z2S + eta_t*y2S*ypa*yqd*zqc + eta_t*y2S*ypb*yqd*z2S + eta_t*y2S*ypb*yqd*zqc + y1D2*y2S*z2S + y1D2*y2S*zqc + y1D2*yqd*z2S + y1D2*yqd*zqc + y1D*y2S*ypa*z2S + y1D*y2S*ypa*zqc + y1D*y2S*ypb*z2S + y1D*y2S*ypb*zqc + y1D*ypa*yqd*z2S + y1D*ypa*yqd*zqc + y1D*ypb*yqd*z2S + y1D*ypb*yqd*zqc + y2S*ypa*ypb*z2S + y2S*ypa*ypb*zqc + ypa*ypb*yqd*z2S + ypa*ypb*yqd*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2233) ! | py  py  pz  pz   ( 176) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - 2*eta_p2*eta_q*eta_t*y2S*ypq*zpq2 - 2*eta_p2*eta_q*y1D*ypq*zpq2 - eta_p2*eta_q*ypa*ypq*zpq2 - eta_p2*eta_q*ypb*ypq*zpq2 + eta_p2*eta_t2*y2S2*zpq2 + 2*eta_p2*eta_t*y1D*y2S*zpq2 + eta_p2*eta_t*y2S*ypa*zpq2 + eta_p2*eta_t*y2S*ypb*zpq2 + eta_p2*y1D2*zpq2 + eta_p2*y1D*ypa*zpq2 + eta_p2*y1D*ypb*zpq2 + eta_p2*ypa*ypb*zpq2 + 2*eta_p*eta_q2*ypq2*z2S*zpq + eta_p*eta_q2*ypq2*zpq*zqc + eta_p*eta_q2*ypq2*zpq*zqd - 4*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq*zqc - 2*eta_p*eta_q*eta_t*y2S*ypq*zpq*zqd - 4*eta_p*eta_q*y1D*ypq*z2S*zpq - 2*eta_p*eta_q*y1D*ypq*zpq*zqc - 2*eta_p*eta_q*y1D*ypq*zpq*zqd - 2*eta_p*eta_q*ypa*ypq*z2S*zpq - eta_p*eta_q*ypa*ypq*zpq*zqc - eta_p*eta_q*ypa*ypq*zpq*zqd - 2*eta_p*eta_q*ypb*ypq*z2S*zpq - eta_p*eta_q*ypb*ypq*zpq*zqc - eta_p*eta_q*ypb*ypq*zpq*zqd + 2*eta_p*eta_t2*y2S2*z2S*zpq + eta_p*eta_t2*y2S2*zpq*zqc + eta_p*eta_t2*y2S2*zpq*zqd + 4*eta_p*eta_t*y1D*y2S*z2S*zpq + 2*eta_p*eta_t*y1D*y2S*zpq*zqc + 2*eta_p*eta_t*y1D*y2S*zpq*zqd + 2*eta_p*eta_t*y2S*ypa*z2S*zpq + eta_p*eta_t*y2S*ypa*zpq*zqc + eta_p*eta_t*y2S*ypa*zpq*zqd + 2*eta_p*eta_t*y2S*ypb*z2S*zpq + eta_p*eta_t*y2S*ypb*zpq*zqc + eta_p*eta_t*y2S*ypb*zpq*zqd + 2*eta_p*y1D2*z2S*zpq + eta_p*y1D2*zpq*zqc + eta_p*y1D2*zpq*zqd + 2*eta_p*y1D*ypa*z2S*zpq + eta_p*y1D*ypa*zpq*zqc + eta_p*y1D*ypa*zpq*zqd + 2*eta_p*y1D*ypb*z2S*zpq + eta_p*y1D*ypb*zpq*zqc + eta_p*y1D*ypb*zpq*zqd + 2*eta_p*ypa*ypb*z2S*zpq + eta_p*ypa*ypb*zpq*zqc + eta_p*ypa*ypb*zpq*zqd + eta_q2*ypq2*z2S2 + eta_q2*ypq2*z2S*zqc + eta_q2*ypq2*z2S*zqd + eta_q2*ypq2*zqc*zqd - 2*eta_q*eta_t*y2S*ypq*z2S2 - 2*eta_q*eta_t*y2S*ypq*z2S*zqc - 2*eta_q*eta_t*y2S*ypq*z2S*zqd - 2*eta_q*eta_t*y2S*ypq*zqc*zqd - 2*eta_q*y1D*ypq*z2S2 - 2*eta_q*y1D*ypq*z2S*zqc - 2*eta_q*y1D*ypq*z2S*zqd - 2*eta_q*y1D*ypq*zqc*zqd - eta_q*ypa*ypq*z2S2 - eta_q*ypa*ypq*z2S*zqc - eta_q*ypa*ypq*z2S*zqd - eta_q*ypa*ypq*zqc*zqd - eta_q*ypb*ypq*z2S2 - eta_q*ypb*ypq*z2S*zqc - eta_q*ypb*ypq*z2S*zqd - eta_q*ypb*ypq*zqc*zqd + eta_t2*y2S2*z2S2 + eta_t2*y2S2*z2S*zqc + eta_t2*y2S2*z2S*zqd + eta_t2*y2S2*zqc*zqd + 2*eta_t*y1D*y2S*z2S2 + 2*eta_t*y1D*y2S*z2S*zqc + 2*eta_t*y1D*y2S*z2S*zqd + 2*eta_t*y1D*y2S*zqc*zqd + eta_t*y2S*ypa*z2S2 + eta_t*y2S*ypa*z2S*zqc + eta_t*y2S*ypa*z2S*zqd + eta_t*y2S*ypa*zqc*zqd + eta_t*y2S*ypb*z2S2 + eta_t*y2S*ypb*z2S*zqc + eta_t*y2S*ypb*z2S*zqd + eta_t*y2S*ypb*zqc*zqd + y1D2*z2S2 + y1D2*z2S*zqc + y1D2*z2S*zqd + y1D2*zqc*zqd + y1D*ypa*z2S2 + y1D*ypa*z2S*zqc + y1D*ypa*z2S*zqd + y1D*ypa*zqc*zqd + y1D*ypb*z2S2 + y1D*ypb*z2S*zqc + y1D*ypb*z2S*zqd + y1D*ypb*zqc*zqd + ypa*ypb*z2S2 + ypa*ypb*z2S*zqc + ypa*ypb*z2S*zqd + ypa*ypb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2300) ! | py  pz  s   s    ( 177) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypa*zpq - eta_q*ypq*z1D - eta_q*ypq*zpb + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpb + eta_t*ypa*z2S + y1D*z1D + y1D*zpb + ypa*z1D + ypa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2301) ! | py  pz  s   px   ( 178) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypa*zpq - eta_q*ypq*z1D - eta_q*ypq*zpb + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpb + eta_t*ypa*z2S + y1D*z1D + y1D*zpb + ypa*z1D + ypa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2302) ! | py  pz  s   py   ( 179) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpb + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpb + eta_p*eta_t*ypa*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpb + eta_p*ypa*ypq*z1D + eta_p*ypa*ypq*zpb + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqd*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqd*zpq - eta_q*eta_t*ypq*yqd*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqd*zpq - eta_q*y2S*ypa*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpb - eta_q*ypa*yqd*zpq - eta_q*ypq*yqd*z1D - eta_q*ypq*yqd*zpb + eta_t2*y2S2*z2S + eta_t2*y2S*yqd*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqd*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpb + eta_t*y2S*ypa*z2S + eta_t*y2S*yqd*z1D + eta_t*y2S*yqd*zpb + eta_t*ypa*yqd*z2S + y1D*y2S*z1D + y1D*y2S*zpb + y1D*yqd*z1D + y1D*yqd*zpb + y2S*ypa*z1D + y2S*ypa*zpb + ypa*yqd*z1D + ypa*yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2303) ! | py  pz  s   pz   ( 180) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypa*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpb*zpq + eta_p*eta_t*ypa*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpb*zpq + eta_p*ypa*z1D*zpq + eta_p*ypa*zpb*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqd - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqd - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqd - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqd - eta_q*ypa*z2S*zpq - eta_q*ypa*zpq*zqd - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqd - eta_q*ypq*z2S*zpb - eta_q*ypq*zpb*zqd + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqd + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqd + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqd + eta_t*y2S*z2S*zpb + eta_t*y2S*zpb*zqd + eta_t*ypa*z2S2 + eta_t*ypa*z2S*zqd + y1D*z1D*z2S + y1D*z1D*zqd + y1D*z2S*zpb + y1D*zpb*zqd + ypa*z1D*z2S + ypa*z1D*zqd + ypa*z2S*zpb + ypa*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2310) ! | py  pz  px  s    ( 181) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypa*zpq - eta_q*ypq*z1D - eta_q*ypq*zpb + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpb + eta_t*ypa*z2S + y1D*z1D + y1D*zpb + ypa*z1D + ypa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2311) ! | py  pz  px  px   ( 182) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypa*zpq - eta_q*ypq*z1D - eta_q*ypq*zpb + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpb + eta_t*ypa*z2S + y1D*z1D + y1D*zpb + ypa*z1D + ypa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2312) ! | py  pz  px  py   ( 183) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpb + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpb + eta_p*eta_t*ypa*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpb + eta_p*ypa*ypq*z1D + eta_p*ypa*ypq*zpb + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqd*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqd*zpq - eta_q*eta_t*ypq*yqd*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqd*zpq - eta_q*y2S*ypa*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpb - eta_q*ypa*yqd*zpq - eta_q*ypq*yqd*z1D - eta_q*ypq*yqd*zpb + eta_t2*y2S2*z2S + eta_t2*y2S*yqd*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqd*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpb + eta_t*y2S*ypa*z2S + eta_t*y2S*yqd*z1D + eta_t*y2S*yqd*zpb + eta_t*ypa*yqd*z2S + y1D*y2S*z1D + y1D*y2S*zpb + y1D*yqd*z1D + y1D*yqd*zpb + y2S*ypa*z1D + y2S*ypa*zpb + ypa*yqd*z1D + ypa*yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2313) ! | py  pz  px  pz   ( 184) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypa*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpb*zpq + eta_p*eta_t*ypa*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpb*zpq + eta_p*ypa*z1D*zpq + eta_p*ypa*zpb*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqd - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqd - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqd - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqd - eta_q*ypa*z2S*zpq - eta_q*ypa*zpq*zqd - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqd - eta_q*ypq*z2S*zpb - eta_q*ypq*zpb*zqd + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqd + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqd + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqd + eta_t*y2S*z2S*zpb + eta_t*y2S*zpb*zqd + eta_t*ypa*z2S2 + eta_t*ypa*z2S*zqd + y1D*z1D*z2S + y1D*z1D*zqd + y1D*z2S*zpb + y1D*zpb*zqd + ypa*z1D*z2S + ypa*z1D*zqd + ypa*z2S*zpb + ypa*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2320) ! | py  pz  py  s    ( 185) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpb + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpb + eta_p*eta_t*ypa*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpb + eta_p*ypa*ypq*z1D + eta_p*ypa*ypq*zpb + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqc*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqc*zpq - eta_q*eta_t*ypq*yqc*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqc*zpq - eta_q*y2S*ypa*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpb - eta_q*ypa*yqc*zpq - eta_q*ypq*yqc*z1D - eta_q*ypq*yqc*zpb + eta_t2*y2S2*z2S + eta_t2*y2S*yqc*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqc*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpb + eta_t*y2S*ypa*z2S + eta_t*y2S*yqc*z1D + eta_t*y2S*yqc*zpb + eta_t*ypa*yqc*z2S + y1D*y2S*z1D + y1D*y2S*zpb + y1D*yqc*z1D + y1D*yqc*zpb + y2S*ypa*z1D + y2S*ypa*zpb + ypa*yqc*z1D + ypa*yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2321) ! | py  pz  py  px   ( 186) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypa*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpb + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpb + eta_p*eta_t*ypa*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpb + eta_p*ypa*ypq*z1D + eta_p*ypa*ypq*zpb + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqc*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqc*zpq - eta_q*eta_t*ypq*yqc*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqc*zpq - eta_q*y2S*ypa*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpb - eta_q*ypa*yqc*zpq - eta_q*ypq*yqc*z1D - eta_q*ypq*yqc*zpb + eta_t2*y2S2*z2S + eta_t2*y2S*yqc*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqc*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpb + eta_t*y2S*ypa*z2S + eta_t*y2S*yqc*z1D + eta_t*y2S*yqc*zpb + eta_t*ypa*yqc*z2S + y1D*y2S*z1D + y1D*y2S*zpb + y1D*yqc*z1D + y1D*yqc*zpb + y2S*ypa*z1D + y2S*ypa*zpb + ypa*yqc*z1D + ypa*yqc*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2322) ! | py  pz  py  py   ( 187) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq3*zpq - eta_p2*eta_q*eta_t*y2S*ypq2*zpq - eta_p2*eta_q*eta_t*ypq3*z2S - eta_p2*eta_q*y1D*ypq2*zpq - eta_p2*eta_q*ypa*ypq2*zpq - eta_p2*eta_q*ypq3*z1D - eta_p2*eta_q*ypq3*zpb + eta_p2*eta_t2*y2S*ypq2*z2S + eta_p2*eta_t*y1D*ypq2*z2S + eta_p2*eta_t*y2S*ypq2*z1D + eta_p2*eta_t*y2S*ypq2*zpb + eta_p2*eta_t*ypa*ypq2*z2S + eta_p2*y1D*ypq2*z1D + eta_p2*y1D*ypq2*zpb + eta_p2*ypa*ypq2*z1D + eta_p2*ypa*ypq2*zpb + 2*eta_p*eta_q2*y2S*ypq2*zpq + eta_p*eta_q2*ypq2*yqc*zpq + eta_p*eta_q2*ypq2*yqd*zpq - 2*eta_p*eta_q*eta_t*y2S2*ypq*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq2*z2S - eta_p*eta_q*eta_t*y2S*ypq*yqc*zpq - eta_p*eta_q*eta_t*y2S*ypq*yqd*zpq - eta_p*eta_q*eta_t*ypq2*yqc*z2S - eta_p*eta_q*eta_t*ypq2*yqd*z2S - 2*eta_p*eta_q*y1D*y2S*ypq*zpq - eta_p*eta_q*y1D*ypq*yqc*zpq - eta_p*eta_q*y1D*ypq*yqd*zpq - 2*eta_p*eta_q*y2S*ypa*ypq*zpq - 2*eta_p*eta_q*y2S*ypq2*z1D - 2*eta_p*eta_q*y2S*ypq2*zpb - eta_p*eta_q*ypa*ypq*yqc*zpq - eta_p*eta_q*ypa*ypq*yqd*zpq - eta_p*eta_q*ypq2*yqc*z1D - eta_p*eta_q*ypq2*yqc*zpb - eta_p*eta_q*ypq2*yqd*z1D - eta_p*eta_q*ypq2*yqd*zpb + 2*eta_p*eta_t2*y2S2*ypq*z2S + eta_p*eta_t2*y2S*ypq*yqc*z2S + eta_p*eta_t2*y2S*ypq*yqd*z2S + 2*eta_p*eta_t*y1D*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*yqc*z2S + eta_p*eta_t*y1D*ypq*yqd*z2S + 2*eta_p*eta_t*y2S2*ypq*z1D + 2*eta_p*eta_t*y2S2*ypq*zpb + 2*eta_p*eta_t*y2S*ypa*ypq*z2S + eta_p*eta_t*y2S*ypq*yqc*z1D + eta_p*eta_t*y2S*ypq*yqc*zpb + eta_p*eta_t*y2S*ypq*yqd*z1D + eta_p*eta_t*y2S*ypq*yqd*zpb + eta_p*eta_t*ypa*ypq*yqc*z2S + eta_p*eta_t*ypa*ypq*yqd*z2S + 2*eta_p*y1D*y2S*ypq*z1D + 2*eta_p*y1D*y2S*ypq*zpb + eta_p*y1D*ypq*yqc*z1D + eta_p*y1D*ypq*yqc*zpb + eta_p*y1D*ypq*yqd*z1D + eta_p*y1D*ypq*yqd*zpb + 2*eta_p*y2S*ypa*ypq*z1D + 2*eta_p*y2S*ypa*ypq*zpb + eta_p*ypa*ypq*yqc*z1D + eta_p*ypa*ypq*yqc*zpb + eta_p*ypa*ypq*yqd*z1D + eta_p*ypa*ypq*yqd*zpb + eta_q2*y2S2*ypq*zpq + eta_q2*y2S*ypq*yqc*zpq + eta_q2*y2S*ypq*yqd*zpq + eta_q2*ypq*yqc*yqd*zpq - eta_q*eta_t*y2S3*zpq - eta_q*eta_t*y2S2*ypq*z2S - eta_q*eta_t*y2S2*yqc*zpq - eta_q*eta_t*y2S2*yqd*zpq - eta_q*eta_t*y2S*ypq*yqc*z2S - eta_q*eta_t*y2S*ypq*yqd*z2S - eta_q*eta_t*y2S*yqc*yqd*zpq - eta_q*eta_t*ypq*yqc*yqd*z2S - eta_q*y1D*y2S2*zpq - eta_q*y1D*y2S*yqc*zpq - eta_q*y1D*y2S*yqd*zpq - eta_q*y1D*yqc*yqd*zpq - eta_q*y2S2*ypa*zpq - eta_q*y2S2*ypq*z1D - eta_q*y2S2*ypq*zpb - eta_q*y2S*ypa*yqc*zpq - eta_q*y2S*ypa*yqd*zpq - eta_q*y2S*ypq*yqc*z1D - eta_q*y2S*ypq*yqc*zpb - eta_q*y2S*ypq*yqd*z1D - eta_q*y2S*ypq*yqd*zpb - eta_q*ypa*yqc*yqd*zpq - eta_q*ypq*yqc*yqd*z1D - eta_q*ypq*yqc*yqd*zpb + eta_t2*y2S3*z2S + eta_t2*y2S2*yqc*z2S + eta_t2*y2S2*yqd*z2S + eta_t2*y2S*yqc*yqd*z2S + eta_t*y1D*y2S2*z2S + eta_t*y1D*y2S*yqc*z2S + eta_t*y1D*y2S*yqd*z2S + eta_t*y1D*yqc*yqd*z2S + eta_t*y2S3*z1D + eta_t*y2S3*zpb + eta_t*y2S2*ypa*z2S + eta_t*y2S2*yqc*z1D + eta_t*y2S2*yqc*zpb + eta_t*y2S2*yqd*z1D + eta_t*y2S2*yqd*zpb + eta_t*y2S*ypa*yqc*z2S + eta_t*y2S*ypa*yqd*z2S + eta_t*y2S*yqc*yqd*z1D + eta_t*y2S*yqc*yqd*zpb + eta_t*ypa*yqc*yqd*z2S + y1D*y2S2*z1D + y1D*y2S2*zpb + y1D*y2S*yqc*z1D + y1D*y2S*yqc*zpb + y1D*y2S*yqd*z1D + y1D*y2S*yqd*zpb + y1D*yqc*yqd*z1D + y1D*yqc*yqd*zpb + y2S2*ypa*z1D + y2S2*ypa*zpb + y2S*ypa*yqc*z1D + y2S*ypa*yqc*zpb + y2S*ypa*yqd*z1D + y2S*ypa*yqd*zpb + ypa*yqc*yqd*z1D + ypa*yqc*yqd*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2323) ! | py  pz  py  pz   ( 188) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - eta_p2*eta_q*eta_t*y2S*ypq*zpq2 - eta_p2*eta_q*eta_t*ypq2*z2S*zpq - eta_p2*eta_q*y1D*ypq*zpq2 - eta_p2*eta_q*ypa*ypq*zpq2 - eta_p2*eta_q*ypq2*z1D*zpq - eta_p2*eta_q*ypq2*zpb*zpq + eta_p2*eta_t2*y2S*ypq*z2S*zpq + eta_p2*eta_t*y1D*ypq*z2S*zpq + eta_p2*eta_t*y2S*ypq*z1D*zpq + eta_p2*eta_t*y2S*ypq*zpb*zpq + eta_p2*eta_t*ypa*ypq*z2S*zpq + eta_p2*y1D*ypq*z1D*zpq + eta_p2*y1D*ypq*zpb*zpq + eta_p2*ypa*ypq*z1D*zpq + eta_p2*ypa*ypq*zpb*zpq + eta_p*eta_q2*y2S*ypq*zpq2 + eta_p*eta_q2*ypq2*z2S*zpq + eta_p*eta_q2*ypq2*zpq*zqd + eta_p*eta_q2*ypq*yqc*zpq2 - eta_p*eta_q*eta_t*y2S2*zpq2 - 2*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq*zqd - eta_p*eta_q*eta_t*y2S*yqc*zpq2 - eta_p*eta_q*eta_t*ypq2*z2S2 - eta_p*eta_q*eta_t*ypq2*z2S*zqd - eta_p*eta_q*eta_t*ypq*yqc*z2S*zpq - eta_p*eta_q*y1D*y2S*zpq2 - eta_p*eta_q*y1D*ypq*z2S*zpq - eta_p*eta_q*y1D*ypq*zpq*zqd - eta_p*eta_q*y1D*yqc*zpq2 - eta_p*eta_q*y2S*ypa*zpq2 - eta_p*eta_q*y2S*ypq*z1D*zpq - eta_p*eta_q*y2S*ypq*zpb*zpq - eta_p*eta_q*ypa*ypq*z2S*zpq - eta_p*eta_q*ypa*ypq*zpq*zqd - eta_p*eta_q*ypa*yqc*zpq2 - eta_p*eta_q*ypq2*z1D*z2S - eta_p*eta_q*ypq2*z1D*zqd - eta_p*eta_q*ypq2*z2S*zpb - eta_p*eta_q*ypq2*zpb*zqd - eta_p*eta_q*ypq*yqc*z1D*zpq - eta_p*eta_q*ypq*yqc*zpb*zpq + eta_p*eta_t2*y2S2*z2S*zpq + eta_p*eta_t2*y2S*ypq*z2S2 + eta_p*eta_t2*y2S*ypq*z2S*zqd + eta_p*eta_t2*y2S*yqc*z2S*zpq + eta_p*eta_t*y1D*y2S*z2S*zpq + eta_p*eta_t*y1D*ypq*z2S2 + eta_p*eta_t*y1D*ypq*z2S*zqd + eta_p*eta_t*y1D*yqc*z2S*zpq + eta_p*eta_t*y2S2*z1D*zpq + eta_p*eta_t*y2S2*zpb*zpq + eta_p*eta_t*y2S*ypa*z2S*zpq + eta_p*eta_t*y2S*ypq*z1D*z2S + eta_p*eta_t*y2S*ypq*z1D*zqd + eta_p*eta_t*y2S*ypq*z2S*zpb + eta_p*eta_t*y2S*ypq*zpb*zqd + eta_p*eta_t*y2S*yqc*z1D*zpq + eta_p*eta_t*y2S*yqc*zpb*zpq + eta_p*eta_t*ypa*ypq*z2S2 + eta_p*eta_t*ypa*ypq*z2S*zqd + eta_p*eta_t*ypa*yqc*z2S*zpq + eta_p*y1D*y2S*z1D*zpq + eta_p*y1D*y2S*zpb*zpq + eta_p*y1D*ypq*z1D*z2S + eta_p*y1D*ypq*z1D*zqd + eta_p*y1D*ypq*z2S*zpb + eta_p*y1D*ypq*zpb*zqd + eta_p*y1D*yqc*z1D*zpq + eta_p*y1D*yqc*zpb*zpq + eta_p*y2S*ypa*z1D*zpq + eta_p*y2S*ypa*zpb*zpq + eta_p*ypa*ypq*z1D*z2S + eta_p*ypa*ypq*z1D*zqd + eta_p*ypa*ypq*z2S*zpb + eta_p*ypa*ypq*zpb*zqd + eta_p*ypa*yqc*z1D*zpq + eta_p*ypa*yqc*zpb*zpq + eta_q2*y2S*ypq*z2S*zpq + eta_q2*y2S*ypq*zpq*zqd + eta_q2*ypq*yqc*z2S*zpq + eta_q2*ypq*yqc*zpq*zqd - eta_q*eta_t*y2S2*z2S*zpq - eta_q*eta_t*y2S2*zpq*zqd - eta_q*eta_t*y2S*ypq*z2S2 - eta_q*eta_t*y2S*ypq*z2S*zqd - eta_q*eta_t*y2S*yqc*z2S*zpq - eta_q*eta_t*y2S*yqc*zpq*zqd - eta_q*eta_t*ypq*yqc*z2S2 - eta_q*eta_t*ypq*yqc*z2S*zqd - eta_q*y1D*y2S*z2S*zpq - eta_q*y1D*y2S*zpq*zqd - eta_q*y1D*yqc*z2S*zpq - eta_q*y1D*yqc*zpq*zqd - eta_q*y2S*ypa*z2S*zpq - eta_q*y2S*ypa*zpq*zqd - eta_q*y2S*ypq*z1D*z2S - eta_q*y2S*ypq*z1D*zqd - eta_q*y2S*ypq*z2S*zpb - eta_q*y2S*ypq*zpb*zqd - eta_q*ypa*yqc*z2S*zpq - eta_q*ypa*yqc*zpq*zqd - eta_q*ypq*yqc*z1D*z2S - eta_q*ypq*yqc*z1D*zqd - eta_q*ypq*yqc*z2S*zpb - eta_q*ypq*yqc*zpb*zqd + eta_t2*y2S2*z2S2 + eta_t2*y2S2*z2S*zqd + eta_t2*y2S*yqc*z2S2 + eta_t2*y2S*yqc*z2S*zqd + eta_t*y1D*y2S*z2S2 + eta_t*y1D*y2S*z2S*zqd + eta_t*y1D*yqc*z2S2 + eta_t*y1D*yqc*z2S*zqd + eta_t*y2S2*z1D*z2S + eta_t*y2S2*z1D*zqd + eta_t*y2S2*z2S*zpb + eta_t*y2S2*zpb*zqd + eta_t*y2S*ypa*z2S2 + eta_t*y2S*ypa*z2S*zqd + eta_t*y2S*yqc*z1D*z2S + eta_t*y2S*yqc*z1D*zqd + eta_t*y2S*yqc*z2S*zpb + eta_t*y2S*yqc*zpb*zqd + eta_t*ypa*yqc*z2S2 + eta_t*ypa*yqc*z2S*zqd + y1D*y2S*z1D*z2S + y1D*y2S*z1D*zqd + y1D*y2S*z2S*zpb + y1D*y2S*zpb*zqd + y1D*yqc*z1D*z2S + y1D*yqc*z1D*zqd + y1D*yqc*z2S*zpb + y1D*yqc*zpb*zqd + y2S*ypa*z1D*z2S + y2S*ypa*z1D*zqd + y2S*ypa*z2S*zpb + y2S*ypa*zpb*zqd + ypa*yqc*z1D*z2S + ypa*yqc*z1D*zqd + ypa*yqc*z2S*zpb + ypa*yqc*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2330) ! | py  pz  pz  s    ( 189) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypa*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpb*zpq + eta_p*eta_t*ypa*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpb*zpq + eta_p*ypa*z1D*zpq + eta_p*ypa*zpb*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqc - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqc - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqc - eta_q*ypa*z2S*zpq - eta_q*ypa*zpq*zqc - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqc - eta_q*ypq*z2S*zpb - eta_q*ypq*zpb*zqc + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqc + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqc + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqc + eta_t*y2S*z2S*zpb + eta_t*y2S*zpb*zqc + eta_t*ypa*z2S2 + eta_t*ypa*z2S*zqc + y1D*z1D*z2S + y1D*z1D*zqc + y1D*z2S*zpb + y1D*zpb*zqc + ypa*z1D*z2S + ypa*z1D*zqc + ypa*z2S*zpb + ypa*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2331) ! | py  pz  pz  px   ( 190) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypa*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpb*zpq + eta_p*eta_t*ypa*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpb*zpq + eta_p*ypa*z1D*zpq + eta_p*ypa*zpb*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqc - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqc - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqc - eta_q*ypa*z2S*zpq - eta_q*ypa*zpq*zqc - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqc - eta_q*ypq*z2S*zpb - eta_q*ypq*zpb*zqc + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqc + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqc + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqc + eta_t*y2S*z2S*zpb + eta_t*y2S*zpb*zqc + eta_t*ypa*z2S2 + eta_t*ypa*z2S*zqc + y1D*z1D*z2S + y1D*z1D*zqc + y1D*z2S*zpb + y1D*zpb*zqc + ypa*z1D*z2S + ypa*z1D*zqc + ypa*z2S*zpb + ypa*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2332) ! | py  pz  pz  py   ( 191) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - eta_p2*eta_q*eta_t*y2S*ypq*zpq2 - eta_p2*eta_q*eta_t*ypq2*z2S*zpq - eta_p2*eta_q*y1D*ypq*zpq2 - eta_p2*eta_q*ypa*ypq*zpq2 - eta_p2*eta_q*ypq2*z1D*zpq - eta_p2*eta_q*ypq2*zpb*zpq + eta_p2*eta_t2*y2S*ypq*z2S*zpq + eta_p2*eta_t*y1D*ypq*z2S*zpq + eta_p2*eta_t*y2S*ypq*z1D*zpq + eta_p2*eta_t*y2S*ypq*zpb*zpq + eta_p2*eta_t*ypa*ypq*z2S*zpq + eta_p2*y1D*ypq*z1D*zpq + eta_p2*y1D*ypq*zpb*zpq + eta_p2*ypa*ypq*z1D*zpq + eta_p2*ypa*ypq*zpb*zpq + eta_p*eta_q2*y2S*ypq*zpq2 + eta_p*eta_q2*ypq2*z2S*zpq + eta_p*eta_q2*ypq2*zpq*zqc + eta_p*eta_q2*ypq*yqd*zpq2 - eta_p*eta_q*eta_t*y2S2*zpq2 - 2*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq*zqc - eta_p*eta_q*eta_t*y2S*yqd*zpq2 - eta_p*eta_q*eta_t*ypq2*z2S2 - eta_p*eta_q*eta_t*ypq2*z2S*zqc - eta_p*eta_q*eta_t*ypq*yqd*z2S*zpq - eta_p*eta_q*y1D*y2S*zpq2 - eta_p*eta_q*y1D*ypq*z2S*zpq - eta_p*eta_q*y1D*ypq*zpq*zqc - eta_p*eta_q*y1D*yqd*zpq2 - eta_p*eta_q*y2S*ypa*zpq2 - eta_p*eta_q*y2S*ypq*z1D*zpq - eta_p*eta_q*y2S*ypq*zpb*zpq - eta_p*eta_q*ypa*ypq*z2S*zpq - eta_p*eta_q*ypa*ypq*zpq*zqc - eta_p*eta_q*ypa*yqd*zpq2 - eta_p*eta_q*ypq2*z1D*z2S - eta_p*eta_q*ypq2*z1D*zqc - eta_p*eta_q*ypq2*z2S*zpb - eta_p*eta_q*ypq2*zpb*zqc - eta_p*eta_q*ypq*yqd*z1D*zpq - eta_p*eta_q*ypq*yqd*zpb*zpq + eta_p*eta_t2*y2S2*z2S*zpq + eta_p*eta_t2*y2S*ypq*z2S2 + eta_p*eta_t2*y2S*ypq*z2S*zqc + eta_p*eta_t2*y2S*yqd*z2S*zpq + eta_p*eta_t*y1D*y2S*z2S*zpq + eta_p*eta_t*y1D*ypq*z2S2 + eta_p*eta_t*y1D*ypq*z2S*zqc + eta_p*eta_t*y1D*yqd*z2S*zpq + eta_p*eta_t*y2S2*z1D*zpq + eta_p*eta_t*y2S2*zpb*zpq + eta_p*eta_t*y2S*ypa*z2S*zpq + eta_p*eta_t*y2S*ypq*z1D*z2S + eta_p*eta_t*y2S*ypq*z1D*zqc + eta_p*eta_t*y2S*ypq*z2S*zpb + eta_p*eta_t*y2S*ypq*zpb*zqc + eta_p*eta_t*y2S*yqd*z1D*zpq + eta_p*eta_t*y2S*yqd*zpb*zpq + eta_p*eta_t*ypa*ypq*z2S2 + eta_p*eta_t*ypa*ypq*z2S*zqc + eta_p*eta_t*ypa*yqd*z2S*zpq + eta_p*y1D*y2S*z1D*zpq + eta_p*y1D*y2S*zpb*zpq + eta_p*y1D*ypq*z1D*z2S + eta_p*y1D*ypq*z1D*zqc + eta_p*y1D*ypq*z2S*zpb + eta_p*y1D*ypq*zpb*zqc + eta_p*y1D*yqd*z1D*zpq + eta_p*y1D*yqd*zpb*zpq + eta_p*y2S*ypa*z1D*zpq + eta_p*y2S*ypa*zpb*zpq + eta_p*ypa*ypq*z1D*z2S + eta_p*ypa*ypq*z1D*zqc + eta_p*ypa*ypq*z2S*zpb + eta_p*ypa*ypq*zpb*zqc + eta_p*ypa*yqd*z1D*zpq + eta_p*ypa*yqd*zpb*zpq + eta_q2*y2S*ypq*z2S*zpq + eta_q2*y2S*ypq*zpq*zqc + eta_q2*ypq*yqd*z2S*zpq + eta_q2*ypq*yqd*zpq*zqc - eta_q*eta_t*y2S2*z2S*zpq - eta_q*eta_t*y2S2*zpq*zqc - eta_q*eta_t*y2S*ypq*z2S2 - eta_q*eta_t*y2S*ypq*z2S*zqc - eta_q*eta_t*y2S*yqd*z2S*zpq - eta_q*eta_t*y2S*yqd*zpq*zqc - eta_q*eta_t*ypq*yqd*z2S2 - eta_q*eta_t*ypq*yqd*z2S*zqc - eta_q*y1D*y2S*z2S*zpq - eta_q*y1D*y2S*zpq*zqc - eta_q*y1D*yqd*z2S*zpq - eta_q*y1D*yqd*zpq*zqc - eta_q*y2S*ypa*z2S*zpq - eta_q*y2S*ypa*zpq*zqc - eta_q*y2S*ypq*z1D*z2S - eta_q*y2S*ypq*z1D*zqc - eta_q*y2S*ypq*z2S*zpb - eta_q*y2S*ypq*zpb*zqc - eta_q*ypa*yqd*z2S*zpq - eta_q*ypa*yqd*zpq*zqc - eta_q*ypq*yqd*z1D*z2S - eta_q*ypq*yqd*z1D*zqc - eta_q*ypq*yqd*z2S*zpb - eta_q*ypq*yqd*zpb*zqc + eta_t2*y2S2*z2S2 + eta_t2*y2S2*z2S*zqc + eta_t2*y2S*yqd*z2S2 + eta_t2*y2S*yqd*z2S*zqc + eta_t*y1D*y2S*z2S2 + eta_t*y1D*y2S*z2S*zqc + eta_t*y1D*yqd*z2S2 + eta_t*y1D*yqd*z2S*zqc + eta_t*y2S2*z1D*z2S + eta_t*y2S2*z1D*zqc + eta_t*y2S2*z2S*zpb + eta_t*y2S2*zpb*zqc + eta_t*y2S*ypa*z2S2 + eta_t*y2S*ypa*z2S*zqc + eta_t*y2S*yqd*z1D*z2S + eta_t*y2S*yqd*z1D*zqc + eta_t*y2S*yqd*z2S*zpb + eta_t*y2S*yqd*zpb*zqc + eta_t*ypa*yqd*z2S2 + eta_t*ypa*yqd*z2S*zqc + y1D*y2S*z1D*z2S + y1D*y2S*z1D*zqc + y1D*y2S*z2S*zpb + y1D*y2S*zpb*zqc + y1D*yqd*z1D*z2S + y1D*yqd*z1D*zqc + y1D*yqd*z2S*zpb + y1D*yqd*zpb*zqc + y2S*ypa*z1D*z2S + y2S*ypa*z1D*zqc + y2S*ypa*z2S*zpb + y2S*ypa*zpb*zqc + ypa*yqd*z1D*z2S + ypa*yqd*z1D*zqc + ypa*yqd*z2S*zpb + ypa*yqd*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (2333) ! | py  pz  pz  pz   ( 192) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq*zpq3 - eta_p2*eta_q*eta_t*y2S*zpq3 - eta_p2*eta_q*eta_t*ypq*z2S*zpq2 - eta_p2*eta_q*y1D*zpq3 - eta_p2*eta_q*ypa*zpq3 - eta_p2*eta_q*ypq*z1D*zpq2 - eta_p2*eta_q*ypq*zpb*zpq2 + eta_p2*eta_t2*y2S*z2S*zpq2 + eta_p2*eta_t*y1D*z2S*zpq2 + eta_p2*eta_t*y2S*z1D*zpq2 + eta_p2*eta_t*y2S*zpb*zpq2 + eta_p2*eta_t*ypa*z2S*zpq2 + eta_p2*y1D*z1D*zpq2 + eta_p2*y1D*zpb*zpq2 + eta_p2*ypa*z1D*zpq2 + eta_p2*ypa*zpb*zpq2 + 2*eta_p*eta_q2*ypq*z2S*zpq2 + eta_p*eta_q2*ypq*zpq2*zqc + eta_p*eta_q2*ypq*zpq2*zqd - 2*eta_p*eta_q*eta_t*y2S*z2S*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2*zqc - eta_p*eta_q*eta_t*y2S*zpq2*zqd - 2*eta_p*eta_q*eta_t*ypq*z2S2*zpq - eta_p*eta_q*eta_t*ypq*z2S*zpq*zqc - eta_p*eta_q*eta_t*ypq*z2S*zpq*zqd - 2*eta_p*eta_q*y1D*z2S*zpq2 - eta_p*eta_q*y1D*zpq2*zqc - eta_p*eta_q*y1D*zpq2*zqd - 2*eta_p*eta_q*ypa*z2S*zpq2 - eta_p*eta_q*ypa*zpq2*zqc - eta_p*eta_q*ypa*zpq2*zqd - 2*eta_p*eta_q*ypq*z1D*z2S*zpq - eta_p*eta_q*ypq*z1D*zpq*zqc - eta_p*eta_q*ypq*z1D*zpq*zqd - 2*eta_p*eta_q*ypq*z2S*zpb*zpq - eta_p*eta_q*ypq*zpb*zpq*zqc - eta_p*eta_q*ypq*zpb*zpq*zqd + 2*eta_p*eta_t2*y2S*z2S2*zpq + eta_p*eta_t2*y2S*z2S*zpq*zqc + eta_p*eta_t2*y2S*z2S*zpq*zqd + 2*eta_p*eta_t*y1D*z2S2*zpq + eta_p*eta_t*y1D*z2S*zpq*zqc + eta_p*eta_t*y1D*z2S*zpq*zqd + 2*eta_p*eta_t*y2S*z1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq*zqc + eta_p*eta_t*y2S*z1D*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpb*zpq + eta_p*eta_t*y2S*zpb*zpq*zqc + eta_p*eta_t*y2S*zpb*zpq*zqd + 2*eta_p*eta_t*ypa*z2S2*zpq + eta_p*eta_t*ypa*z2S*zpq*zqc + eta_p*eta_t*ypa*z2S*zpq*zqd + 2*eta_p*y1D*z1D*z2S*zpq + eta_p*y1D*z1D*zpq*zqc + eta_p*y1D*z1D*zpq*zqd + 2*eta_p*y1D*z2S*zpb*zpq + eta_p*y1D*zpb*zpq*zqc + eta_p*y1D*zpb*zpq*zqd + 2*eta_p*ypa*z1D*z2S*zpq + eta_p*ypa*z1D*zpq*zqc + eta_p*ypa*z1D*zpq*zqd + 2*eta_p*ypa*z2S*zpb*zpq + eta_p*ypa*zpb*zpq*zqc + eta_p*ypa*zpb*zpq*zqd + eta_q2*ypq*z2S2*zpq + eta_q2*ypq*z2S*zpq*zqc + eta_q2*ypq*z2S*zpq*zqd + eta_q2*ypq*zpq*zqc*zqd - eta_q*eta_t*y2S*z2S2*zpq - eta_q*eta_t*y2S*z2S*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq*zqd - eta_q*eta_t*y2S*zpq*zqc*zqd - eta_q*eta_t*ypq*z2S3 - eta_q*eta_t*ypq*z2S2*zqc - eta_q*eta_t*ypq*z2S2*zqd - eta_q*eta_t*ypq*z2S*zqc*zqd - eta_q*y1D*z2S2*zpq - eta_q*y1D*z2S*zpq*zqc - eta_q*y1D*z2S*zpq*zqd - eta_q*y1D*zpq*zqc*zqd - eta_q*ypa*z2S2*zpq - eta_q*ypa*z2S*zpq*zqc - eta_q*ypa*z2S*zpq*zqd - eta_q*ypa*zpq*zqc*zqd - eta_q*ypq*z1D*z2S2 - eta_q*ypq*z1D*z2S*zqc - eta_q*ypq*z1D*z2S*zqd - eta_q*ypq*z1D*zqc*zqd - eta_q*ypq*z2S2*zpb - eta_q*ypq*z2S*zpb*zqc - eta_q*ypq*z2S*zpb*zqd - eta_q*ypq*zpb*zqc*zqd + eta_t2*y2S*z2S3 + eta_t2*y2S*z2S2*zqc + eta_t2*y2S*z2S2*zqd + eta_t2*y2S*z2S*zqc*zqd + eta_t*y1D*z2S3 + eta_t*y1D*z2S2*zqc + eta_t*y1D*z2S2*zqd + eta_t*y1D*z2S*zqc*zqd + eta_t*y2S*z1D*z2S2 + eta_t*y2S*z1D*z2S*zqc + eta_t*y2S*z1D*z2S*zqd + eta_t*y2S*z1D*zqc*zqd + eta_t*y2S*z2S2*zpb + eta_t*y2S*z2S*zpb*zqc + eta_t*y2S*z2S*zpb*zqd + eta_t*y2S*zpb*zqc*zqd + eta_t*ypa*z2S3 + eta_t*ypa*z2S2*zqc + eta_t*ypa*z2S2*zqd + eta_t*ypa*z2S*zqc*zqd + y1D*z1D*z2S2 + y1D*z1D*z2S*zqc + y1D*z1D*z2S*zqd + y1D*z1D*zqc*zqd + y1D*z2S2*zpb + y1D*z2S*zpb*zqc + y1D*z2S*zpb*zqd + y1D*zpb*zqc*zqd + ypa*z1D*z2S2 + ypa*z1D*z2S*zqc + ypa*z1D*z2S*zqd + ypa*z1D*zqc*zqd + ypa*z2S2*zpb + ypa*z2S*zpb*zqc + ypa*z2S*zpb*zqd + ypa*zpb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3000) ! | pz  s   s   s    ( 193) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3001) ! | pz  s   s   px   ( 194) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3002) ! | pz  s   s   py   ( 195) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpa + yqd*z1D + yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3003) ! | pz  s   s   pz   ( 196) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpa + zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3010) ! | pz  s   px  s    ( 197) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3011) ! | pz  s   px  px   ( 198) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3012) ! | pz  s   px  py   ( 199) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpa + yqd*z1D + yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3013) ! | pz  s   px  pz   ( 200) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpa + zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3020) ! | pz  s   py  s    ( 201) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpa + yqc*z1D + yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3021) ! | pz  s   py  px   ( 202) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpa + yqc*z1D + yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3022) ! | pz  s   py  py   ( 203) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*ypq2*z2S + eta_p2*ypq2*z1D + eta_p2*ypq2*zpa - 2*eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq*yqc*zpq - eta_p*eta_q*ypq*yqd*zpq + 2*eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*ypq*yqc*z2S + eta_p*eta_t*ypq*yqd*z2S + 2*eta_p*y2S*ypq*z1D + 2*eta_p*y2S*ypq*zpa + eta_p*ypq*yqc*z1D + eta_p*ypq*yqc*zpa + eta_p*ypq*yqd*z1D + eta_p*ypq*yqd*zpa - eta_q*y2S2*zpq - eta_q*y2S*yqc*zpq - eta_q*y2S*yqd*zpq - eta_q*yqc*yqd*zpq + eta_t*y2S2*z2S + eta_t*y2S*yqc*z2S + eta_t*y2S*yqd*z2S + eta_t*yqc*yqd*z2S + y2S2*z1D + y2S2*zpa + y2S*yqc*z1D + y2S*yqc*zpa + y2S*yqd*z1D + y2S*yqd*zpa + yqc*yqd*z1D + yqc*yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3023) ! | pz  s   py  pz   ( 204) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpa*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqd - eta_p*eta_q*yqc*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqd + eta_p*eta_t*yqc*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpa*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqd + eta_p*ypq*z2S*zpa + eta_p*ypq*zpa*zqd + eta_p*yqc*z1D*zpq + eta_p*yqc*zpa*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqd - eta_q*yqc*z2S*zpq - eta_q*yqc*zpq*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqd + eta_t*yqc*z2S2 + eta_t*yqc*z2S*zqd + y2S*z1D*z2S + y2S*z1D*zqd + y2S*z2S*zpa + y2S*zpa*zqd + yqc*z1D*z2S + yqc*z1D*zqd + yqc*z2S*zpa + yqc*zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3030) ! | pz  s   pz  s    ( 205) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpa + zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3031) ! | pz  s   pz  px   ( 206) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpa + zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3032) ! | pz  s   pz  py   ( 207) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpa*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*yqd*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqc + eta_p*eta_t*yqd*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpa*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqc + eta_p*ypq*z2S*zpa + eta_p*ypq*zpa*zqc + eta_p*yqd*z1D*zpq + eta_p*yqd*zpa*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqc - eta_q*yqd*z2S*zpq - eta_q*yqd*zpq*zqc + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*yqd*z2S2 + eta_t*yqd*z2S*zqc + y2S*z1D*z2S + y2S*z1D*zqc + y2S*z2S*zpa + y2S*zpa*zqc + yqd*z1D*z2S + yqd*z1D*zqc + yqd*z2S*zpa + yqd*zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3033) ! | pz  s   pz  pz   ( 208) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*zpq3 + eta_p2*eta_t*z2S*zpq2 + eta_p2*z1D*zpq2 + eta_p2*zpa*zpq2 - 2*eta_p*eta_q*z2S*zpq2 - eta_p*eta_q*zpq2*zqc - eta_p*eta_q*zpq2*zqd + 2*eta_p*eta_t*z2S2*zpq + eta_p*eta_t*z2S*zpq*zqc + eta_p*eta_t*z2S*zpq*zqd + 2*eta_p*z1D*z2S*zpq + eta_p*z1D*zpq*zqc + eta_p*z1D*zpq*zqd + 2*eta_p*z2S*zpa*zpq + eta_p*zpa*zpq*zqc + eta_p*zpa*zpq*zqd - eta_q*z2S2*zpq - eta_q*z2S*zpq*zqc - eta_q*z2S*zpq*zqd - eta_q*zpq*zqc*zqd + eta_t*z2S3 + eta_t*z2S2*zqc + eta_t*z2S2*zqd + eta_t*z2S*zqc*zqd + z1D*z2S2 + z1D*z2S*zqc + z1D*z2S*zqd + z1D*zqc*zqd + z2S2*zpa + z2S*zpa*zqc + z2S*zpa*zqd + zpa*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3100) ! | pz  px  s   s    ( 209) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3101) ! | pz  px  s   px   ( 210) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3102) ! | pz  px  s   py   ( 211) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpa + yqd*z1D + yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3103) ! | pz  px  s   pz   ( 212) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpa + zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3110) ! | pz  px  px  s    ( 213) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3111) ! | pz  px  px  px   ( 214) 
      
      n         = 0
      const     = pi2 * D2 * (z1D + zpa + eta_t * z2S - eta_q * zpq) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3112) ! | pz  px  px  py   ( 215) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqd*zpq + eta_t*y2S*z2S + eta_t*yqd*z2S + y2S*z1D + y2S*zpa + yqd*z1D + yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3113) ! | pz  px  px  pz   ( 216) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqd + eta_t*z2S2 + eta_t*z2S*zqd + z1D*z2S + z1D*zqd + z2S*zpa + zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3120) ! | pz  px  py  s    ( 217) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpa + yqc*z1D + yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3121) ! | pz  px  py  px   ( 218) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*ypq*zpq + eta_p*eta_t*ypq*z2S + eta_p*ypq*z1D + eta_p*ypq*zpa - eta_q*y2S*zpq - eta_q*yqc*zpq + eta_t*y2S*z2S + eta_t*yqc*z2S + y2S*z1D + y2S*zpa + yqc*z1D + yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3122) ! | pz  px  py  py   ( 219) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq2*zpq + eta_p2*eta_t*ypq2*z2S + eta_p2*ypq2*z1D + eta_p2*ypq2*zpa - 2*eta_p*eta_q*y2S*ypq*zpq - eta_p*eta_q*ypq*yqc*zpq - eta_p*eta_q*ypq*yqd*zpq + 2*eta_p*eta_t*y2S*ypq*z2S + eta_p*eta_t*ypq*yqc*z2S + eta_p*eta_t*ypq*yqd*z2S + 2*eta_p*y2S*ypq*z1D + 2*eta_p*y2S*ypq*zpa + eta_p*ypq*yqc*z1D + eta_p*ypq*yqc*zpa + eta_p*ypq*yqd*z1D + eta_p*ypq*yqd*zpa - eta_q*y2S2*zpq - eta_q*y2S*yqc*zpq - eta_q*y2S*yqd*zpq - eta_q*yqc*yqd*zpq + eta_t*y2S2*z2S + eta_t*y2S*yqc*z2S + eta_t*y2S*yqd*z2S + eta_t*yqc*yqd*z2S + y2S2*z1D + y2S2*zpa + y2S*yqc*z1D + y2S*yqc*zpa + y2S*yqd*z1D + y2S*yqd*zpa + yqc*yqd*z1D + yqc*yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3123) ! | pz  px  py  pz   ( 220) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpa*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqd - eta_p*eta_q*yqc*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqd + eta_p*eta_t*yqc*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpa*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqd + eta_p*ypq*z2S*zpa + eta_p*ypq*zpa*zqd + eta_p*yqc*z1D*zpq + eta_p*yqc*zpa*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqd - eta_q*yqc*z2S*zpq - eta_q*yqc*zpq*zqd + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqd + eta_t*yqc*z2S2 + eta_t*yqc*z2S*zqd + y2S*z1D*z2S + y2S*z1D*zqd + y2S*z2S*zpa + y2S*zpa*zqd + yqc*z1D*z2S + yqc*z1D*zqd + yqc*z2S*zpa + yqc*zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3130) ! | pz  px  pz  s    ( 221) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpa + zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3131) ! | pz  px  pz  px   ( 222) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p*eta_q*zpq2 + eta_p*eta_t*z2S*zpq + eta_p*z1D*zpq + eta_p*zpa*zpq - eta_q*z2S*zpq - eta_q*zpq*zqc + eta_t*z2S2 + eta_t*z2S*zqc + z1D*z2S + z1D*zqc + z2S*zpa + zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3132) ! | pz  px  pz  py   ( 223) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*ypq*zpq2 + eta_p2*eta_t*ypq*z2S*zpq + eta_p2*ypq*z1D*zpq + eta_p2*ypq*zpa*zpq - eta_p*eta_q*y2S*zpq2 - eta_p*eta_q*ypq*z2S*zpq - eta_p*eta_q*ypq*zpq*zqc - eta_p*eta_q*yqd*zpq2 + eta_p*eta_t*y2S*z2S*zpq + eta_p*eta_t*ypq*z2S2 + eta_p*eta_t*ypq*z2S*zqc + eta_p*eta_t*yqd*z2S*zpq + eta_p*y2S*z1D*zpq + eta_p*y2S*zpa*zpq + eta_p*ypq*z1D*z2S + eta_p*ypq*z1D*zqc + eta_p*ypq*z2S*zpa + eta_p*ypq*zpa*zqc + eta_p*yqd*z1D*zpq + eta_p*yqd*zpa*zpq - eta_q*y2S*z2S*zpq - eta_q*y2S*zpq*zqc - eta_q*yqd*z2S*zpq - eta_q*yqd*zpq*zqc + eta_t*y2S*z2S2 + eta_t*y2S*z2S*zqc + eta_t*yqd*z2S2 + eta_t*yqd*z2S*zqc + y2S*z1D*z2S + y2S*z1D*zqc + y2S*z2S*zpa + y2S*zpa*zqc + yqd*z1D*z2S + yqd*z1D*zqc + yqd*z2S*zpa + yqd*zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3133) ! | pz  px  pz  pz   ( 224) 
      
      n         = 0
      const     = pi2 * D2 * (-eta_p2*eta_q*zpq3 + eta_p2*eta_t*z2S*zpq2 + eta_p2*z1D*zpq2 + eta_p2*zpa*zpq2 - 2*eta_p*eta_q*z2S*zpq2 - eta_p*eta_q*zpq2*zqc - eta_p*eta_q*zpq2*zqd + 2*eta_p*eta_t*z2S2*zpq + eta_p*eta_t*z2S*zpq*zqc + eta_p*eta_t*z2S*zpq*zqd + 2*eta_p*z1D*z2S*zpq + eta_p*z1D*zpq*zqc + eta_p*z1D*zpq*zqd + 2*eta_p*z2S*zpa*zpq + eta_p*zpa*zpq*zqc + eta_p*zpa*zpq*zqd - eta_q*z2S2*zpq - eta_q*z2S*zpq*zqc - eta_q*z2S*zpq*zqd - eta_q*zpq*zqc*zqd + eta_t*z2S3 + eta_t*z2S2*zqc + eta_t*z2S2*zqd + eta_t*z2S*zqc*zqd + z1D*z2S2 + z1D*z2S*zqc + z1D*z2S*zqd + z1D*zqc*zqd + z2S2*zpa + z2S*zpa*zqc + z2S*zpa*zqd + zpa*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = inv_ax * 0.5d0 * (  sxpb * ( I_A(abs(0-1)) + I_A(0+1) ) ) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * ( I_dp * cxpb * ( I_A(abs(n-1)) - I_A(n+1) ) + sxpb * ( I_A(abs(n-1)) + I_A(n+1) ) )
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3200) ! | pz  py  s   s    ( 225) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypb*zpq - eta_q*ypq*z1D - eta_q*ypq*zpa + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpa + eta_t*ypb*z2S + y1D*z1D + y1D*zpa + ypb*z1D + ypb*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3201) ! | pz  py  s   px   ( 226) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypb*zpq - eta_q*ypq*z1D - eta_q*ypq*zpa + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpa + eta_t*ypb*z2S + y1D*z1D + y1D*zpa + ypb*z1D + ypb*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3202) ! | pz  py  s   py   ( 227) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpa + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpa + eta_p*eta_t*ypb*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpa + eta_p*ypb*ypq*z1D + eta_p*ypb*ypq*zpa + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqd*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqd*zpq - eta_q*eta_t*ypq*yqd*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqd*zpq - eta_q*y2S*ypb*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpa - eta_q*ypb*yqd*zpq - eta_q*ypq*yqd*z1D - eta_q*ypq*yqd*zpa + eta_t2*y2S2*z2S + eta_t2*y2S*yqd*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqd*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpa + eta_t*y2S*ypb*z2S + eta_t*y2S*yqd*z1D + eta_t*y2S*yqd*zpa + eta_t*ypb*yqd*z2S + y1D*y2S*z1D + y1D*y2S*zpa + y1D*yqd*z1D + y1D*yqd*zpa + y2S*ypb*z1D + y2S*ypb*zpa + ypb*yqd*z1D + ypb*yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3203) ! | pz  py  s   pz   ( 228) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypb*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpa*zpq + eta_p*eta_t*ypb*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpa*zpq + eta_p*ypb*z1D*zpq + eta_p*ypb*zpa*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqd - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqd - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqd - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqd - eta_q*ypb*z2S*zpq - eta_q*ypb*zpq*zqd - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqd - eta_q*ypq*z2S*zpa - eta_q*ypq*zpa*zqd + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqd + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqd + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqd + eta_t*y2S*z2S*zpa + eta_t*y2S*zpa*zqd + eta_t*ypb*z2S2 + eta_t*ypb*z2S*zqd + y1D*z1D*z2S + y1D*z1D*zqd + y1D*z2S*zpa + y1D*zpa*zqd + ypb*z1D*z2S + ypb*z1D*zqd + ypb*z2S*zpa + ypb*zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3210) ! | pz  py  px  s    ( 229) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypb*zpq - eta_q*ypq*z1D - eta_q*ypq*zpa + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpa + eta_t*ypb*z2S + y1D*z1D + y1D*zpa + ypb*z1D + ypb*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3211) ! | pz  py  px  px   ( 230) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*ypq*zpq - eta_q*eta_t*y2S*zpq - eta_q*eta_t*ypq*z2S - eta_q*y1D*zpq - eta_q*ypb*zpq - eta_q*ypq*z1D - eta_q*ypq*zpa + eta_t2*y2S*z2S + eta_t*y1D*z2S + eta_t*y2S*z1D + eta_t*y2S*zpa + eta_t*ypb*z2S + y1D*z1D + y1D*zpa + ypb*z1D + ypb*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3212) ! | pz  py  px  py   ( 231) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpa + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpa + eta_p*eta_t*ypb*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpa + eta_p*ypb*ypq*z1D + eta_p*ypb*ypq*zpa + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqd*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqd*zpq - eta_q*eta_t*ypq*yqd*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqd*zpq - eta_q*y2S*ypb*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpa - eta_q*ypb*yqd*zpq - eta_q*ypq*yqd*z1D - eta_q*ypq*yqd*zpa + eta_t2*y2S2*z2S + eta_t2*y2S*yqd*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqd*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpa + eta_t*y2S*ypb*z2S + eta_t*y2S*yqd*z1D + eta_t*y2S*yqd*zpa + eta_t*ypb*yqd*z2S + y1D*y2S*z1D + y1D*y2S*zpa + y1D*yqd*z1D + y1D*yqd*zpa + y2S*ypb*z1D + y2S*ypb*zpa + ypb*yqd*z1D + ypb*yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3213) ! | pz  py  px  pz   ( 232) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypb*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpa*zpq + eta_p*eta_t*ypb*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpa*zpq + eta_p*ypb*z1D*zpq + eta_p*ypb*zpa*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqd - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqd - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqd - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqd - eta_q*ypb*z2S*zpq - eta_q*ypb*zpq*zqd - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqd - eta_q*ypq*z2S*zpa - eta_q*ypq*zpa*zqd + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqd + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqd + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqd + eta_t*y2S*z2S*zpa + eta_t*y2S*zpa*zqd + eta_t*ypb*z2S2 + eta_t*ypb*z2S*zqd + y1D*z1D*z2S + y1D*z1D*zqd + y1D*z2S*zpa + y1D*zpa*zqd + ypb*z1D*z2S + ypb*z1D*zqd + ypb*z2S*zpa + ypb*zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3220) ! | pz  py  py  s    ( 233) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpa + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpa + eta_p*eta_t*ypb*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpa + eta_p*ypb*ypq*z1D + eta_p*ypb*ypq*zpa + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqc*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqc*zpq - eta_q*eta_t*ypq*yqc*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqc*zpq - eta_q*y2S*ypb*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpa - eta_q*ypb*yqc*zpq - eta_q*ypq*yqc*z1D - eta_q*ypq*yqc*zpa + eta_t2*y2S2*z2S + eta_t2*y2S*yqc*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqc*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpa + eta_t*y2S*ypb*z2S + eta_t*y2S*yqc*z1D + eta_t*y2S*yqc*zpa + eta_t*ypb*yqc*z2S + y1D*y2S*z1D + y1D*y2S*zpa + y1D*yqc*z1D + y1D*yqc*zpa + y2S*ypb*z1D + y2S*ypb*zpa + ypb*yqc*z1D + ypb*yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3221) ! | pz  py  py  px   ( 234) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq2*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq - eta_p*eta_q*eta_t*ypq2*z2S - eta_p*eta_q*y1D*ypq*zpq - eta_p*eta_q*ypb*ypq*zpq - eta_p*eta_q*ypq2*z1D - eta_p*eta_q*ypq2*zpa + eta_p*eta_t2*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*z2S + eta_p*eta_t*y2S*ypq*z1D + eta_p*eta_t*y2S*ypq*zpa + eta_p*eta_t*ypb*ypq*z2S + eta_p*y1D*ypq*z1D + eta_p*y1D*ypq*zpa + eta_p*ypb*ypq*z1D + eta_p*ypb*ypq*zpa + eta_q2*y2S*ypq*zpq + eta_q2*ypq*yqc*zpq - eta_q*eta_t*y2S2*zpq - eta_q*eta_t*y2S*ypq*z2S - eta_q*eta_t*y2S*yqc*zpq - eta_q*eta_t*ypq*yqc*z2S - eta_q*y1D*y2S*zpq - eta_q*y1D*yqc*zpq - eta_q*y2S*ypb*zpq - eta_q*y2S*ypq*z1D - eta_q*y2S*ypq*zpa - eta_q*ypb*yqc*zpq - eta_q*ypq*yqc*z1D - eta_q*ypq*yqc*zpa + eta_t2*y2S2*z2S + eta_t2*y2S*yqc*z2S + eta_t*y1D*y2S*z2S + eta_t*y1D*yqc*z2S + eta_t*y2S2*z1D + eta_t*y2S2*zpa + eta_t*y2S*ypb*z2S + eta_t*y2S*yqc*z1D + eta_t*y2S*yqc*zpa + eta_t*ypb*yqc*z2S + y1D*y2S*z1D + y1D*y2S*zpa + y1D*yqc*z1D + y1D*yqc*zpa + y2S*ypb*z1D + y2S*ypb*zpa + ypb*yqc*z1D + ypb*yqc*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3222) ! | pz  py  py  py   ( 235) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq3*zpq - eta_p2*eta_q*eta_t*y2S*ypq2*zpq - eta_p2*eta_q*eta_t*ypq3*z2S - eta_p2*eta_q*y1D*ypq2*zpq - eta_p2*eta_q*ypb*ypq2*zpq - eta_p2*eta_q*ypq3*z1D - eta_p2*eta_q*ypq3*zpa + eta_p2*eta_t2*y2S*ypq2*z2S + eta_p2*eta_t*y1D*ypq2*z2S + eta_p2*eta_t*y2S*ypq2*z1D + eta_p2*eta_t*y2S*ypq2*zpa + eta_p2*eta_t*ypb*ypq2*z2S + eta_p2*y1D*ypq2*z1D + eta_p2*y1D*ypq2*zpa + eta_p2*ypb*ypq2*z1D + eta_p2*ypb*ypq2*zpa + 2*eta_p*eta_q2*y2S*ypq2*zpq + eta_p*eta_q2*ypq2*yqc*zpq + eta_p*eta_q2*ypq2*yqd*zpq - 2*eta_p*eta_q*eta_t*y2S2*ypq*zpq - 2*eta_p*eta_q*eta_t*y2S*ypq2*z2S - eta_p*eta_q*eta_t*y2S*ypq*yqc*zpq - eta_p*eta_q*eta_t*y2S*ypq*yqd*zpq - eta_p*eta_q*eta_t*ypq2*yqc*z2S - eta_p*eta_q*eta_t*ypq2*yqd*z2S - 2*eta_p*eta_q*y1D*y2S*ypq*zpq - eta_p*eta_q*y1D*ypq*yqc*zpq - eta_p*eta_q*y1D*ypq*yqd*zpq - 2*eta_p*eta_q*y2S*ypb*ypq*zpq - 2*eta_p*eta_q*y2S*ypq2*z1D - 2*eta_p*eta_q*y2S*ypq2*zpa - eta_p*eta_q*ypb*ypq*yqc*zpq - eta_p*eta_q*ypb*ypq*yqd*zpq - eta_p*eta_q*ypq2*yqc*z1D - eta_p*eta_q*ypq2*yqc*zpa - eta_p*eta_q*ypq2*yqd*z1D - eta_p*eta_q*ypq2*yqd*zpa + 2*eta_p*eta_t2*y2S2*ypq*z2S + eta_p*eta_t2*y2S*ypq*yqc*z2S + eta_p*eta_t2*y2S*ypq*yqd*z2S + 2*eta_p*eta_t*y1D*y2S*ypq*z2S + eta_p*eta_t*y1D*ypq*yqc*z2S + eta_p*eta_t*y1D*ypq*yqd*z2S + 2*eta_p*eta_t*y2S2*ypq*z1D + 2*eta_p*eta_t*y2S2*ypq*zpa + 2*eta_p*eta_t*y2S*ypb*ypq*z2S + eta_p*eta_t*y2S*ypq*yqc*z1D + eta_p*eta_t*y2S*ypq*yqc*zpa + eta_p*eta_t*y2S*ypq*yqd*z1D + eta_p*eta_t*y2S*ypq*yqd*zpa + eta_p*eta_t*ypb*ypq*yqc*z2S + eta_p*eta_t*ypb*ypq*yqd*z2S + 2*eta_p*y1D*y2S*ypq*z1D + 2*eta_p*y1D*y2S*ypq*zpa + eta_p*y1D*ypq*yqc*z1D + eta_p*y1D*ypq*yqc*zpa + eta_p*y1D*ypq*yqd*z1D + eta_p*y1D*ypq*yqd*zpa + 2*eta_p*y2S*ypb*ypq*z1D + 2*eta_p*y2S*ypb*ypq*zpa + eta_p*ypb*ypq*yqc*z1D + eta_p*ypb*ypq*yqc*zpa + eta_p*ypb*ypq*yqd*z1D + eta_p*ypb*ypq*yqd*zpa + eta_q2*y2S2*ypq*zpq + eta_q2*y2S*ypq*yqc*zpq + eta_q2*y2S*ypq*yqd*zpq + eta_q2*ypq*yqc*yqd*zpq - eta_q*eta_t*y2S3*zpq - eta_q*eta_t*y2S2*ypq*z2S - eta_q*eta_t*y2S2*yqc*zpq - eta_q*eta_t*y2S2*yqd*zpq - eta_q*eta_t*y2S*ypq*yqc*z2S - eta_q*eta_t*y2S*ypq*yqd*z2S - eta_q*eta_t*y2S*yqc*yqd*zpq - eta_q*eta_t*ypq*yqc*yqd*z2S - eta_q*y1D*y2S2*zpq - eta_q*y1D*y2S*yqc*zpq - eta_q*y1D*y2S*yqd*zpq - eta_q*y1D*yqc*yqd*zpq - eta_q*y2S2*ypb*zpq - eta_q*y2S2*ypq*z1D - eta_q*y2S2*ypq*zpa - eta_q*y2S*ypb*yqc*zpq - eta_q*y2S*ypb*yqd*zpq - eta_q*y2S*ypq*yqc*z1D - eta_q*y2S*ypq*yqc*zpa - eta_q*y2S*ypq*yqd*z1D - eta_q*y2S*ypq*yqd*zpa - eta_q*ypb*yqc*yqd*zpq - eta_q*ypq*yqc*yqd*z1D - eta_q*ypq*yqc*yqd*zpa + eta_t2*y2S3*z2S + eta_t2*y2S2*yqc*z2S + eta_t2*y2S2*yqd*z2S + eta_t2*y2S*yqc*yqd*z2S + eta_t*y1D*y2S2*z2S + eta_t*y1D*y2S*yqc*z2S + eta_t*y1D*y2S*yqd*z2S + eta_t*y1D*yqc*yqd*z2S + eta_t*y2S3*z1D + eta_t*y2S3*zpa + eta_t*y2S2*ypb*z2S + eta_t*y2S2*yqc*z1D + eta_t*y2S2*yqc*zpa + eta_t*y2S2*yqd*z1D + eta_t*y2S2*yqd*zpa + eta_t*y2S*ypb*yqc*z2S + eta_t*y2S*ypb*yqd*z2S + eta_t*y2S*yqc*yqd*z1D + eta_t*y2S*yqc*yqd*zpa + eta_t*ypb*yqc*yqd*z2S + y1D*y2S2*z1D + y1D*y2S2*zpa + y1D*y2S*yqc*z1D + y1D*y2S*yqc*zpa + y1D*y2S*yqd*z1D + y1D*y2S*yqd*zpa + y1D*yqc*yqd*z1D + y1D*yqc*yqd*zpa + y2S2*ypb*z1D + y2S2*ypb*zpa + y2S*ypb*yqc*z1D + y2S*ypb*yqc*zpa + y2S*ypb*yqd*z1D + y2S*ypb*yqd*zpa + ypb*yqc*yqd*z1D + ypb*yqc*yqd*zpa) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3223) ! | pz  py  py  pz   ( 236) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - eta_p2*eta_q*eta_t*y2S*ypq*zpq2 - eta_p2*eta_q*eta_t*ypq2*z2S*zpq - eta_p2*eta_q*y1D*ypq*zpq2 - eta_p2*eta_q*ypb*ypq*zpq2 - eta_p2*eta_q*ypq2*z1D*zpq - eta_p2*eta_q*ypq2*zpa*zpq + eta_p2*eta_t2*y2S*ypq*z2S*zpq + eta_p2*eta_t*y1D*ypq*z2S*zpq + eta_p2*eta_t*y2S*ypq*z1D*zpq + eta_p2*eta_t*y2S*ypq*zpa*zpq + eta_p2*eta_t*ypb*ypq*z2S*zpq + eta_p2*y1D*ypq*z1D*zpq + eta_p2*y1D*ypq*zpa*zpq + eta_p2*ypb*ypq*z1D*zpq + eta_p2*ypb*ypq*zpa*zpq + eta_p*eta_q2*y2S*ypq*zpq2 + eta_p*eta_q2*ypq2*z2S*zpq + eta_p*eta_q2*ypq2*zpq*zqd + eta_p*eta_q2*ypq*yqc*zpq2 - eta_p*eta_q*eta_t*y2S2*zpq2 - 2*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq*zqd - eta_p*eta_q*eta_t*y2S*yqc*zpq2 - eta_p*eta_q*eta_t*ypq2*z2S2 - eta_p*eta_q*eta_t*ypq2*z2S*zqd - eta_p*eta_q*eta_t*ypq*yqc*z2S*zpq - eta_p*eta_q*y1D*y2S*zpq2 - eta_p*eta_q*y1D*ypq*z2S*zpq - eta_p*eta_q*y1D*ypq*zpq*zqd - eta_p*eta_q*y1D*yqc*zpq2 - eta_p*eta_q*y2S*ypb*zpq2 - eta_p*eta_q*y2S*ypq*z1D*zpq - eta_p*eta_q*y2S*ypq*zpa*zpq - eta_p*eta_q*ypb*ypq*z2S*zpq - eta_p*eta_q*ypb*ypq*zpq*zqd - eta_p*eta_q*ypb*yqc*zpq2 - eta_p*eta_q*ypq2*z1D*z2S - eta_p*eta_q*ypq2*z1D*zqd - eta_p*eta_q*ypq2*z2S*zpa - eta_p*eta_q*ypq2*zpa*zqd - eta_p*eta_q*ypq*yqc*z1D*zpq - eta_p*eta_q*ypq*yqc*zpa*zpq + eta_p*eta_t2*y2S2*z2S*zpq + eta_p*eta_t2*y2S*ypq*z2S2 + eta_p*eta_t2*y2S*ypq*z2S*zqd + eta_p*eta_t2*y2S*yqc*z2S*zpq + eta_p*eta_t*y1D*y2S*z2S*zpq + eta_p*eta_t*y1D*ypq*z2S2 + eta_p*eta_t*y1D*ypq*z2S*zqd + eta_p*eta_t*y1D*yqc*z2S*zpq + eta_p*eta_t*y2S2*z1D*zpq + eta_p*eta_t*y2S2*zpa*zpq + eta_p*eta_t*y2S*ypb*z2S*zpq + eta_p*eta_t*y2S*ypq*z1D*z2S + eta_p*eta_t*y2S*ypq*z1D*zqd + eta_p*eta_t*y2S*ypq*z2S*zpa + eta_p*eta_t*y2S*ypq*zpa*zqd + eta_p*eta_t*y2S*yqc*z1D*zpq + eta_p*eta_t*y2S*yqc*zpa*zpq + eta_p*eta_t*ypb*ypq*z2S2 + eta_p*eta_t*ypb*ypq*z2S*zqd + eta_p*eta_t*ypb*yqc*z2S*zpq + eta_p*y1D*y2S*z1D*zpq + eta_p*y1D*y2S*zpa*zpq + eta_p*y1D*ypq*z1D*z2S + eta_p*y1D*ypq*z1D*zqd + eta_p*y1D*ypq*z2S*zpa + eta_p*y1D*ypq*zpa*zqd + eta_p*y1D*yqc*z1D*zpq + eta_p*y1D*yqc*zpa*zpq + eta_p*y2S*ypb*z1D*zpq + eta_p*y2S*ypb*zpa*zpq + eta_p*ypb*ypq*z1D*z2S + eta_p*ypb*ypq*z1D*zqd + eta_p*ypb*ypq*z2S*zpa + eta_p*ypb*ypq*zpa*zqd + eta_p*ypb*yqc*z1D*zpq + eta_p*ypb*yqc*zpa*zpq + eta_q2*y2S*ypq*z2S*zpq + eta_q2*y2S*ypq*zpq*zqd + eta_q2*ypq*yqc*z2S*zpq + eta_q2*ypq*yqc*zpq*zqd - eta_q*eta_t*y2S2*z2S*zpq - eta_q*eta_t*y2S2*zpq*zqd - eta_q*eta_t*y2S*ypq*z2S2 - eta_q*eta_t*y2S*ypq*z2S*zqd - eta_q*eta_t*y2S*yqc*z2S*zpq - eta_q*eta_t*y2S*yqc*zpq*zqd - eta_q*eta_t*ypq*yqc*z2S2 - eta_q*eta_t*ypq*yqc*z2S*zqd - eta_q*y1D*y2S*z2S*zpq - eta_q*y1D*y2S*zpq*zqd - eta_q*y1D*yqc*z2S*zpq - eta_q*y1D*yqc*zpq*zqd - eta_q*y2S*ypb*z2S*zpq - eta_q*y2S*ypb*zpq*zqd - eta_q*y2S*ypq*z1D*z2S - eta_q*y2S*ypq*z1D*zqd - eta_q*y2S*ypq*z2S*zpa - eta_q*y2S*ypq*zpa*zqd - eta_q*ypb*yqc*z2S*zpq - eta_q*ypb*yqc*zpq*zqd - eta_q*ypq*yqc*z1D*z2S - eta_q*ypq*yqc*z1D*zqd - eta_q*ypq*yqc*z2S*zpa - eta_q*ypq*yqc*zpa*zqd + eta_t2*y2S2*z2S2 + eta_t2*y2S2*z2S*zqd + eta_t2*y2S*yqc*z2S2 + eta_t2*y2S*yqc*z2S*zqd + eta_t*y1D*y2S*z2S2 + eta_t*y1D*y2S*z2S*zqd + eta_t*y1D*yqc*z2S2 + eta_t*y1D*yqc*z2S*zqd + eta_t*y2S2*z1D*z2S + eta_t*y2S2*z1D*zqd + eta_t*y2S2*z2S*zpa + eta_t*y2S2*zpa*zqd + eta_t*y2S*ypb*z2S2 + eta_t*y2S*ypb*z2S*zqd + eta_t*y2S*yqc*z1D*z2S + eta_t*y2S*yqc*z1D*zqd + eta_t*y2S*yqc*z2S*zpa + eta_t*y2S*yqc*zpa*zqd + eta_t*ypb*yqc*z2S2 + eta_t*ypb*yqc*z2S*zqd + y1D*y2S*z1D*z2S + y1D*y2S*z1D*zqd + y1D*y2S*z2S*zpa + y1D*y2S*zpa*zqd + y1D*yqc*z1D*z2S + y1D*yqc*z1D*zqd + y1D*yqc*z2S*zpa + y1D*yqc*zpa*zqd + y2S*ypb*z1D*z2S + y2S*ypb*z1D*zqd + y2S*ypb*z2S*zpa + y2S*ypb*zpa*zqd + ypb*yqc*z1D*z2S + ypb*yqc*z1D*zqd + ypb*yqc*z2S*zpa + ypb*yqc*zpa*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3230) ! | pz  py  pz  s    ( 237) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypb*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpa*zpq + eta_p*eta_t*ypb*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpa*zpq + eta_p*ypb*z1D*zpq + eta_p*ypb*zpa*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqc - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqc - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqc - eta_q*ypb*z2S*zpq - eta_q*ypb*zpq*zqc - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqc - eta_q*ypq*z2S*zpa - eta_q*ypq*zpa*zqc + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqc + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqc + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqc + eta_t*y2S*z2S*zpa + eta_t*y2S*zpa*zqc + eta_t*ypb*z2S2 + eta_t*ypb*z2S*zqc + y1D*z1D*z2S + y1D*z1D*zqc + y1D*z2S*zpa + y1D*zpa*zqc + ypb*z1D*z2S + ypb*z1D*zqc + ypb*z2S*zpa + ypb*zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3231) ! | pz  py  pz  px   ( 238) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2 - eta_p*eta_q*eta_t*ypq*z2S*zpq - eta_p*eta_q*y1D*zpq2 - eta_p*eta_q*ypb*zpq2 - eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq + eta_p*eta_t2*y2S*z2S*zpq + eta_p*eta_t*y1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq + eta_p*eta_t*y2S*zpa*zpq + eta_p*eta_t*ypb*z2S*zpq + eta_p*y1D*z1D*zpq + eta_p*y1D*zpa*zpq + eta_p*ypb*z1D*zpq + eta_p*ypb*zpa*zpq + eta_q2*ypq*z2S*zpq + eta_q2*ypq*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq - eta_q*eta_t*y2S*zpq*zqc - eta_q*eta_t*ypq*z2S2 - eta_q*eta_t*ypq*z2S*zqc - eta_q*y1D*z2S*zpq - eta_q*y1D*zpq*zqc - eta_q*ypb*z2S*zpq - eta_q*ypb*zpq*zqc - eta_q*ypq*z1D*z2S - eta_q*ypq*z1D*zqc - eta_q*ypq*z2S*zpa - eta_q*ypq*zpa*zqc + eta_t2*y2S*z2S2 + eta_t2*y2S*z2S*zqc + eta_t*y1D*z2S2 + eta_t*y1D*z2S*zqc + eta_t*y2S*z1D*z2S + eta_t*y2S*z1D*zqc + eta_t*y2S*z2S*zpa + eta_t*y2S*zpa*zqc + eta_t*ypb*z2S2 + eta_t*ypb*z2S*zqc + y1D*z1D*z2S + y1D*z1D*zqc + y1D*z2S*zpa + y1D*zpa*zqc + ypb*z1D*z2S + ypb*z1D*zqc + ypb*z2S*zpa + ypb*zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3232) ! | pz  py  pz  py   ( 239) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - eta_p2*eta_q*eta_t*y2S*ypq*zpq2 - eta_p2*eta_q*eta_t*ypq2*z2S*zpq - eta_p2*eta_q*y1D*ypq*zpq2 - eta_p2*eta_q*ypb*ypq*zpq2 - eta_p2*eta_q*ypq2*z1D*zpq - eta_p2*eta_q*ypq2*zpa*zpq + eta_p2*eta_t2*y2S*ypq*z2S*zpq + eta_p2*eta_t*y1D*ypq*z2S*zpq + eta_p2*eta_t*y2S*ypq*z1D*zpq + eta_p2*eta_t*y2S*ypq*zpa*zpq + eta_p2*eta_t*ypb*ypq*z2S*zpq + eta_p2*y1D*ypq*z1D*zpq + eta_p2*y1D*ypq*zpa*zpq + eta_p2*ypb*ypq*z1D*zpq + eta_p2*ypb*ypq*zpa*zpq + eta_p*eta_q2*y2S*ypq*zpq2 + eta_p*eta_q2*ypq2*z2S*zpq + eta_p*eta_q2*ypq2*zpq*zqc + eta_p*eta_q2*ypq*yqd*zpq2 - eta_p*eta_q*eta_t*y2S2*zpq2 - 2*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - eta_p*eta_q*eta_t*y2S*ypq*zpq*zqc - eta_p*eta_q*eta_t*y2S*yqd*zpq2 - eta_p*eta_q*eta_t*ypq2*z2S2 - eta_p*eta_q*eta_t*ypq2*z2S*zqc - eta_p*eta_q*eta_t*ypq*yqd*z2S*zpq - eta_p*eta_q*y1D*y2S*zpq2 - eta_p*eta_q*y1D*ypq*z2S*zpq - eta_p*eta_q*y1D*ypq*zpq*zqc - eta_p*eta_q*y1D*yqd*zpq2 - eta_p*eta_q*y2S*ypb*zpq2 - eta_p*eta_q*y2S*ypq*z1D*zpq - eta_p*eta_q*y2S*ypq*zpa*zpq - eta_p*eta_q*ypb*ypq*z2S*zpq - eta_p*eta_q*ypb*ypq*zpq*zqc - eta_p*eta_q*ypb*yqd*zpq2 - eta_p*eta_q*ypq2*z1D*z2S - eta_p*eta_q*ypq2*z1D*zqc - eta_p*eta_q*ypq2*z2S*zpa - eta_p*eta_q*ypq2*zpa*zqc - eta_p*eta_q*ypq*yqd*z1D*zpq - eta_p*eta_q*ypq*yqd*zpa*zpq + eta_p*eta_t2*y2S2*z2S*zpq + eta_p*eta_t2*y2S*ypq*z2S2 + eta_p*eta_t2*y2S*ypq*z2S*zqc + eta_p*eta_t2*y2S*yqd*z2S*zpq + eta_p*eta_t*y1D*y2S*z2S*zpq + eta_p*eta_t*y1D*ypq*z2S2 + eta_p*eta_t*y1D*ypq*z2S*zqc + eta_p*eta_t*y1D*yqd*z2S*zpq + eta_p*eta_t*y2S2*z1D*zpq + eta_p*eta_t*y2S2*zpa*zpq + eta_p*eta_t*y2S*ypb*z2S*zpq + eta_p*eta_t*y2S*ypq*z1D*z2S + eta_p*eta_t*y2S*ypq*z1D*zqc + eta_p*eta_t*y2S*ypq*z2S*zpa + eta_p*eta_t*y2S*ypq*zpa*zqc + eta_p*eta_t*y2S*yqd*z1D*zpq + eta_p*eta_t*y2S*yqd*zpa*zpq + eta_p*eta_t*ypb*ypq*z2S2 + eta_p*eta_t*ypb*ypq*z2S*zqc + eta_p*eta_t*ypb*yqd*z2S*zpq + eta_p*y1D*y2S*z1D*zpq + eta_p*y1D*y2S*zpa*zpq + eta_p*y1D*ypq*z1D*z2S + eta_p*y1D*ypq*z1D*zqc + eta_p*y1D*ypq*z2S*zpa + eta_p*y1D*ypq*zpa*zqc + eta_p*y1D*yqd*z1D*zpq + eta_p*y1D*yqd*zpa*zpq + eta_p*y2S*ypb*z1D*zpq + eta_p*y2S*ypb*zpa*zpq + eta_p*ypb*ypq*z1D*z2S + eta_p*ypb*ypq*z1D*zqc + eta_p*ypb*ypq*z2S*zpa + eta_p*ypb*ypq*zpa*zqc + eta_p*ypb*yqd*z1D*zpq + eta_p*ypb*yqd*zpa*zpq + eta_q2*y2S*ypq*z2S*zpq + eta_q2*y2S*ypq*zpq*zqc + eta_q2*ypq*yqd*z2S*zpq + eta_q2*ypq*yqd*zpq*zqc - eta_q*eta_t*y2S2*z2S*zpq - eta_q*eta_t*y2S2*zpq*zqc - eta_q*eta_t*y2S*ypq*z2S2 - eta_q*eta_t*y2S*ypq*z2S*zqc - eta_q*eta_t*y2S*yqd*z2S*zpq - eta_q*eta_t*y2S*yqd*zpq*zqc - eta_q*eta_t*ypq*yqd*z2S2 - eta_q*eta_t*ypq*yqd*z2S*zqc - eta_q*y1D*y2S*z2S*zpq - eta_q*y1D*y2S*zpq*zqc - eta_q*y1D*yqd*z2S*zpq - eta_q*y1D*yqd*zpq*zqc - eta_q*y2S*ypb*z2S*zpq - eta_q*y2S*ypb*zpq*zqc - eta_q*y2S*ypq*z1D*z2S - eta_q*y2S*ypq*z1D*zqc - eta_q*y2S*ypq*z2S*zpa - eta_q*y2S*ypq*zpa*zqc - eta_q*ypb*yqd*z2S*zpq - eta_q*ypb*yqd*zpq*zqc - eta_q*ypq*yqd*z1D*z2S - eta_q*ypq*yqd*z1D*zqc - eta_q*ypq*yqd*z2S*zpa - eta_q*ypq*yqd*zpa*zqc + eta_t2*y2S2*z2S2 + eta_t2*y2S2*z2S*zqc + eta_t2*y2S*yqd*z2S2 + eta_t2*y2S*yqd*z2S*zqc + eta_t*y1D*y2S*z2S2 + eta_t*y1D*y2S*z2S*zqc + eta_t*y1D*yqd*z2S2 + eta_t*y1D*yqd*z2S*zqc + eta_t*y2S2*z1D*z2S + eta_t*y2S2*z1D*zqc + eta_t*y2S2*z2S*zpa + eta_t*y2S2*zpa*zqc + eta_t*y2S*ypb*z2S2 + eta_t*y2S*ypb*z2S*zqc + eta_t*y2S*yqd*z1D*z2S + eta_t*y2S*yqd*z1D*zqc + eta_t*y2S*yqd*z2S*zpa + eta_t*y2S*yqd*zpa*zqc + eta_t*ypb*yqd*z2S2 + eta_t*ypb*yqd*z2S*zqc + y1D*y2S*z1D*z2S + y1D*y2S*z1D*zqc + y1D*y2S*z2S*zpa + y1D*y2S*zpa*zqc + y1D*yqd*z1D*z2S + y1D*yqd*z1D*zqc + y1D*yqd*z2S*zpa + y1D*yqd*zpa*zqc + y2S*ypb*z1D*z2S + y2S*ypb*z1D*zqc + y2S*ypb*z2S*zpa + y2S*ypb*zpa*zqc + ypb*yqd*z1D*z2S + ypb*yqd*z1D*zqc + ypb*yqd*z2S*zpa + ypb*yqd*zpa*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3233) ! | pz  py  pz  pz   ( 240) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq*zpq3 - eta_p2*eta_q*eta_t*y2S*zpq3 - eta_p2*eta_q*eta_t*ypq*z2S*zpq2 - eta_p2*eta_q*y1D*zpq3 - eta_p2*eta_q*ypb*zpq3 - eta_p2*eta_q*ypq*z1D*zpq2 - eta_p2*eta_q*ypq*zpa*zpq2 + eta_p2*eta_t2*y2S*z2S*zpq2 + eta_p2*eta_t*y1D*z2S*zpq2 + eta_p2*eta_t*y2S*z1D*zpq2 + eta_p2*eta_t*y2S*zpa*zpq2 + eta_p2*eta_t*ypb*z2S*zpq2 + eta_p2*y1D*z1D*zpq2 + eta_p2*y1D*zpa*zpq2 + eta_p2*ypb*z1D*zpq2 + eta_p2*ypb*zpa*zpq2 + 2*eta_p*eta_q2*ypq*z2S*zpq2 + eta_p*eta_q2*ypq*zpq2*zqc + eta_p*eta_q2*ypq*zpq2*zqd - 2*eta_p*eta_q*eta_t*y2S*z2S*zpq2 - eta_p*eta_q*eta_t*y2S*zpq2*zqc - eta_p*eta_q*eta_t*y2S*zpq2*zqd - 2*eta_p*eta_q*eta_t*ypq*z2S2*zpq - eta_p*eta_q*eta_t*ypq*z2S*zpq*zqc - eta_p*eta_q*eta_t*ypq*z2S*zpq*zqd - 2*eta_p*eta_q*y1D*z2S*zpq2 - eta_p*eta_q*y1D*zpq2*zqc - eta_p*eta_q*y1D*zpq2*zqd - 2*eta_p*eta_q*ypb*z2S*zpq2 - eta_p*eta_q*ypb*zpq2*zqc - eta_p*eta_q*ypb*zpq2*zqd - 2*eta_p*eta_q*ypq*z1D*z2S*zpq - eta_p*eta_q*ypq*z1D*zpq*zqc - eta_p*eta_q*ypq*z1D*zpq*zqd - 2*eta_p*eta_q*ypq*z2S*zpa*zpq - eta_p*eta_q*ypq*zpa*zpq*zqc - eta_p*eta_q*ypq*zpa*zpq*zqd + 2*eta_p*eta_t2*y2S*z2S2*zpq + eta_p*eta_t2*y2S*z2S*zpq*zqc + eta_p*eta_t2*y2S*z2S*zpq*zqd + 2*eta_p*eta_t*y1D*z2S2*zpq + eta_p*eta_t*y1D*z2S*zpq*zqc + eta_p*eta_t*y1D*z2S*zpq*zqd + 2*eta_p*eta_t*y2S*z1D*z2S*zpq + eta_p*eta_t*y2S*z1D*zpq*zqc + eta_p*eta_t*y2S*z1D*zpq*zqd + 2*eta_p*eta_t*y2S*z2S*zpa*zpq + eta_p*eta_t*y2S*zpa*zpq*zqc + eta_p*eta_t*y2S*zpa*zpq*zqd + 2*eta_p*eta_t*ypb*z2S2*zpq + eta_p*eta_t*ypb*z2S*zpq*zqc + eta_p*eta_t*ypb*z2S*zpq*zqd + 2*eta_p*y1D*z1D*z2S*zpq + eta_p*y1D*z1D*zpq*zqc + eta_p*y1D*z1D*zpq*zqd + 2*eta_p*y1D*z2S*zpa*zpq + eta_p*y1D*zpa*zpq*zqc + eta_p*y1D*zpa*zpq*zqd + 2*eta_p*ypb*z1D*z2S*zpq + eta_p*ypb*z1D*zpq*zqc + eta_p*ypb*z1D*zpq*zqd + 2*eta_p*ypb*z2S*zpa*zpq + eta_p*ypb*zpa*zpq*zqc + eta_p*ypb*zpa*zpq*zqd + eta_q2*ypq*z2S2*zpq + eta_q2*ypq*z2S*zpq*zqc + eta_q2*ypq*z2S*zpq*zqd + eta_q2*ypq*zpq*zqc*zqd - eta_q*eta_t*y2S*z2S2*zpq - eta_q*eta_t*y2S*z2S*zpq*zqc - eta_q*eta_t*y2S*z2S*zpq*zqd - eta_q*eta_t*y2S*zpq*zqc*zqd - eta_q*eta_t*ypq*z2S3 - eta_q*eta_t*ypq*z2S2*zqc - eta_q*eta_t*ypq*z2S2*zqd - eta_q*eta_t*ypq*z2S*zqc*zqd - eta_q*y1D*z2S2*zpq - eta_q*y1D*z2S*zpq*zqc - eta_q*y1D*z2S*zpq*zqd - eta_q*y1D*zpq*zqc*zqd - eta_q*ypb*z2S2*zpq - eta_q*ypb*z2S*zpq*zqc - eta_q*ypb*z2S*zpq*zqd - eta_q*ypb*zpq*zqc*zqd - eta_q*ypq*z1D*z2S2 - eta_q*ypq*z1D*z2S*zqc - eta_q*ypq*z1D*z2S*zqd - eta_q*ypq*z1D*zqc*zqd - eta_q*ypq*z2S2*zpa - eta_q*ypq*z2S*zpa*zqc - eta_q*ypq*z2S*zpa*zqd - eta_q*ypq*zpa*zqc*zqd + eta_t2*y2S*z2S3 + eta_t2*y2S*z2S2*zqc + eta_t2*y2S*z2S2*zqd + eta_t2*y2S*z2S*zqc*zqd + eta_t*y1D*z2S3 + eta_t*y1D*z2S2*zqc + eta_t*y1D*z2S2*zqd + eta_t*y1D*z2S*zqc*zqd + eta_t*y2S*z1D*z2S2 + eta_t*y2S*z1D*z2S*zqc + eta_t*y2S*z1D*z2S*zqd + eta_t*y2S*z1D*zqc*zqd + eta_t*y2S*z2S2*zpa + eta_t*y2S*z2S*zpa*zqc + eta_t*y2S*z2S*zpa*zqd + eta_t*y2S*zpa*zqc*zqd + eta_t*ypb*z2S3 + eta_t*ypb*z2S2*zqc + eta_t*ypb*z2S2*zqd + eta_t*ypb*z2S*zqc*zqd + y1D*z1D*z2S2 + y1D*z1D*z2S*zqc + y1D*z1D*z2S*zqd + y1D*z1D*zqc*zqd + y1D*z2S2*zpa + y1D*z2S*zpa*zqc + y1D*z2S*zpa*zqd + y1D*zpa*zqc*zqd + ypb*z1D*z2S2 + ypb*z1D*z2S*zqc + ypb*z1D*z2S*zqd + ypb*z1D*zqc*zqd + ypb*z2S2*zpa + ypb*z2S*zpa*zqc + ypb*z2S*zpa*zqd + ypb*zpa*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3300) ! | pz  pz  s   s    ( 241) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*zpq2 - 2*eta_q*eta_t*z2S*zpq - 2*eta_q*z1D*zpq - eta_q*zpa*zpq - eta_q*zpb*zpq + eta_t2*z2S2 + 2*eta_t*z1D*z2S + eta_t*z2S*zpa + eta_t*z2S*zpb + z1D2 + z1D*zpa + z1D*zpb + zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3301) ! | pz  pz  s   px   ( 242) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*zpq2 - 2*eta_q*eta_t*z2S*zpq - 2*eta_q*z1D*zpq - eta_q*zpa*zpq - eta_q*zpb*zpq + eta_t2*z2S2 + 2*eta_t*z1D*z2S + eta_t*z2S*zpa + eta_t*z2S*zpb + z1D2 + z1D*zpa + z1D*zpb + zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3302) ! | pz  pz  s   py   ( 243) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*ypq*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S + eta_p*eta_t*ypq*z2S*zpa + eta_p*eta_t*ypq*z2S*zpb + eta_p*ypq*z1D2 + eta_p*ypq*z1D*zpa + eta_p*ypq*z1D*zpb + eta_p*ypq*zpa*zpb + eta_q2*y2S*zpq2 + eta_q2*yqd*zpq2 - 2*eta_q*eta_t*y2S*z2S*zpq - 2*eta_q*eta_t*yqd*z2S*zpq - 2*eta_q*y2S*z1D*zpq - eta_q*y2S*zpa*zpq - eta_q*y2S*zpb*zpq - 2*eta_q*yqd*z1D*zpq - eta_q*yqd*zpa*zpq - eta_q*yqd*zpb*zpq + eta_t2*y2S*z2S2 + eta_t2*yqd*z2S2 + 2*eta_t*y2S*z1D*z2S + eta_t*y2S*z2S*zpa + eta_t*y2S*z2S*zpb + 2*eta_t*yqd*z1D*z2S + eta_t*yqd*z2S*zpa + eta_t*yqd*z2S*zpb + y2S*z1D2 + y2S*z1D*zpa + y2S*z1D*zpb + y2S*zpa*zpb + yqd*z1D2 + yqd*z1D*zpa + yqd*z1D*zpb + yqd*zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3303) ! | pz  pz  s   pz   ( 244) 

      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*zpq3 - 2*eta_p*eta_q*eta_t*z2S*zpq2 - 2*eta_p*eta_q*z1D*zpq2 - eta_p*eta_q*zpa*zpq2 - eta_p*eta_q*zpb*zpq2 + eta_p*eta_t2*z2S2*zpq + 2*eta_p*eta_t*z1D*z2S*zpq + eta_p*eta_t*z2S*zpa*zpq + eta_p*eta_t*z2S*zpb*zpq + eta_p*z1D2*zpq + eta_p*z1D*zpa*zpq + eta_p*z1D*zpb*zpq + eta_p*zpa*zpb*zpq + eta_q2*z2S*zpq2 + eta_q2*zpq2*zqd - 2*eta_q*eta_t*z2S2*zpq - 2*eta_q*eta_t*z2S*zpq*zqd - 2*eta_q*z1D*z2S*zpq - 2*eta_q*z1D*zpq*zqd - eta_q*z2S*zpa*zpq - eta_q*z2S*zpb*zpq - eta_q*zpa*zpq*zqd - eta_q*zpb*zpq*zqd + eta_t2*z2S3 + eta_t2*z2S2*zqd + 2*eta_t*z1D*z2S2 + 2*eta_t*z1D*z2S*zqd + eta_t*z2S2*zpa + eta_t*z2S2*zpb + eta_t*z2S*zpa*zqd + eta_t*z2S*zpb*zqd + z1D2*z2S + z1D2*zqd + z1D*z2S*zpa + z1D*z2S*zpb + z1D*zpa*zqd + z1D*zpb*zqd + z2S*zpa*zpb + zpa*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3310) ! | pz  pz  px  s    ( 245) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*zpq2 - 2*eta_q*eta_t*z2S*zpq - 2*eta_q*z1D*zpq - eta_q*zpa*zpq - eta_q*zpb*zpq + eta_t2*z2S2 + 2*eta_t*z1D*z2S + eta_t*z2S*zpa + eta_t*z2S*zpb + z1D2 + z1D*zpa + z1D*zpb + zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3311) ! | pz  pz  px  px   ( 246) 
      
      n         = 0
      const     = pi2 * D2 * (eta_q2*zpq2 - 2*eta_q*eta_t*z2S*zpq - 2*eta_q*z1D*zpq - eta_q*zpa*zpq - eta_q*zpb*zpq + eta_t2*z2S2 + 2*eta_t*z1D*z2S + eta_t*z2S*zpa + eta_t*z2S*zpb + z1D2 + z1D*zpa + z1D*zpb + zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) )
      sum       = I_A(0) * inv_ax2 * (cxqc * cxqd * I_B(0) - c2xqcd * (0.25d0 * (I_B(abs(0-2)) + 2.d0 * I_B(0) + I_B(abs(0+2)) ) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        if (dabs(B) < 1.d-300) then
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) )
        else
          termBn  = inv_ax2 * (cxqc * cxqd * I_B(n) - c2xqcd * (0.25d0 * (I_B(abs(n-2)) + 2.d0 * I_B(n) + I_B(abs(n+2)) ) ) - I_dp/B * n * s2xqcd * (0.5d0 * ( I_B(abs(n-1)) + I_B(abs(n+1)) ) ) )
        end if
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3312) ! | pz  pz  px  py   ( 247) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*ypq*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S + eta_p*eta_t*ypq*z2S*zpa + eta_p*eta_t*ypq*z2S*zpb + eta_p*ypq*z1D2 + eta_p*ypq*z1D*zpa + eta_p*ypq*z1D*zpb + eta_p*ypq*zpa*zpb + eta_q2*y2S*zpq2 + eta_q2*yqd*zpq2 - 2*eta_q*eta_t*y2S*z2S*zpq - 2*eta_q*eta_t*yqd*z2S*zpq - 2*eta_q*y2S*z1D*zpq - eta_q*y2S*zpa*zpq - eta_q*y2S*zpb*zpq - 2*eta_q*yqd*z1D*zpq - eta_q*yqd*zpa*zpq - eta_q*yqd*zpb*zpq + eta_t2*y2S*z2S2 + eta_t2*yqd*z2S2 + 2*eta_t*y2S*z1D*z2S + eta_t*y2S*z2S*zpa + eta_t*y2S*z2S*zpb + 2*eta_t*yqd*z1D*z2S + eta_t*yqd*z2S*zpa + eta_t*yqd*z2S*zpb + y2S*z1D2 + y2S*z1D*zpa + y2S*z1D*zpb + y2S*zpa*zpb + yqd*z1D2 + yqd*z1D*zpa + yqd*z1D*zpb + yqd*zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3313) ! | pz  pz  px  pz   ( 248) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*zpq3 - 2*eta_p*eta_q*eta_t*z2S*zpq2 - 2*eta_p*eta_q*z1D*zpq2 - eta_p*eta_q*zpa*zpq2 - eta_p*eta_q*zpb*zpq2 + eta_p*eta_t2*z2S2*zpq + 2*eta_p*eta_t*z1D*z2S*zpq + eta_p*eta_t*z2S*zpa*zpq + eta_p*eta_t*z2S*zpb*zpq + eta_p*z1D2*zpq + eta_p*z1D*zpa*zpq + eta_p*z1D*zpb*zpq + eta_p*zpa*zpb*zpq + eta_q2*z2S*zpq2 + eta_q2*zpq2*zqd - 2*eta_q*eta_t*z2S2*zpq - 2*eta_q*eta_t*z2S*zpq*zqd - 2*eta_q*z1D*z2S*zpq - 2*eta_q*z1D*zpq*zqd - eta_q*z2S*zpa*zpq - eta_q*z2S*zpb*zpq - eta_q*zpa*zpq*zqd - eta_q*zpb*zpq*zqd + eta_t2*z2S3 + eta_t2*z2S2*zqd + 2*eta_t*z1D*z2S2 + 2*eta_t*z1D*z2S*zqd + eta_t*z2S2*zpa + eta_t*z2S2*zpb + eta_t*z2S*zpa*zqd + eta_t*z2S*zpb*zqd + z1D2*z2S + z1D2*zqd + z1D*z2S*zpa + z1D*z2S*zpb + z1D*zpa*zqd + z1D*zpb*zqd + z2S*zpa*zpb + zpa*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqc * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqc * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3320) ! | pz  pz  py  s    ( 249) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*ypq*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S + eta_p*eta_t*ypq*z2S*zpa + eta_p*eta_t*ypq*z2S*zpb + eta_p*ypq*z1D2 + eta_p*ypq*z1D*zpa + eta_p*ypq*z1D*zpb + eta_p*ypq*zpa*zpb + eta_q2*y2S*zpq2 + eta_q2*yqc*zpq2 - 2*eta_q*eta_t*y2S*z2S*zpq - 2*eta_q*eta_t*yqc*z2S*zpq - 2*eta_q*y2S*z1D*zpq - eta_q*y2S*zpa*zpq - eta_q*y2S*zpb*zpq - 2*eta_q*yqc*z1D*zpq - eta_q*yqc*zpa*zpq - eta_q*yqc*zpb*zpq + eta_t2*y2S*z2S2 + eta_t2*yqc*z2S2 + 2*eta_t*y2S*z1D*z2S + eta_t*y2S*z2S*zpa + eta_t*y2S*z2S*zpb + 2*eta_t*yqc*z1D*z2S + eta_t*yqc*z2S*zpa + eta_t*yqc*z2S*zpb + y2S*z1D2 + y2S*z1D*zpa + y2S*z1D*zpb + y2S*zpa*zpb + yqc*z1D2 + yqc*z1D*zpa + yqc*z1D*zpb + yqc*zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3321) ! | pz  pz  py  px   ( 250) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*ypq*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq - eta_p*eta_q*ypq*zpa*zpq - eta_p*eta_q*ypq*zpb*zpq + eta_p*eta_t2*ypq*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S + eta_p*eta_t*ypq*z2S*zpa + eta_p*eta_t*ypq*z2S*zpb + eta_p*ypq*z1D2 + eta_p*ypq*z1D*zpa + eta_p*ypq*z1D*zpb + eta_p*ypq*zpa*zpb + eta_q2*y2S*zpq2 + eta_q2*yqc*zpq2 - 2*eta_q*eta_t*y2S*z2S*zpq - 2*eta_q*eta_t*yqc*z2S*zpq - 2*eta_q*y2S*z1D*zpq - eta_q*y2S*zpa*zpq - eta_q*y2S*zpb*zpq - 2*eta_q*yqc*z1D*zpq - eta_q*yqc*zpa*zpq - eta_q*yqc*zpb*zpq + eta_t2*y2S*z2S2 + eta_t2*yqc*z2S2 + 2*eta_t*y2S*z1D*z2S + eta_t*y2S*z2S*zpa + eta_t*y2S*z2S*zpb + 2*eta_t*yqc*z1D*z2S + eta_t*yqc*z2S*zpa + eta_t*yqc*z2S*zpb + y2S*z1D2 + y2S*z1D*zpa + y2S*z1D*zpb + y2S*zpa*zpb + yqc*z1D2 + yqc*z1D*zpa + yqc*z1D*zpb + yqc*zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3322) ! | pz  pz  py  py   ( 251) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq2*zpq2 - 2*eta_p2*eta_q*eta_t*ypq2*z2S*zpq - 2*eta_p2*eta_q*ypq2*z1D*zpq - eta_p2*eta_q*ypq2*zpa*zpq - eta_p2*eta_q*ypq2*zpb*zpq + eta_p2*eta_t2*ypq2*z2S2 + 2*eta_p2*eta_t*ypq2*z1D*z2S + eta_p2*eta_t*ypq2*z2S*zpa + eta_p2*eta_t*ypq2*z2S*zpb + eta_p2*ypq2*z1D2 + eta_p2*ypq2*z1D*zpa + eta_p2*ypq2*z1D*zpb + eta_p2*ypq2*zpa*zpb + 2*eta_p*eta_q2*y2S*ypq*zpq2 + eta_p*eta_q2*ypq*yqc*zpq2 + eta_p*eta_q2*ypq*yqd*zpq2 - 4*eta_p*eta_q*eta_t*y2S*ypq*z2S*zpq - 2*eta_p*eta_q*eta_t*ypq*yqc*z2S*zpq - 2*eta_p*eta_q*eta_t*ypq*yqd*z2S*zpq - 4*eta_p*eta_q*y2S*ypq*z1D*zpq - 2*eta_p*eta_q*y2S*ypq*zpa*zpq - 2*eta_p*eta_q*y2S*ypq*zpb*zpq - 2*eta_p*eta_q*ypq*yqc*z1D*zpq - eta_p*eta_q*ypq*yqc*zpa*zpq - eta_p*eta_q*ypq*yqc*zpb*zpq - 2*eta_p*eta_q*ypq*yqd*z1D*zpq - eta_p*eta_q*ypq*yqd*zpa*zpq - eta_p*eta_q*ypq*yqd*zpb*zpq + 2*eta_p*eta_t2*y2S*ypq*z2S2 + eta_p*eta_t2*ypq*yqc*z2S2 + eta_p*eta_t2*ypq*yqd*z2S2 + 4*eta_p*eta_t*y2S*ypq*z1D*z2S + 2*eta_p*eta_t*y2S*ypq*z2S*zpa + 2*eta_p*eta_t*y2S*ypq*z2S*zpb + 2*eta_p*eta_t*ypq*yqc*z1D*z2S + eta_p*eta_t*ypq*yqc*z2S*zpa + eta_p*eta_t*ypq*yqc*z2S*zpb + 2*eta_p*eta_t*ypq*yqd*z1D*z2S + eta_p*eta_t*ypq*yqd*z2S*zpa + eta_p*eta_t*ypq*yqd*z2S*zpb + 2*eta_p*y2S*ypq*z1D2 + 2*eta_p*y2S*ypq*z1D*zpa + 2*eta_p*y2S*ypq*z1D*zpb + 2*eta_p*y2S*ypq*zpa*zpb + eta_p*ypq*yqc*z1D2 + eta_p*ypq*yqc*z1D*zpa + eta_p*ypq*yqc*z1D*zpb + eta_p*ypq*yqc*zpa*zpb + eta_p*ypq*yqd*z1D2 + eta_p*ypq*yqd*z1D*zpa + eta_p*ypq*yqd*z1D*zpb + eta_p*ypq*yqd*zpa*zpb + eta_q2*y2S2*zpq2 + eta_q2*y2S*yqc*zpq2 + eta_q2*y2S*yqd*zpq2 + eta_q2*yqc*yqd*zpq2 - 2*eta_q*eta_t*y2S2*z2S*zpq - 2*eta_q*eta_t*y2S*yqc*z2S*zpq - 2*eta_q*eta_t*y2S*yqd*z2S*zpq - 2*eta_q*eta_t*yqc*yqd*z2S*zpq - 2*eta_q*y2S2*z1D*zpq - eta_q*y2S2*zpa*zpq - eta_q*y2S2*zpb*zpq - 2*eta_q*y2S*yqc*z1D*zpq - eta_q*y2S*yqc*zpa*zpq - eta_q*y2S*yqc*zpb*zpq - 2*eta_q*y2S*yqd*z1D*zpq - eta_q*y2S*yqd*zpa*zpq - eta_q*y2S*yqd*zpb*zpq - 2*eta_q*yqc*yqd*z1D*zpq - eta_q*yqc*yqd*zpa*zpq - eta_q*yqc*yqd*zpb*zpq + eta_t2*y2S2*z2S2 + eta_t2*y2S*yqc*z2S2 + eta_t2*y2S*yqd*z2S2 + eta_t2*yqc*yqd*z2S2 + 2*eta_t*y2S2*z1D*z2S + eta_t*y2S2*z2S*zpa + eta_t*y2S2*z2S*zpb + 2*eta_t*y2S*yqc*z1D*z2S + eta_t*y2S*yqc*z2S*zpa + eta_t*y2S*yqc*z2S*zpb + 2*eta_t*y2S*yqd*z1D*z2S + eta_t*y2S*yqd*z2S*zpa + eta_t*y2S*yqd*z2S*zpb + 2*eta_t*yqc*yqd*z1D*z2S + eta_t*yqc*yqd*z2S*zpa + eta_t*yqc*yqd*z2S*zpb + y2S2*z1D2 + y2S2*z1D*zpa + y2S2*z1D*zpb + y2S2*zpa*zpb + y2S*yqc*z1D2 + y2S*yqc*z1D*zpa + y2S*yqc*z1D*zpb + y2S*yqc*zpa*zpb + y2S*yqd*z1D2 + y2S*yqd*z1D*zpa + y2S*yqd*z1D*zpb + y2S*yqd*zpa*zpb + yqc*yqd*z1D2 + yqc*yqd*z1D*zpa + yqc*yqd*z1D*zpb + yqc*yqd*zpa*zpb) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3323) ! | pz  pz  py  pz   ( 252) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq*zpq3 - 2*eta_p2*eta_q*eta_t*ypq*z2S*zpq2 - 2*eta_p2*eta_q*ypq*z1D*zpq2 - eta_p2*eta_q*ypq*zpa*zpq2 - eta_p2*eta_q*ypq*zpb*zpq2 + eta_p2*eta_t2*ypq*z2S2*zpq + 2*eta_p2*eta_t*ypq*z1D*z2S*zpq + eta_p2*eta_t*ypq*z2S*zpa*zpq + eta_p2*eta_t*ypq*z2S*zpb*zpq + eta_p2*ypq*z1D2*zpq + eta_p2*ypq*z1D*zpa*zpq + eta_p2*ypq*z1D*zpb*zpq + eta_p2*ypq*zpa*zpb*zpq + eta_p*eta_q2*y2S*zpq3 + eta_p*eta_q2*ypq*z2S*zpq2 + eta_p*eta_q2*ypq*zpq2*zqd + eta_p*eta_q2*yqc*zpq3 - 2*eta_p*eta_q*eta_t*y2S*z2S*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S2*zpq - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq*zqd - 2*eta_p*eta_q*eta_t*yqc*z2S*zpq2 - 2*eta_p*eta_q*y2S*z1D*zpq2 - eta_p*eta_q*y2S*zpa*zpq2 - eta_p*eta_q*y2S*zpb*zpq2 - 2*eta_p*eta_q*ypq*z1D*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq*zqd - eta_p*eta_q*ypq*z2S*zpa*zpq - eta_p*eta_q*ypq*z2S*zpb*zpq - eta_p*eta_q*ypq*zpa*zpq*zqd - eta_p*eta_q*ypq*zpb*zpq*zqd - 2*eta_p*eta_q*yqc*z1D*zpq2 - eta_p*eta_q*yqc*zpa*zpq2 - eta_p*eta_q*yqc*zpb*zpq2 + eta_p*eta_t2*y2S*z2S2*zpq + eta_p*eta_t2*ypq*z2S3 + eta_p*eta_t2*ypq*z2S2*zqd + eta_p*eta_t2*yqc*z2S2*zpq + 2*eta_p*eta_t*y2S*z1D*z2S*zpq + eta_p*eta_t*y2S*z2S*zpa*zpq + eta_p*eta_t*y2S*z2S*zpb*zpq + 2*eta_p*eta_t*ypq*z1D*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S*zqd + eta_p*eta_t*ypq*z2S2*zpa + eta_p*eta_t*ypq*z2S2*zpb + eta_p*eta_t*ypq*z2S*zpa*zqd + eta_p*eta_t*ypq*z2S*zpb*zqd + 2*eta_p*eta_t*yqc*z1D*z2S*zpq + eta_p*eta_t*yqc*z2S*zpa*zpq + eta_p*eta_t*yqc*z2S*zpb*zpq + eta_p*y2S*z1D2*zpq + eta_p*y2S*z1D*zpa*zpq + eta_p*y2S*z1D*zpb*zpq + eta_p*y2S*zpa*zpb*zpq + eta_p*ypq*z1D2*z2S + eta_p*ypq*z1D2*zqd + eta_p*ypq*z1D*z2S*zpa + eta_p*ypq*z1D*z2S*zpb + eta_p*ypq*z1D*zpa*zqd + eta_p*ypq*z1D*zpb*zqd + eta_p*ypq*z2S*zpa*zpb + eta_p*ypq*zpa*zpb*zqd + eta_p*yqc*z1D2*zpq + eta_p*yqc*z1D*zpa*zpq + eta_p*yqc*z1D*zpb*zpq + eta_p*yqc*zpa*zpb*zpq + eta_q2*y2S*z2S*zpq2 + eta_q2*y2S*zpq2*zqd + eta_q2*yqc*z2S*zpq2 + eta_q2*yqc*zpq2*zqd - 2*eta_q*eta_t*y2S*z2S2*zpq - 2*eta_q*eta_t*y2S*z2S*zpq*zqd - 2*eta_q*eta_t*yqc*z2S2*zpq - 2*eta_q*eta_t*yqc*z2S*zpq*zqd - 2*eta_q*y2S*z1D*z2S*zpq - 2*eta_q*y2S*z1D*zpq*zqd - eta_q*y2S*z2S*zpa*zpq - eta_q*y2S*z2S*zpb*zpq - eta_q*y2S*zpa*zpq*zqd - eta_q*y2S*zpb*zpq*zqd - 2*eta_q*yqc*z1D*z2S*zpq - 2*eta_q*yqc*z1D*zpq*zqd - eta_q*yqc*z2S*zpa*zpq - eta_q*yqc*z2S*zpb*zpq - eta_q*yqc*zpa*zpq*zqd - eta_q*yqc*zpb*zpq*zqd + eta_t2*y2S*z2S3 + eta_t2*y2S*z2S2*zqd + eta_t2*yqc*z2S3 + eta_t2*yqc*z2S2*zqd + 2*eta_t*y2S*z1D*z2S2 + 2*eta_t*y2S*z1D*z2S*zqd + eta_t*y2S*z2S2*zpa + eta_t*y2S*z2S2*zpb + eta_t*y2S*z2S*zpa*zqd + eta_t*y2S*z2S*zpb*zqd + 2*eta_t*yqc*z1D*z2S2 + 2*eta_t*yqc*z1D*z2S*zqd + eta_t*yqc*z2S2*zpa + eta_t*yqc*z2S2*zpb + eta_t*yqc*z2S*zpa*zqd + eta_t*yqc*z2S*zpb*zqd + y2S*z1D2*z2S + y2S*z1D2*zqd + y2S*z1D*z2S*zpa + y2S*z1D*z2S*zpb + y2S*z1D*zpa*zqd + y2S*z1D*zpb*zqd + y2S*z2S*zpa*zpb + y2S*zpa*zpb*zqd + yqc*z1D2*z2S + yqc*z1D2*zqd + yqc*z1D*z2S*zpa + yqc*z1D*z2S*zpb + yqc*z1D*zpa*zqd + yqc*z1D*zpb*zqd + yqc*z2S*zpa*zpb + yqc*zpa*zpb*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3330) ! | pz  pz  pz  s    ( 253) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*zpq3 - 2*eta_p*eta_q*eta_t*z2S*zpq2 - 2*eta_p*eta_q*z1D*zpq2 - eta_p*eta_q*zpa*zpq2 - eta_p*eta_q*zpb*zpq2 + eta_p*eta_t2*z2S2*zpq + 2*eta_p*eta_t*z1D*z2S*zpq + eta_p*eta_t*z2S*zpa*zpq + eta_p*eta_t*z2S*zpb*zpq + eta_p*z1D2*zpq + eta_p*z1D*zpa*zpq + eta_p*z1D*zpb*zpq + eta_p*zpa*zpb*zpq + eta_q2*z2S*zpq2 + eta_q2*zpq2*zqc - 2*eta_q*eta_t*z2S2*zpq - 2*eta_q*eta_t*z2S*zpq*zqc - 2*eta_q*z1D*z2S*zpq - 2*eta_q*z1D*zpq*zqc - eta_q*z2S*zpa*zpq - eta_q*z2S*zpb*zpq - eta_q*zpa*zpq*zqc - eta_q*zpb*zpq*zqc + eta_t2*z2S3 + eta_t2*z2S2*zqc + 2*eta_t*z1D*z2S2 + 2*eta_t*z1D*z2S*zqc + eta_t*z2S2*zpa + eta_t*z2S2*zpb + eta_t*z2S*zpa*zqc + eta_t*z2S*zpb*zqc + z1D2*z2S + z1D2*zqc + z1D*z2S*zpa + z1D*z2S*zpb + z1D*zpa*zqc + z1D*zpb*zqc + z2S*zpa*zpb + zpa*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3331) ! | pz  pz  pz  px   ( 254) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p*eta_q2*zpq3 - 2*eta_p*eta_q*eta_t*z2S*zpq2 - 2*eta_p*eta_q*z1D*zpq2 - eta_p*eta_q*zpa*zpq2 - eta_p*eta_q*zpb*zpq2 + eta_p*eta_t2*z2S2*zpq + 2*eta_p*eta_t*z1D*z2S*zpq + eta_p*eta_t*z2S*zpa*zpq + eta_p*eta_t*z2S*zpb*zpq + eta_p*z1D2*zpq + eta_p*z1D*zpa*zpq + eta_p*z1D*zpb*zpq + eta_p*zpa*zpb*zpq + eta_q2*z2S*zpq2 + eta_q2*zpq2*zqc - 2*eta_q*eta_t*z2S2*zpq - 2*eta_q*eta_t*z2S*zpq*zqc - 2*eta_q*z1D*z2S*zpq - 2*eta_q*z1D*zpq*zqc - eta_q*z2S*zpa*zpq - eta_q*z2S*zpb*zpq - eta_q*zpa*zpq*zqc - eta_q*zpb*zpq*zqc + eta_t2*z2S3 + eta_t2*z2S2*zqc + 2*eta_t*z1D*z2S2 + 2*eta_t*z1D*z2S*zqc + eta_t*z2S2*zpa + eta_t*z2S2*zpb + eta_t*z2S*zpa*zqc + eta_t*z2S*zpb*zqc + z1D2*z2S + z1D2*zqc + z1D*z2S*zpa + z1D*z2S*zpb + z1D*zpa*zqc + z1D*zpb*zqc + z2S*zpa*zpb + zpa*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * inv_ax * 0.5d0 * ( sxqd * ( I_B(abs(0-1)) + I_B(0+1) ) ) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * ( I_B(abs(n-1)) - I_B(n+1) ) + sxqd * ( I_B(abs(n-1)) + I_B(n+1) ) )
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3332) ! | pz  pz  pz  py   ( 255) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*ypq*zpq3 - 2*eta_p2*eta_q*eta_t*ypq*z2S*zpq2 - 2*eta_p2*eta_q*ypq*z1D*zpq2 - eta_p2*eta_q*ypq*zpa*zpq2 - eta_p2*eta_q*ypq*zpb*zpq2 + eta_p2*eta_t2*ypq*z2S2*zpq + 2*eta_p2*eta_t*ypq*z1D*z2S*zpq + eta_p2*eta_t*ypq*z2S*zpa*zpq + eta_p2*eta_t*ypq*z2S*zpb*zpq + eta_p2*ypq*z1D2*zpq + eta_p2*ypq*z1D*zpa*zpq + eta_p2*ypq*z1D*zpb*zpq + eta_p2*ypq*zpa*zpb*zpq + eta_p*eta_q2*y2S*zpq3 + eta_p*eta_q2*ypq*z2S*zpq2 + eta_p*eta_q2*ypq*zpq2*zqc + eta_p*eta_q2*yqd*zpq3 - 2*eta_p*eta_q*eta_t*y2S*z2S*zpq2 - 2*eta_p*eta_q*eta_t*ypq*z2S2*zpq - 2*eta_p*eta_q*eta_t*ypq*z2S*zpq*zqc - 2*eta_p*eta_q*eta_t*yqd*z2S*zpq2 - 2*eta_p*eta_q*y2S*z1D*zpq2 - eta_p*eta_q*y2S*zpa*zpq2 - eta_p*eta_q*y2S*zpb*zpq2 - 2*eta_p*eta_q*ypq*z1D*z2S*zpq - 2*eta_p*eta_q*ypq*z1D*zpq*zqc - eta_p*eta_q*ypq*z2S*zpa*zpq - eta_p*eta_q*ypq*z2S*zpb*zpq - eta_p*eta_q*ypq*zpa*zpq*zqc - eta_p*eta_q*ypq*zpb*zpq*zqc - 2*eta_p*eta_q*yqd*z1D*zpq2 - eta_p*eta_q*yqd*zpa*zpq2 - eta_p*eta_q*yqd*zpb*zpq2 + eta_p*eta_t2*y2S*z2S2*zpq + eta_p*eta_t2*ypq*z2S3 + eta_p*eta_t2*ypq*z2S2*zqc + eta_p*eta_t2*yqd*z2S2*zpq + 2*eta_p*eta_t*y2S*z1D*z2S*zpq + eta_p*eta_t*y2S*z2S*zpa*zpq + eta_p*eta_t*y2S*z2S*zpb*zpq + 2*eta_p*eta_t*ypq*z1D*z2S2 + 2*eta_p*eta_t*ypq*z1D*z2S*zqc + eta_p*eta_t*ypq*z2S2*zpa + eta_p*eta_t*ypq*z2S2*zpb + eta_p*eta_t*ypq*z2S*zpa*zqc + eta_p*eta_t*ypq*z2S*zpb*zqc + 2*eta_p*eta_t*yqd*z1D*z2S*zpq + eta_p*eta_t*yqd*z2S*zpa*zpq + eta_p*eta_t*yqd*z2S*zpb*zpq + eta_p*y2S*z1D2*zpq + eta_p*y2S*z1D*zpa*zpq + eta_p*y2S*z1D*zpb*zpq + eta_p*y2S*zpa*zpb*zpq + eta_p*ypq*z1D2*z2S + eta_p*ypq*z1D2*zqc + eta_p*ypq*z1D*z2S*zpa + eta_p*ypq*z1D*z2S*zpb + eta_p*ypq*z1D*zpa*zqc + eta_p*ypq*z1D*zpb*zqc + eta_p*ypq*z2S*zpa*zpb + eta_p*ypq*zpa*zpb*zqc + eta_p*yqd*z1D2*zpq + eta_p*yqd*z1D*zpa*zpq + eta_p*yqd*z1D*zpb*zpq + eta_p*yqd*zpa*zpb*zpq + eta_q2*y2S*z2S*zpq2 + eta_q2*y2S*zpq2*zqc + eta_q2*yqd*z2S*zpq2 + eta_q2*yqd*zpq2*zqc - 2*eta_q*eta_t*y2S*z2S2*zpq - 2*eta_q*eta_t*y2S*z2S*zpq*zqc - 2*eta_q*eta_t*yqd*z2S2*zpq - 2*eta_q*eta_t*yqd*z2S*zpq*zqc - 2*eta_q*y2S*z1D*z2S*zpq - 2*eta_q*y2S*z1D*zpq*zqc - eta_q*y2S*z2S*zpa*zpq - eta_q*y2S*z2S*zpb*zpq - eta_q*y2S*zpa*zpq*zqc - eta_q*y2S*zpb*zpq*zqc - 2*eta_q*yqd*z1D*z2S*zpq - 2*eta_q*yqd*z1D*zpq*zqc - eta_q*yqd*z2S*zpa*zpq - eta_q*yqd*z2S*zpb*zpq - eta_q*yqd*zpa*zpq*zqc - eta_q*yqd*zpb*zpq*zqc + eta_t2*y2S*z2S3 + eta_t2*y2S*z2S2*zqc + eta_t2*yqd*z2S3 + eta_t2*yqd*z2S2*zqc + 2*eta_t*y2S*z1D*z2S2 + 2*eta_t*y2S*z1D*z2S*zqc + eta_t*y2S*z2S2*zpa + eta_t*y2S*z2S2*zpb + eta_t*y2S*z2S*zpa*zqc + eta_t*y2S*z2S*zpb*zqc + 2*eta_t*yqd*z1D*z2S2 + 2*eta_t*yqd*z1D*z2S*zqc + eta_t*yqd*z2S2*zpa + eta_t*yqd*z2S2*zpb + eta_t*yqd*z2S*zpa*zqc + eta_t*yqd*z2S*zpb*zqc + y2S*z1D2*z2S + y2S*z1D2*zqc + y2S*z1D*z2S*zpa + y2S*z1D*z2S*zpb + y2S*z1D*zpa*zqc + y2S*z1D*zpb*zqc + y2S*z2S*zpa*zpb + y2S*zpa*zpb*zqc + yqd*z1D2*z2S + yqd*z1D2*zqc + yqd*z1D*z2S*zpa + yqd*z1D*z2S*zpb + yqd*z1D*zpa*zqc + yqd*z1D*zpb*zqc + yqd*z2S*zpa*zpb + yqd*zpa*zpb*zqc) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case (3333) ! | pz  pz  pz  pz   ( 256) 
      
      n         = 0
      const     = pi2 * D2 * (eta_p2*eta_q2*zpq4 - 2*eta_p2*eta_q*eta_t*z2S*zpq3 - 2*eta_p2*eta_q*z1D*zpq3 - eta_p2*eta_q*zpa*zpq3 - eta_p2*eta_q*zpb*zpq3 + eta_p2*eta_t2*z2S2*zpq2 + 2*eta_p2*eta_t*z1D*z2S*zpq2 + eta_p2*eta_t*z2S*zpa*zpq2 + eta_p2*eta_t*z2S*zpb*zpq2 + eta_p2*z1D2*zpq2 + eta_p2*z1D*zpa*zpq2 + eta_p2*z1D*zpb*zpq2 + eta_p2*zpa*zpb*zpq2 + 2*eta_p*eta_q2*z2S*zpq3 + eta_p*eta_q2*zpq3*zqc + eta_p*eta_q2*zpq3*zqd - 4*eta_p*eta_q*eta_t*z2S2*zpq2 - 2*eta_p*eta_q*eta_t*z2S*zpq2*zqc - 2*eta_p*eta_q*eta_t*z2S*zpq2*zqd - 4*eta_p*eta_q*z1D*z2S*zpq2 - 2*eta_p*eta_q*z1D*zpq2*zqc - 2*eta_p*eta_q*z1D*zpq2*zqd - 2*eta_p*eta_q*z2S*zpa*zpq2 - 2*eta_p*eta_q*z2S*zpb*zpq2 - eta_p*eta_q*zpa*zpq2*zqc - eta_p*eta_q*zpa*zpq2*zqd - eta_p*eta_q*zpb*zpq2*zqc - eta_p*eta_q*zpb*zpq2*zqd + 2*eta_p*eta_t2*z2S3*zpq + eta_p*eta_t2*z2S2*zpq*zqc + eta_p*eta_t2*z2S2*zpq*zqd + 4*eta_p*eta_t*z1D*z2S2*zpq + 2*eta_p*eta_t*z1D*z2S*zpq*zqc + 2*eta_p*eta_t*z1D*z2S*zpq*zqd + 2*eta_p*eta_t*z2S2*zpa*zpq + 2*eta_p*eta_t*z2S2*zpb*zpq + eta_p*eta_t*z2S*zpa*zpq*zqc + eta_p*eta_t*z2S*zpa*zpq*zqd + eta_p*eta_t*z2S*zpb*zpq*zqc + eta_p*eta_t*z2S*zpb*zpq*zqd + 2*eta_p*z1D2*z2S*zpq + eta_p*z1D2*zpq*zqc + eta_p*z1D2*zpq*zqd + 2*eta_p*z1D*z2S*zpa*zpq + 2*eta_p*z1D*z2S*zpb*zpq + eta_p*z1D*zpa*zpq*zqc + eta_p*z1D*zpa*zpq*zqd + eta_p*z1D*zpb*zpq*zqc + eta_p*z1D*zpb*zpq*zqd + 2*eta_p*z2S*zpa*zpb*zpq + eta_p*zpa*zpb*zpq*zqc + eta_p*zpa*zpb*zpq*zqd + eta_q2*z2S2*zpq2 + eta_q2*z2S*zpq2*zqc + eta_q2*z2S*zpq2*zqd + eta_q2*zpq2*zqc*zqd - 2*eta_q*eta_t*z2S3*zpq - 2*eta_q*eta_t*z2S2*zpq*zqc - 2*eta_q*eta_t*z2S2*zpq*zqd - 2*eta_q*eta_t*z2S*zpq*zqc*zqd - 2*eta_q*z1D*z2S2*zpq - 2*eta_q*z1D*z2S*zpq*zqc - 2*eta_q*z1D*z2S*zpq*zqd - 2*eta_q*z1D*zpq*zqc*zqd - eta_q*z2S2*zpa*zpq - eta_q*z2S2*zpb*zpq - eta_q*z2S*zpa*zpq*zqc - eta_q*z2S*zpa*zpq*zqd - eta_q*z2S*zpb*zpq*zqc - eta_q*z2S*zpb*zpq*zqd - eta_q*zpa*zpq*zqc*zqd - eta_q*zpb*zpq*zqc*zqd + eta_t2*z2S4 + eta_t2*z2S3*zqc + eta_t2*z2S3*zqd + eta_t2*z2S2*zqc*zqd + 2*eta_t*z1D*z2S3 + 2*eta_t*z1D*z2S2*zqc + 2*eta_t*z1D*z2S2*zqd + 2*eta_t*z1D*z2S*zqc*zqd + eta_t*z2S3*zpa + eta_t*z2S3*zpb + eta_t*z2S2*zpa*zqc + eta_t*z2S2*zpa*zqd + eta_t*z2S2*zpb*zqc + eta_t*z2S2*zpb*zqd + eta_t*z2S*zpa*zqc*zqd + eta_t*z2S*zpb*zqc*zqd + z1D2*z2S2 + z1D2*z2S*zqc + z1D2*z2S*zqd + z1D2*zqc*zqd + z1D*z2S2*zpa + z1D*z2S2*zpb + z1D*z2S*zpa*zqc + z1D*z2S*zpa*zqd + z1D*z2S*zpb*zqc + z1D*z2S*zpb*zqd + z1D*zpa*zqc*zqd + z1D*zpb*zqc*zqd + z2S2*zpa*zpb + z2S*zpa*zpb*zqc + z2S*zpa*zpb*zqd + zpa*zpb*zqc*zqd) * dexp(-p * q * t2 * D2 * ( ypq*ypq + zpq*zpq ) ) 
      sum       = I_A(0) * I_B(0) * I_C(0) * const
      do n      = 1 , Nmax
        termAn  = I_A(n)
        termBn  = I_B(n)
        termc   = I_C(n) 
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
          sum = sum + 2.d0 * real(term) * const
          current_term = current_term * expo_term
      end do

      case default
        sum = 0.0d0
      end select

end function S


end subroutine integrate_ERI_sum



subroutine bessel_I_scaled_backward(Nmax, x, I_scaled)

      use bessel_functions

      implicit none
      integer, intent(in)           :: Nmax
      double precision, intent(in)  :: x
      double precision, intent(out) :: I_scaled(0:Nmax)
      integer                       :: n
      double precision              :: scale_factor , I_0_normalized

      ! ============================================================
      ! SPECIAL CASE: x = 0
      ! I_0(0) = 1, I_n(0) = 0 for n > 0
      ! Scaled versions: exp(-0) * I_n(0) = I_n(0)
      ! ============================================================
      if (x < 1.d-10) then
        I_scaled(0) = 1.d0
        do n = 1, Nmax
          I_scaled(n) = 0.d0
        end do
        return
      end if      
      ! ============================================================
      ! NORMAL CASE: x > 0
      ! ============================================================
      
      ! Initialize
      I_scaled(Nmax)    = iv_scaled(Nmax, x)
      I_scaled(Nmax-1)  = iv_scaled(Nmax-1, x)
      I_0_normalized    = iv_scaled(0, x)
      
      ! Backward recursion
      do n = Nmax-1, 1, -1
        I_scaled(n-1) = I_scaled(n+1) + (2.d0 * dble(n) / x) * I_scaled(n)
      end do
      
      ! Normalize
      scale_factor = I_0_normalized / I_scaled(0)
      
      ! Apply normalization
      I_scaled(0:Nmax) = I_scaled(0:Nmax) * scale_factor  ! Array notation

      
end subroutine bessel_I_scaled_backward