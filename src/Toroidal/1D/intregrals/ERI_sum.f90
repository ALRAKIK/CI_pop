subroutine integrate_ERI_sum(pattern_id,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,result)
      
      use quadpack , only : dqagi , dqags , dqag
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod
      use bessel_functions
      use files

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

      double precision,parameter         :: epsabs = 1.0d-10 , epsrel = 1.0d-8
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      
      integer                            :: Nmax
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi

      double precision                   :: inv_ax 
      double precision                   :: inv_ax2

      double precision                   :: cxpa   , sxpa
      double precision                   :: cxpb   , sxpb
      double precision                   :: cxqc   , sxqc
      double precision                   :: cxqd   , sxqd
      double precision                   :: c2xpab , s2xpab
      double precision                   :: c2xqcd , s2xqcd

      double precision                   :: A , B , C
      integer                            :: Peak

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


      ! --------------------------------------------------------------- !
      
      ! --------------------------------------------------------------- !

      t2 = t  * t 
      t3 = t2 * t 
      t4 = t3 * t 

      C   = 2.d0  * t * t * inv_ax2

      D   = 1.d0/(dsqrt(p*q+(p+q)*t2))

      D2  = D * D

      termAn = 0.d0 
      termBn = 0.d0 
      termC  = 0.d0

      !               the peak of the sum              !
      Peak         = ceiling(min(A,B,C))
      !         the maximum of terms in the sum        !
      Nmax         = Peak+10
      !                    Phase term                  !
      expo_term    = exp(I_dp*phi)
      !                   Initial term                 !
      current_term = expo_term
      
      select case(pattern_id)
      
      ! /////////////////////////////////////////////////////////////// !


      case (0000) ! | s   s   s   s    ( 1 ) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C)
        term    = current_term * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0001) ! | s   s   s   px   ( 2 ) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C)
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0010) ! | s   s   px  s    ( 5 ) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0011) ! | s   s   px  px   ( 6 ) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if 
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0022) ! | s   s   py  py   ( 11) 
      n           = 0
      const       = ( 0.5d0 * (p+t2) * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
      current_term = current_term * expo_term
      end do


      case (0033) ! | s   s   pz  pz   ( 16) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * (p+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0100) ! | s   px  s   s    ( 17) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0101) ! | s   px  s   px   ( 18) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (0110) ! | s   px  px  s    ( 21) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0111) ! | s   px  px  px   ( 22) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if 
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0122) ! | s   px  py  py   ( 27) 
      n           = 0
      const       = ( 0.5d0 * (p+t2) * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
      current_term = current_term * expo_term
      end do

      case (0133) ! | s   px  pz  pz   ( 32) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * (p+t2) * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0202) ! | s   py  s   py   ( 35) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0212) ! | s   py  px  py   ( 39) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (0220) ! | s   py  py  s    ( 41) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0221) ! | s   py  py  px   ( 42) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0303) ! | s   pz  s   pz   ( 52) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (0313) ! | s   pz  px  pz   ( 56) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0330) ! | s   pz  pz  s    ( 61) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (0331) ! | s   pz  pz  px   ( 62) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1000) ! | px  s   s   s    ( 65) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1001) ! | px  s   s   px   ( 66) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (1010) ! | px  s   px  s    ( 69) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1011) ! | px  s   px  px   ( 70) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if 
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1022) ! | px  s   py  py   ( 75) 
      n           = 0
      const       = ( 0.5d0 * (p+t2) * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1033) ! | px  s   pz  pz   ( 80) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * (p+t2) * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1100) ! | px  px  s   s    ( 81) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then 
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1101) ! | px  px  s   px   ( 82) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then 
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1110) ! | px  px  px  s    ( 85) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then 
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1111) ! | px  px  px  px   ( 86) 
      n           = 0
      const       =  (pi * D)  *  (pi * D)
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then  
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc * cxqd * iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc * cxqd * iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1122) ! | px  px  py  py   ( 91) 
      n           = 0
      const       = ( 0.5d0 * (p+t2) * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then 
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1133) ! | px  px  pz  pz   ( 96) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * (p+t2) * D2 * pi * D )
      sum         = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        !if (dabs(p_x) < 1.d-15) then 
        !termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A))))
        !else 
        termAn  = inv_ax2 * (cxpa * cxpb * iv_scaled(n,A)-c2xpab*(0.25d0*(iv_scaled(n-2,A)+2.d0*iv_scaled(n,A)+iv_scaled(n+2,A)))+I_dp/A*n * s2xpab * (0.5d0*(iv_scaled(n-1,A)+iv_scaled(n+1,A))))
        !end if 
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1202) ! | px  py  s   py   ( 99) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (1212) ! | px  py  px  py   ( 103) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (1220) ! | px  py  py  s    ( 105) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1221) ! | px  py  py  px   ( 106) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1303) ! | px  pz  s   pz   ( 116) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1313) ! | px  pz  px  pz   ( 120) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1330) ! | px  pz  pz  s    ( 125) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (1331) ! | px  pz  pz  px   ( 126) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpa * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2002) ! | py  s   s   py   ( 131) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2012) ! | py  s   px  py   ( 135) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2020) ! | py  s   py  s    ( 137) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2021) ! | py  s   py  px   ( 138) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2102) ! | py  px  s   py   ( 147) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2112) ! | py  px  px  py   ( 151) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2120) ! | py  px  py  s    ( 153) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2121) ! | py  px  py  px   ( 154) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) *  (pi * D)
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (2200) ! | py  py  s   s    ( 161) 
      n           = 0
      const       = ( 0.5d0  * (q+t2) * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2201) ! | py  py  s   px   ( 162) 
      n           = 0
      const       = ( 0.5d0  * (q+t2) * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2210) ! | py  py  px  s    ( 165) 
      n           = 0
      const       = ( 0.5d0  * (q+t2) * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2211) ! | py  py  px  px   ( 166) 
      n           = 0
      const       = ( 0.5d0  * (q+t2) * D2 * pi * D ) *  (pi * D)
      sum         = iv_scaled(n, A) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if 
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (2222) ! | py  py  py  py   ( 171) 
      n           = 0
      const       =   (0.25d0 * ( 1.d0 + 3.d0 * t4 * D2 ) * D2 * pi * D )  *  (pi * D)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2233) ! | py  py  pz  pz   ( 176) 
      n           = 0
      const       = ( 0.5d0  * (q+t2) * D2 * pi * D ) * ( 0.5d0 * (p+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2323) ! | py  pz  py  pz   ( 188) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (2332) ! | py  pz  pz  py   ( 191)
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3003) ! | pz  s   s   pz   ( 196) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case (3013) ! | pz  s   px  pz   ( 200) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3030) ! | pz  s   pz  s    ( 205) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3031) ! | pz  s   pz  px   ( 206) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3103) ! | pz  px  s   pz   ( 212) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3113) ! | pz  px  px  pz   ( 216) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3130) ! | pz  px  pz  s    ( 221) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3131) ! | pz  px  pz  px   ( 222) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = inv_ax * 0.5d0 * (sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A))) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (iv_scaled(n-1,A) - iv_scaled(n+1,A)) + sxpb * (iv_scaled(n-1,A)+iv_scaled(n+1,A)))
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3223) ! | pz  py  py  pz   ( 236) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3232) ! | pz  py  pz  py   ( 239) 
      n           = 0
      const       = ( 0.5d0 * t2 * D2 * pi * D ) * ( 0.5d0 * t2 * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3300) ! | pz  pz  s   s    ( 241) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t2) * D2 * pi * D )  !* exp(A+B-2.d0*(p+q)*inv_ax2)
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3301) ! | pz  pz  s   px   ( 242) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqd * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3310) ! | pz  pz  px  s    ( 245) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax * 0.5d0 * (sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (iv_scaled(n-1,B)-iv_scaled(n+1,B)) + sxqc * (iv_scaled(n-1,B)+iv_scaled(n+1,B)))
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3311) ! | pz  pz  px  px   ( 246) 
      n           = 0
      const       =  (pi * D)  * ( 0.5d0  * (q+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        !if (dabs(q_x) < 1.d-15) then 
        !termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B))))
        !else 
        termBn  = inv_ax2 * (cxqc*cxqd*iv_scaled(n,B)-c2xqcd*(0.25d0*(iv_scaled(n-2,B)+2.d0*iv_scaled(n,B)+iv_scaled(n+2,B)))-I_dp/B*n * s2xqcd * (0.5d0*(iv_scaled(n-1,B)+iv_scaled(n+1,B))))
        !end if 
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3322) ! | pz  pz  py  py   ( 251) 
      n           = 0
      const       = ( 0.5d0 * (p+t2) * D2 * pi * D ) * ( 0.5d0  * (q+t2) * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do

      case (3333) ! | pz  pz  pz  pz   ( 256) 
      n           = 0
      const       =  (pi * D)  *   (0.25d0 * ( 1.d0 + 3.d0 * t4 * D2 ) * D2 * pi * D )
      sum         = iv_scaled(n, A) * iv_scaled(n, B) * iv_scaled(n, C) * const
      do n = 1 , Nmax
        termAn  = iv_scaled(n, A)
        termBn  = iv_scaled(n, B)
        termC   = iv_scaled(n, C) 
        term    = exp(I_dp*dble(n)*phi) * termC * termAn * termBn
        if (abs(term) < eps * dabs(sum) ) exit
        sum     = sum + 2.d0 * real(term) * const
        current_term = current_term * expo_term
      end do


      case default
        sum = 0.0d0
      end select

end function S


end subroutine integrate_ERI_sum