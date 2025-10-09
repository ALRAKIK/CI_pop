subroutine integrate_ERI_integral_mod(pattern_id,px_count,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,result)

      use quadpack , only : dqag , dqags, dqagp
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
      integer         , intent(in)       :: px_count 
      

      ! Output parameters

      double precision, intent(out)      :: result


      ! integration parameter ! 

      double precision,parameter         :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier,last, neval , iwork(limit)
      double precision                   :: abserr, work(lenw)
      integer, parameter                 :: key  = 6
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi

      ! ------------------------- local ------------------------------- !
      
      integer                            :: o1,o2,o3,o4 
      double precision                   :: A, A2, A3, A4, A5, A6, A7 
      double precision                   :: B, B2, B3, B4, B5, B6, B7 
      double precision                   :: z, z2, z3, z4, z5, z6, z7, z8
      double precision                   :: bz0, bz1, bz2, bz3, bz4 
      double precision                   :: ax2, ax3, ax4
      double precision                   :: term
      double precision                   :: const
      double precision                   :: erfcx
      double precision,parameter         :: eta = 1e-40

      double precision             :: spa   , cpa  , spb   , cpb
      double precision             :: sqc   , cqc  , sqd   , cqd  
      double precision             :: spab  , cpab , sqcd  , cqcd
      double precision             :: cab   , ccd
      double precision             :: psi 

      ! ------------------- integral variables ------------------------ !

      double precision             :: der_t , integral_t 
      double precision             :: int_xxxx , int_yyxx
      double precision             :: int_yxyx , int_xxyy
      double precision             :: int_yyyy , int_yyzz
      double precision             :: int_yzyz 
      double precision             :: beta, beta2, p2, q2, p3, q3 
      double precision             :: term1 , term2 ,term3
      double precision             :: inv   , Iinv 
      double precision             :: ppq , ppq2, ppq3 , ppq4 

      ! ------------------- derivative variables ---------------------- !

      double precision             :: d1_T    , d2_TT
      double precision             :: d1_A    , d2_AA 
      double precision             :: d1_B    , d2_BB
      double precision             :: d2_AB   , d2_AT  , d2_BT 
      double precision             :: d3_AAA  , d3_BBB , d3_TTT 
      double precision             :: d3_AAB  , d3_ABB , d3_AAT , d3_ATT
      double precision             :: d3_BBT  , d3_BTT , d3_ABT

      double precision             :: d4_ABTT    , d4_ABBT    , d4_AABT 
      double precision             :: d4_AABB 

      double precision             ::  st ,  ct 
      double precision             :: st2 , ct2 
      double precision             :: st3 , ct3 
      double precision             :: st4 , ct4 

      ! --------------------------------------------------------------- !

      o1 = pattern_id / 1000
      o2 = MOD(pattern_id / 100, 10)
      o3 = MOD(pattern_id / 10, 10)
      o4 = MOD(pattern_id, 10)

      ! --------------------------------------------------------------- !

      ax2 = ax  * ax
      ax3 = ax2 * ax
      ax4 = ax3 * ax


      select case (pattern_id)
        

      ! --------------------- only s function ------------------------- ! 

        case (0000) ! | s   s   s   s    ( 1 ) 

          call dqag(f0000, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


      ! --------------------- one  p function ------------------------- !

        case (0001) ! | s   s   s   px   ( 2 ) 
        
          call dqag(f0001, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

        case (0010) ! | s   s   px  s    ( 5 ) 
          call dqag(f0010, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

        case (0100) ! | s   px  s   s    ( 17)  
          call dqag(f0100, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

        case (1000) ! | px  s   s   s    ( 65) 
          call dqag(f1000, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

      ! --------------------- two  p function ------------------------- ! 

        case (0011) ! | s   s   px  px   ( 6 ) 
          call dqag(f0011, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        case (0101) ! | s   px  s   px   ( 18) 
          call dqag(f0101, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        case (0110) ! | s   px  px  s    ( 21)
          call dqag(f0110, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        case (1001) ! | px  s   s   px   ( 66)
          call dqag(f1001, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        case (1010) ! | px  s   px  s    ( 69)
          call dqag(f1010, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
        case (1100) ! | px  px  s   s    ( 81)
          call dqag(f1100, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

        ! --------------------- three  p function --------------------- !
        
        case (0111) ! | s   px  px  px   ( 22) 
            call dqag(f0111, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)

        case (1011) ! | px  s   px  px   ( 70) 
            call dqag(f1011, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)

        case (1101) ! | px  px  s   px   ( 82) 
              call dqag(f1101, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                       abserr, neval, ier, limit, lenw, last, &
                       iwork, work)

         case (1110) ! | px  px  px  s    ( 85) 
             call dqag(f1110, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                       abserr, neval, ier, limit, lenw, last, &
                       iwork, work)
          
        ! --------------------- four   p function --------------------- !

        case (1111) ! | px  px  px  px   ( 86)
            call dqag(f1111, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)




        











        
        

        case default 
          result = 0.d0 

      end select 









      contains 

      double precision function f0000(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)
      
      ! - derivative part - !

      integral_t  = int_xxxx
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0000

      double precision function f0001(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz1   = bessi_scaled(1,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      sqd  = dsin(xqd)
      cqd  = dcos(xqd)
      
      d1_B  = (A * ct + B) * bz1 * z / (z2+eta)
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)

      integral_t  = int_xxxx
      der_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
      f           = der_t * integral_t

      end function f0001


      double precision function f0010(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz1   = bessi_scaled(1,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      sqc  = dsin(xqc)
      cqc  = dcos(xqc)
      
      d1_B  = (A * ct + B) * bz1 * z / (z2+eta)
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)

      integral_t  = int_xxxx
      der_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
      f           = der_t * integral_t

      end function f0010

      double precision function f0100(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz1   = bessi_scaled(1,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      spb  = dsin(xpb)
      cpb  = dcos(xpb)
      
      d1_A  = (A + B * ct) * bz1 * z / (z2+eta)
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)

      integral_t  = int_xxxx
      der_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
      f           = der_t * integral_t

      end function f0100

      double precision function f1000(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz1   = bessi_scaled(1,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      spa  = dsin(xpa)
      cpa  = dcos(xpa)
      
      d1_A  = (A + B * ct) * bz1 * z / (z2+eta)
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)

      integral_t  = int_xxxx
      der_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
      f           = der_t * integral_t

      end function f1000


      double precision function f0011(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))
      
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)
      d2_BB =   ( 2.d0 * A2  *     st2      * bz1 + z * (A*ct + B) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (ccd*bz0 - cqcd*d2_BB - d1_T*sqcd*B/(B2*B+eta) + d2_BT*sqcd*B/(B*B+eta))/ax2 
      f           = der_t * integral_t

      end function f0011

      double precision function f0110(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      

      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqc = dsin(xqc)
      cqc = dcos(xqc)

      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))
      
      
      d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (A*B*d2_AB*spb*sqc - A*cqc*d2_AT*spb + B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax2)
      f           = der_t * integral_t

      end function f0110


      double precision function f1010(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqc = dsin(xqc)
      cqc = dcos(xqc)
      
      
      d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (A*B*d2_AB*spa*sqc - A*cqc*d2_AT*spa + B*cpa*d2_BT*sqc - cpa*cqc*d2_TT)/(A*B*ax2)
      f           = der_t * integral_t

      end function f1010

      
      double precision function f1100(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z
      z3 = z2 * z
      z4 = z3 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      cab  = dcos(xpa) *  dcos(xpb)
      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      
      d1_T  = -A * B * st  * bz1 * z / (z2+eta)
      d2_AA =   ( 2.d0 * B2  *     st2      * bz1 + z * (A + B*ct) * (A + B*ct) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (cab*bz0 - cpab*d2_AA - d1_T*spab*A/(A2*A+eta) + d2_AT*spab*A/(A*A+eta)) / ax2
      f           = der_t * integral_t

      end function f1100

      double precision function f1001(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqd = dsin(xqd)
      cqd = dcos(xqd)
      
      d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (A*B*d2_AB*spa*sqd - A*cqd*d2_AT*spa + B*cpa*d2_BT*sqd - cpa*cqd*d2_TT)/(A*B*ax**2)
      f           = der_t * integral_t

      end function f1001

      double precision function f0101(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqd = dsin(xqd)
      cqd = dcos(xqd)
      
      d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      integral_t  = int_xxxx
      der_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
      f           = der_t * integral_t

      end function f0101

      double precision function f0111(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      A2          = A  * A
      A3          = A2 * A
      A4          = A3 * A
      B           = 2.d0 * q_x / (ax2)
      B2          = B  * B
      B3          = B2 * B
      B4          = B3 * B 
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z
      z6 = z5 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)

      st  = dsin(theta)       ;    ct  = dcos(theta)
      st2 = st * st           ;    ct2 = ct * ct  
      
      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqd = dsin(xqd)
      cqd = dcos(xqd)

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      
      d1_A = (A + B * ct) * bz1 * z / (z2+eta)
      d1_T = -A * B * st  * bz1 * z / (z2+eta)

      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      d3_ABB = (-4.d0 * A * B * st2 * z * (A*ct + B) * (bz0 + bz2) + 4.d0 * A * st2 * (-A2 + A*B*ct + 2.d0*B2)*bz1 + z * (A + B * ct)*(2.d0 * A2 * st2 * (bz0 + bz2) + z * (A*ct + B)*(A*ct + B) * (3.d0 * bz1 + bz3))) * z / ( 4.d0 * (z6+eta))
      d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
      d3_BBT = -A * st * (2.d0 * A3 * ct * z * bz0 - 4.d0 * A3 * ct * bz1 + A2 * B * ct2 * z2 * bz1 + A2 * B * ct2 * z * bz0 - A2 * B * ct2 * bz1 + A2 * B * st2 * z * bz0 - A2 * B * st2 * bz1 + 2.d0 * A2 * B * z * bz0 - 5.d0 * A2 * B * bz1 + 2.d0 * A * B2 * ct * z2 * bz1 + B3 * z2 * bz1 - B3 * z * bz0 + 2.d0 * B3 * bz1) * z /(z6+eta)
      d3_BTT =  A * (-A4 * ct * bz1 - A3 * B * ct2 * z * bz0 + 2.d0 * A3 * B * st2 * z * bz0 - 2.d0 * A3 * B * st2 * bz1 - 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - 2.d0 * A * B3 * ct2 * z * bz0 + A * B3 * st2 * z2 * bz1 - A * B3 * st2 * z * bz0 + 2.d0 * A * B3 * st2 * bz1 - A * B3 * z * bz0 + 2.d0 * A * B3 * bz1 - B4 * ct * z * bz0 + B4 * ct * bz1) * z /(z6+eta)
      

      integral_t  = int_xxxx
      der_t       = (A*B*B2*spb*(ccd*d1_A - cqcd*d3_ABB) - A*B*d2_AT*spb*sqcd   & 
      &              + A*B2*d3_ABT*spb*sqcd + B*B2*cpb*(ccd*d1_T - cqcd*d3_BBT) &
                     - B*cpb*d2_TT*sqcd     + B2*cpb*d3_BTT*sqcd)*(A*B*B2)/((A*B*B2*A*B*B2+eta)*ax3)
      f           = der_t * integral_t

      end function f0111

      double precision function f1011(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      A2          = A  * A
      A3          = A2 * A
      A4          = A3 * A
      B           = 2.d0 * q_x / (ax2)
      B2          = B  * B
      B3          = B2 * B
      B4          = B3 * B 
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z
      z6 = z5 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)

      st  = dsin(theta)       ;    ct  = dcos(theta)
      st2 = st * st           ;    ct2 = ct * ct  
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqd = dsin(xqd)
      cqd = dcos(xqd)

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      
      d1_A = (A + B * ct) * bz1 * z / (z2+eta)
      d1_T = -A * B * st  * bz1 * z / (z2+eta)

      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      d3_ABB = (-4.d0 * A * B * st2 * z * (A*ct + B) * (bz0 + bz2) + 4.d0 * A * st2 * (-A2 + A*B*ct + 2.d0*B2)*bz1 + z * (A + B * ct)*(2.d0 * A2 * st2 * (bz0 + bz2) + z * (A*ct + B)*(A*ct + B) * (3.d0 * bz1 + bz3))) * z / ( 4.d0 * (z6+eta))
      d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
      d3_BBT = -A * st * (2.d0 * A3 * ct * z * bz0 - 4.d0 * A3 * ct * bz1 + A2 * B * ct2 * z2 * bz1 + A2 * B * ct2 * z * bz0 - A2 * B * ct2 * bz1 + A2 * B * st2 * z * bz0 - A2 * B * st2 * bz1 + 2.d0 * A2 * B * z * bz0 - 5.d0 * A2 * B * bz1 + 2.d0 * A * B2 * ct * z2 * bz1 + B3 * z2 * bz1 - B3 * z * bz0 + 2.d0 * B3 * bz1) * z /(z6+eta)
      d3_BTT =  A * (-A4 * ct * bz1 - A3 * B * ct2 * z * bz0 + 2.d0 * A3 * B * st2 * z * bz0 - 2.d0 * A3 * B * st2 * bz1 - 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - 2.d0 * A * B3 * ct2 * z * bz0 + A * B3 * st2 * z2 * bz1 - A * B3 * st2 * z * bz0 + 2.d0 * A * B3 * st2 * bz1 - A * B3 * z * bz0 + 2.d0 * A * B3 * bz1 - B4 * ct * z * bz0 + B4 * ct * bz1) * z /(z6+eta)
      

      integral_t  = int_xxxx
      der_t       = (A*B*B2*spa*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spa*sqcd &
      &             +A*B2*d3_ABT*spa*sqcd + 1.0*B*B2*cpa*(ccd*d1_T - cqcd*d3_BBT) &
      &             -B*cpa*d2_TT*sqcd + 1.0*B2*cpa*d3_BTT*sqcd)/(A*B*B2*ax3)
      f           = der_t * integral_t

      end function f1011


      double precision function f1101(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      A2          = A  * A
      A3          = A2 * A
      A4          = A3 * A
      B           = 2.d0 * q_x / (ax2)
      B2          = B  * B
      B3          = B2 * B
      B4          = B3 * B 
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z
      z6 = z5 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)

      st  = dsin(theta)       ;    ct  = dcos(theta)
      st2 = st * st           ;    ct2 = ct * ct  
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqd = dsin(xqd)
      cqd = dcos(xqd)

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      
      d1_B = (A * ct + B) * bz1 * z / (z2+eta)
      d1_T = -A * B * st  * bz1 * z / (z2+eta)

      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      d3_AAB = (-4.d0 * A * B * st2 * z * (A + B*ct)     *  (bz0 + bz2) + 2*B2*st2*z*(A*ct + B)*(bz0 + bz2) + 4*B*st2*(2*A2 + A*B*ct - B2)*bz1 + z2*(A + B*ct)**2*(A*ct + B)*(3*bz1 + bz3))*z/(4*(z6+eta))
      d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
      d3_AAT = - B * st * (A3 * z2 * bz1 - A3 * z * bz0 + 2.d0 * A3 * bz1 + 2.d0 * A2 * B * ct * z2 * bz1 + A * B2 * ct2 * z2 * bz1 + A * B2 * ct2 * z * bz0 - A * B2 * ct2 * bz1 + A * B2 * st2 * z * bz0 - A * B2 * st2 * bz1 + 2.d0 * A * B2 * z * bz0 - 5.d0 * A * B2 * bz1 + 2.d0 * B3 * ct * z * bz0 - 4.d0 * B3 * ct * bz1) * z /(z6+eta)
      d3_ATT = B * (-A4 * ct * z * bz0 + A4 * ct * bz1 - 2.d0 * A3 * B * ct2 * z * bz0 + A3 * B * st2 * z2 * bz1 - A3 * B * st2 * z * bz0 + 2.d0 * A3 * B * st2 * bz1 - A3 * B * z * bz0 + 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - A * B3 * ct2 * z * bz0 + 2.d0 * A * B3 * st2 * z * bz0 - 2.d0 * A * B3 * st2 * bz1 - 2.d0 * A * B3 * bz1 - B4 * ct * bz1) * z/(z6+eta)
      

      integral_t  = int_xxxx
      der_t       = (A*A2*B*sqd*(cab*d1_B - cpab*d3_AAB) + A*A2*cqd*(-cab*d1_T + cpab*d3_AAT) & 
      &             -A*B*d2_BT*spab*sqd + A*cqd*d2_TT*spab &
      &             +A2*B*d3_ABT*spab*sqd - A2*cqd*d3_ATT*spab)/(A*A2*B*ax3)
      f           = der_t * integral_t

      end function f1101

      double precision function f1110(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      A2          = A  * A
      A3          = A2 * A
      A4          = A3 * A
      B           = 2.d0 * q_x / (ax2)
      B2          = B  * B
      B3          = B2 * B
      B4          = B3 * B 
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z
      z6 = z5 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)

      st  = dsin(theta)       ;    ct  = dcos(theta)
      st2 = st * st           ;    ct2 = ct * ct  
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqc = dsin(xqc)
      cqc = dcos(xqc)

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      
      d1_B = (A * ct + B) * bz1 * z / (z2+eta)
      d1_T = -A * B * st  * bz1 * z / (z2+eta)

      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

      d3_AAB = (-4.d0 * A * B * st2 * z * (A + B*ct)     *  (bz0 + bz2) + 2*B2*st2*z*(A*ct + B)*(bz0 + bz2) + 4*B*st2*(2*A2 + A*B*ct - B2)*bz1 + z2*(A + B*ct)**2*(A*ct + B)*(3*bz1 + bz3))*z/(4*(z6+eta))
      d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
      d3_AAT = - B * st * (A3 * z2 * bz1 - A3 * z * bz0 + 2.d0 * A3 * bz1 + 2.d0 * A2 * B * ct * z2 * bz1 + A * B2 * ct2 * z2 * bz1 + A * B2 * ct2 * z * bz0 - A * B2 * ct2 * bz1 + A * B2 * st2 * z * bz0 - A * B2 * st2 * bz1 + 2.d0 * A * B2 * z * bz0 - 5.d0 * A * B2 * bz1 + 2.d0 * B3 * ct * z * bz0 - 4.d0 * B3 * ct * bz1) * z /(z6+eta)
      d3_ATT = B * (-A4 * ct * z * bz0 + A4 * ct * bz1 - 2.d0 * A3 * B * ct2 * z * bz0 + A3 * B * st2 * z2 * bz1 - A3 * B * st2 * z * bz0 + 2.d0 * A3 * B * st2 * bz1 - A3 * B * z * bz0 + 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - A * B3 * ct2 * z * bz0 + 2.d0 * A * B3 * st2 * z * bz0 - 2.d0 * A * B3 * st2 * bz1 - 2.d0 * A * B3 * bz1 - B4 * ct * bz1) * z/(z6+eta)
      

      integral_t  = int_xxxx
      der_t       = (A*A2*B*sqc*(cab*d1_B - cpab*d3_AAB) + A*A2*cqc*(-cab*d1_T + cpab*d3_AAT) &
      &             -A*B*d2_BT*spab*sqc + A*cqc*d2_TT*spab &
      &             +A2*B*d3_ABT*spab*sqc - A2*cqc*d3_ATT*spab)/(A*A2*B*ax3)
      f           = der_t * integral_t

      end function f1110

      double precision function f1111(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      A2          = A  * A
      A3          = A2 * A
      A4          = A3 * A
      A5          = A4 * A
      A6          = A5 * A 
      A7          = A6 * A
      B           = 2.d0 * q_x / (ax2)
      B2          = B  * B
      B3          = B2 * B
      B4          = B3 * B
      B5          = B4 * B
      B6          = B5 * B
      B7          = B6 * B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z
      z6 = z5 * z
      z7 = z6 * z 
      z8 = z7 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)        ;    ct  = dcos(theta)
      st2 = st  * st           ;    ct2 = ct  * ct
      st3 = st2 * st           ;    ct3 = ct2 * ct
      st4 = st3 * st           ;    ct4 = ct3 * ct
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqc = dsin(xqc)
      cqc = dcos(xqc)

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      
      d1_T = -A * B * st  * bz1 * z / (z2+eta)

      d2_AA =   ( 2.d0 * B2  *     st2      * bz1 + z * (A + B*ct) * (A + B*ct) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_BB =   ( 2.d0 * A2  *     st2      * bz1 + z * (A*ct + B) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
      d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)
      d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
      d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)

      d3_AAT = - B * st * (A3 * z2 * bz1 - A3 * z * bz0 + 2.d0 * A3 * bz1 + 2.d0 * A2 * B * ct * z2 * bz1 + A * B2 * ct2 * z2 * bz1 + A * B2 * ct2 * z * bz0 - A * B2 * ct2 * bz1 + A * B2 * st2 * z * bz0 - A * B2 * st2 * bz1 + 2.d0 * A * B2 * z * bz0 - 5.d0 * A * B2 * bz1 + 2.d0 * B3 * ct * z * bz0 - 4.d0 * B3 * ct * bz1) * z /(z6+eta)
      d3_ATT = B * (-A4 * ct * z * bz0 + A4 * ct * bz1 - 2.d0 * A3 * B * ct2 * z * bz0 + A3 * B * st2 * z2 * bz1 - A3 * B * st2 * z * bz0 + 2.d0 * A3 * B * st2 * bz1 - A3 * B * z * bz0 + 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - A * B3 * ct2 * z * bz0 + 2.d0 * A * B3 * st2 * z * bz0 - 2.d0 * A * B3 * st2 * bz1 - 2.d0 * A * B3 * bz1 - B4 * ct * bz1) * z/(z6+eta)
      d3_BBT = -A * st * (2.d0 * A3 * ct * z * bz0 - 4.d0 * A3 * ct * bz1 + A2 * B * ct2 * z2 * bz1 + A2 * B * ct2 * z * bz0 - A2 * B * ct2 * bz1 + A2 * B * st2 * z * bz0 - A2 * B * st2 * bz1 + 2.d0 * A2 * B * z * bz0 - 5.d0 * A2 * B * bz1 + 2.d0 * A * B2 * ct * z2 * bz1 + B3 * z2 * bz1 - B3 * z * bz0 + 2.d0 * B3 * bz1) * z /(z6+eta)
      d3_BTT =  A * (-A4 * ct * bz1 - A3 * B * ct2 * z * bz0 + 2.d0 * A3 * B * st2 * z * bz0 - 2.d0 * A3 * B * st2 * bz1 - 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - 2.d0 * A * B3 * ct2 * z * bz0 + A * B3 * st2 * z2 * bz1 - A * B3 * st2 * z * bz0 + 2.d0 * A * B3 * st2 * bz1 - A * B3 * z * bz0 + 2.d0 * A * B3 * bz1 - B4 * ct * z * bz0 + B4 * ct * bz1) * z /(z6+eta)
      
      d4_AABT = -st * (A5 * z2 * bz1 - A5 * z * bz0 + 2.d0 * A5 * bz1 + A4 * B * ct * z3 * bz0 + A4 * B * ct * z2 * bz1 + 4.d0 * A4 * B * ct * z * bz0 - 8.d0 * A4 * B * ct * bz1 + 2.d0 * A3 * B2 * ct2 * z3 * bz0 + A3 * B2 * ct2 * z2 * bz1 + 3.d0 * A3 * B2 * ct2 * z * bz0 - 5.d0 * A3 * B2 * ct2 * bz1 - 2 * A3 * B2 * st2 * z2 * bz1 + 3.d0 * A3 * B2 * st2 * z * bz0 - 5.d0 * A3 * B2 * st2 * bz1 + A3 * B2 * z3 * bz0 - 2.d0 * A3 * B2 * z2 * bz1 + 11.d0 * A3 * B2 * z * bz0 - 23.d0 * A3 * B2 * bz1 + A2 * B3 * ct3 * z3 * bz0 + A2 * B3 * ct3 * z2 * bz1 - A2 * B3 * ct3 * z * bz0 - A2 * B3 * ct * st2 * z2 * bz1 - A2 * B3 * ct * st2 * z * bz0 + 2.d0 * A2 * B3 * ct * z3 * bz0 + 9.d0 * A2 * B3 * ct * z * bz0 - 16.d0 * A2 * B3 * ct * bz1 + A * B4 * ct2 * z3 * bz0 + 2.d0 * A * B4 * ct2 * z2 * bz1 - 3.d0 * A * B4 * ct2 * bz1 + A * B4 * st2 * z2 * bz1 - 4.d0 * A * B4 * st2 * z * bz0 + 5.d0 * A * B4 * st2 * bz1 + 2.d0 * A * B4 * z2 * bz1 - 5.d0 * A * B4 * z * bz0 + 13.d0 * A * B4 * bz1 + 2.d0 * B5 * ct * z2 * bz1 - 4.d0 * B5 * ct * z * bz0 + 8.d0 * B5 * ct * bz1)*z/(z8+eta)
      d4_ABBT = st * ( -2.d0 * A7 * ct * bz1 - A6 * B * ct2 * z * bz0 - 5.d0 * A6 * B * ct2 * bz1 - 3.d0 * A6 * B * bz1 - 3.d0 * A5 * B2 * ct3 * z * bz0 - 4.d0 * A5 * B2 * ct3 * bz1 - 2.d0 * A5 * B2 * ct * z * bz0 - 7.d0 * A5 * B2 * ct * bz1 + 4.d0 * A5 * ct * z * bz0 - 8.d0 * A5 * ct * bz1 - 2.d0 * A4 * B3 * ct4 * z * bz0 - 4.d0 * A4 * B3 * ct4 * bz1 - 7.d0 * A4 * B3 * ct2 * z * bz0 - 2.d0 * A4 * B3 * ct2 * bz1 - A4 * B3 * z * bz0 + A4 * B3 * bz1 - 4.d0 * A4 * B * ct2 * z * bz0 + 8.d0 * A4 * B * ct2 * bz1 + 9.d0 * A4 * B * z * bz0 - 18.d0 * A4 * B * bz1 - 5.d0 * A3 * B4 * ct3 * z * bz0 - 8.d0 * A3 * B4 * ct3 * bz1 - 5.d0 * A3 * B4 * ct * z * bz0 + 8.d0 * A3 * B4 * ct * bz1 - 8.d0 * A3 * B2 * ct * z * bz0 + 16.d0 * A3 * B2 * ct * bz1 - 4.d0 * A2 * B5 * ct2 * z * bz0 - 5.d0 * A2 * B5 * ct2 * bz1 - A2 * B5 * z * bz0 + 3.d0 * A2 * B5 * bz1 - 14.d0 * A2 * B3 * z * bz0 + 28.d0 * A2 * B3 * bz1 - A * B6 * ct * z * bz0 - 3.d0 * A * B6 * ct * bz1 - 4.d0 * A * B4 * ct * z * bz0 + 8.d0  * A * B4 * ct * bz1 - B7 * bz1 + B5 * z * bz0 - 2.d0 * B5 * bz1)*z/(z8+eta)
      d4_AABB = (12.d0 * A2 * B2 * st4 * z * (bz0 + bz2) - 8.d0  * A * B * st2 * z2 * ( A + B * ct) * (A * ct + B) * ( 3.d0 * bz1 + bz3) + 8.d0 * A * st2 * z *(A + B*ct) * (bz0 + bz2) * (- A2 + A*B*ct + 2.d0 * B2) + 2.d0 * B2 * st2 * z2 * (A * ct + B)**2 * (3.d0 * bz1 + bz3) + 8.d0 * B * st2 * z * (A * ct + B) * (bz0 + bz2) * (2.d0 * A2 + A * B * ct - B2) + 8.d0 * st2 *(2.d0 * A4 - 4.d0 * A3 * B * ct + A2 * B2 * st2 - 12.d0 * A2 * B2 - 4.d0 * A * B3 * ct + 2.d0 * B4)*bz1 + z2 * (A + B * ct)**2 *(2.d0 * A2 * st2 * (bz1 + bz3) + 4.d0 * A2 * st2 * bz1 + 2.d0 * z * (A*ct + B)**2*(bz0 + bz2) + z * (A*ct + B)**2 * (bz0 + 2.d0*bz2 + bz4))) * z /( 8.d0 * (z8+eta))
      d4_ABTT = ( - A6 * ct * z * bz0 + A6 * ct * bz1 - A5 * B * ct2 * z2 * bz1 - A5 * B * ct2 * z * bz0 - 4.d0 * A5 * B * ct2 * bz1 + 2.d0 * A5 * B * st2 * z2 * bz1 - 2.d0 * A5 * B * st2 * z * bz0 + 4.d0 * A5 * B * st2 * bz1 - 2.d0 * A5 * B * z * bz0 + 4.d0 * A5 * B * bz1 - 2.d0 * A4 * B2 * ct3 * z2 * bz1 - 2.d0 * A4 * B2 * ct3 * z * bz0 - A4 * B2 * ct3 * bz1 + A4 * B2 * ct * st2 * z3 * bz0 + 2.d0 * A4 * B2 * ct * st2 * z2 * bz1 + 6.d0 * A4 * B2 * ct * st2 * z * bz0 - 5.d0 * A4 * B2 * ct * st2 * bz1 - 2.d0 * A4 * B2 * ct * z2 * bz1 - A4 * B2 * ct * z * bz0 - 8.d0 * A4 * B **2*ct*bz1 - A3 * B3 * ct4 * z2 * bz1 - A3 * B3 * ct4 * z * bz0 - A3 * B3 * ct4 * bz1 + A3 * B3 * ct2 * st2 * z3 * bz0 + 2.d0 * A3 * B3 * ct2 * st2 * z2 * bz1 + 2.d0 * A3 * B3 * ct2 * st2 * z * bz0 - A3 * B3 * ct2 * st2 * bz1 - 4.d0 * A3 * B3 * ct2 * z2 * bz1 - 4.d0 * A3 * B3 * ct2 * z * bz0 - A3 * B3 * st4 * z2 * bz1 - A3 * B3 * st4 * z * bz0 + A3 * B3 * st2 * z3 * bz0 - 2.d0 * A3 * B3 * st2 * z2 * bz1 + 14.d0 * A3 * B3 * st*2 * z * bz0 - 17.d0 * A3 * B3 * st2 * bz1 - A3 * B3 * z2 * bz1 + 3.d0 * A3 * B3 * z * bz0 - 15.d0 * A3 * B3 * bz1 - 2.d0 * A2 * B4 * ct3 * z2 * bz1 - 2.d0 * A2 * B4 * ct3 * z * bz0 - A2 * B4 * ct3 * bz1 + A2 * B4 * ct * st2 * z3 * bz0 + 2.d0 * A2 * B4 * ct * st2 * z2 * bz1 + 6.d0 * A2 * B4 * ct * st2 * z * bz0 - 5.d0 * A2 * B4 * ct * st2 * bz1 - 2.d0 * A2 * B4 * ct * z2 * bz1 - A2 * B4 * ct * z * bz0 - 8.d0 * A2 * B4 * ct * bz1 - A * B5 * ct2 * z2 * bz1 - A * B5 * ct2 * z * bz0 - 4.d0 * A * B5 * ct2 * bz1 + 2.d0 * A * B5 * st2 * z2 * bz1 - 2.d0 * A * B5 * st2 * z * bz0 + 4.d0 * A * B5 * st2 * bz1 - 2.d0 * A * B5 * z * bz0 + 4.d0 * A * B5 * bz1 - B6 * ct * z * bz0 + B6 * ct * bz1)*z/(z8+eta)

      integral_t  = int_xxxx
      der_t       = (A*A2*B*B2*(cab*ccd*bz0 - cab*cqcd*d2_BB - ccd*cpab*d2_AA + cpab*cqcd*d4_AABB) &
      &             +A*A2*B*sqcd*(-cab*d1_T + cpab*d3_AAT) + A*A2*B2*sqcd*(cab*d2_BT - cpab*d4_AABT) &
      &             +A*B*B2*spab*(-ccd*d1_T + cqcd*d3_BBT) + A*B*d2_TT*spab*sqcd &
      &             -A*B2*d3_BTT*spab*sqcd + A2*B*B2*spab*(ccd*d2_AT - cqcd*d4_ABBT) &
      &             -A2*B*d3_ATT*spab*sqcd + A2*B2*d4_ABTT*spab*sqcd)/(A*A2*B*B2*ax4)
      f           = der_t * integral_t

      end function f1111






end subroutine integrate_ERI_integral_mod