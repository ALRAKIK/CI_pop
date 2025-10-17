subroutine integrate_ERI_integral(pattern_id,px_count,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,result)

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

      double precision,parameter         :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier,last, neval , iwork(limit)
      double precision                   :: abserr, work(lenw)
      integer, parameter                 :: key  = 6
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi




      call dqag(f, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
            abserr, neval, ier, limit, lenw, last, &
            iwork, work)


      !call dqags(f, 0.d0, 2.d0*pi , Epsabs, Epsrel, Result, Abserr, Neval, Ier, &
      !               Limit, Lenw, Last, Iwork, Work)

      !if (ier /= 0) then
      !  write(*,'(A,I8,A)') 'Error code = ', ier
      !end if

      contains

      double precision function f(theta)

      use gsl_bessel_mod
      use torus_init

      implicit none

      double precision, intent(in) :: theta
      
      ! local !

      integer                      :: o1,o2,o3,o4 
      double precision             :: A, A2, A3, A4, A5, A6, A7 
      double precision             :: B, B2, B3, B4, B5, B6, B7 
      double precision             :: z, z2, z3, z4, z5, z6, z7, z8
      double precision             :: bz0, bz1, bz2, bz3, bz4 
      double precision             :: ax2, ax3, ax4
      double precision             :: term
      double precision             :: const
      double precision             :: erfcx
      double precision,parameter   :: eta = 1e-40

      double precision             :: spa   , cpa  , spb   , cpb
      double precision             :: sqc   , cqc  , sqd   , cqd  
      double precision             :: spab  , cpab , sqcd  , cqcd
      double precision             :: cab   , ccd
      double precision             :: psi 

      ! ------------------- integral variables ------------------------ !

      double precision             :: const_t , integral_t 
      double precision             :: int_xxxx , int_yyxx
      double precision             :: int_yxyx , int_xxyy
      double precision             :: int_yyyy , int_yyzz
      double precision             :: int_yzyz 
      double precision             :: beta, beta2, p2, q2, p3, q3 
      double precision             :: term1 , term2 ,term3
      double precision             :: inv
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


      ax2 = ax  * ax
      ax3 = ax2 * ax 
      ax4 = ax3 * ax 

      A = 2.d0 * p_x / (ax2)
      
      A2 = A  * A
      A3 = A2 * A 
      A4 = A3 * A 
      A5 = A4 * A 
      A6 = A5 * A
      A7 = A6 * A 

      B = 2.d0 * q_x / (ax2)

      B2 = B  * B 
      B3 = B2 * B 
      B4 = B3 * B 
      B5 = B4 * B  
      B6 = B5 * B
      B7 = B6 * B
      
      term  = dsqrt( dabs( 2.d0*p*q/(p+q)* (1.d0 - dcos( phi - theta ) ) ) ) / ax

      z     = dsqrt(  dabs ( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
  
      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z 
      z5 = z4 * z
      z6 = z5 * z 
      z7 = z6 * z 
      z8 = z7 * z 

      spa = dsin(xpa); spb = dsin(xpb); sqc = dsin(xqc); sqd  = dsin(xqd)
      cpa = dcos(xpa); cpb = dcos(xpb); cqc = dcos(xqc); cqd  = dcos(xqd)
      
      spab = dsin(ax*(2.d0*xp-xa-xb)) ; sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cpab = dcos(ax*(2.d0*xp-xa-xb)) ; cqcd = dcos(ax*(2.d0*xq-xc-xd))
      
      cab =  dcos(xpa) * dcos(xpb)    ; ccd  = dcos(xqc) *  dcos(xqd)

      st = dsin(theta)               
      st2 = st  * st                  
      st3 = st2 * st                 
      st4 = st3 * st                

      ct  = dcos(theta)
      ct2 = ct  * ct
      ct3 = ct2 * ct
      ct4 = ct3 * ct

      psi = phi - theta 

      const = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2)

      ! -------------------- general integral ------------------------- !

      beta     = 2.d0 * (1.d0-dcos(psi)) / ax2 
      beta2    = beta * beta


      q2 = q  * q 
      q3 = q2 * q
      p2 = p  * p 
      p3 = p2 * p


      ppq      = p+q 
      ppq2     = ppq  * ppq 
      ppq3     = ppq2 * ppq
      ppq4     = ppq3 * ppq

      inv      = ppq/(p*q)

      ! xxxx !

      int_xxxx = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)   

      ! xxyy ! xxzz ! 

      term1    = 2.d0 * p * ppq * dsqrt(beta)

      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))

      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / ( q * ppq3 )  

      ! yyxx ! zzxx ! 

      term1    = 2.d0 * p * q2 * dsqrt(beta)

      term2    = erfcx(term) * dsqrt(pi)  * (2.d0*p2 + q2 + p * q * (3.d0-2.d0*q*beta)) / dsqrt(inv)

      int_yyxx = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / (p2*q*ppq2) 

      ! yxyx ! zxzx ! 

      term1    = -2.d0 * dsqrt(pi*beta) * ppq

      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )

      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! yyyy ! zzzz ! 

      term1 = dexp(term*term) * pi * ( ppq2 * (4.d0 * p + q ) * ( p + 4.d0 * q ) - 4.d0 * p * ( p - 2.d0 * q) * ( 2.d0 * p - q ) * q * ppq * beta  + 12.d0 * p3 * q3 * beta2 )

      term2 = 2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * dsqrt(beta) * (4.d0 * p3 - 3.d0 * p * q2 + 4.d0 * q3 - 3.d0 * p2 * q * (1.d0+2.d0*q*beta) ) 

      term3 = dexp(term*term) * pi * ( - (p+q) * (p+q) * (4.d0*p+q) * (4.d0*q+p) + 4.d0 * p * (p-2.d0*q) * (2.d0*p-q) * q * (p+q) * beta - 12.d0 * p3 * q3 * beta2) * erf(term)

      int_yyyy = const * 0.015625d0 * ( term1 + term2 + term3 ) / ( p2 * dsqrt(inv) * q2 * ppq4 )

      ! yyzz ! zzyy !

      term1 = dexp(term*term) * pi * ( ppq2 * ( 4.d0 * p2 + 11.d0 * p * q + 4.d0 * q2 ) - 4.d0 * p * q * ppq * ( 2.d0 * p2 + p * q + 2.d0 * q2 ) * beta + 4.d0 * p3 * q3 * beta2 )

      term2 = 2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * dsqrt(beta) * ( 4.d0 * p3 + 7.d0 * p * q2 + 4.d0 * q3 + p2 * q * (7.d0 - 2.d0 * q * beta) )

      term3 = dexp(term*term) * pi * ( - ppq2 * ( 4.d0 * p2 + 11.d0 * p * q + 4.d0 * q2 ) + 4.d0 * p * q * ppq * ( 2.d0 * p2 + p * q + 2.d0 * q2 ) * beta - 4.d0 * p3 * q3 * beta2 ) * erf(term)

      int_yyzz = const * 0.015625d0 * (term1 + term2 + term3) / (p2 * dsqrt(inv) * q2 * ppq4)

      ! yzyz ! zyzy ! 

      term1 =   2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * sqrt(beta) * ( 5.d0 * ppq + 2.d0 * p * q * beta )

      term2 = - dexp( term * term ) * pi * (3.d0 * ppq2 + 12.d0 * p * q * ppq * beta + 4.d0 * p2 * q2 * beta2)
      
      term3 =   dexp( term * term ) * pi * (3.d0 * ppq2 + 12.d0 * p * q * ppq * beta + 4.d0 * p2 * q2 * beta2) * erf(term)

      int_yzyz =  - const * 0.015625d0 * (term1 + term2 + term3) / (p * dsqrt(inv) * q * ppq4)

      ! --------------------------------------------------------------- !

      o1 = pattern_id / 1000
      o2 = MOD(pattern_id / 100, 10)
      o3 = MOD(pattern_id / 10, 10)
      o4 = MOD(pattern_id, 10)

      ! -------------------------- Derivative ------------------------- !

      select case (px_count)
        case (0)

          ! Handle no px orbitals

        case (1) 

          ! Handle exactly one px orbital

          bz1 = bessi_scaled(1,z)
          
          d1_A = (A + B * ct) * bz1 * z / (z2+eta)
          d1_B = (A * ct + B) * bz1 * z / (z2+eta)
          d1_T = -A * B * st  * bz1 * z / (z2+eta)

        case (2)
          
          ! Handle two px orbitals

          bz0 = bessi_scaled(0,z)
          bz1 = bessi_scaled(1,z)
          bz2 = bessi_scaled(2,z)
          
          d1_A = (A + B * ct) * bz1 * z / (z2+eta)
          d1_B = (A * ct + B) * bz1 * z / (z2+eta)
          d1_T = -A * B * st  * bz1 * z / (z2+eta)

          d2_AA =   ( 2.d0 * B2  *     st2      * bz1 + z * (A + B*ct) * (A + B*ct) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_BB =   ( 2.d0 * A2  *     st2      * bz1 + z * (A*ct + B) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
          d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
          d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

        case (3)

          ! Handle three px orbitals

          bz0 = bessi_scaled(0,z)
          bz1 = bessi_scaled(1,z)
          bz2 = bessi_scaled(2,z)
          bz3 = bessi_scaled(3,z)

          d1_A = (A + B * ct) * bz1 * z / (z2+eta)
          d1_B = (A * ct + B) * bz1 * z / (z2+eta)
          d1_T = -A * B * st  * bz1 * z / (z2+eta)

          d2_AA =   ( 2.d0 * B2  *     st2      * bz1 + z * (A + B*ct) * (A + B*ct) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_BB =   ( 2.d0 * A2  *     st2      * bz1 + z * (A*ct + B) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
          d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
          d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

          d3_AAA = (A + B * ct ) * ( 6.d0   * B2 * st2 * z   *  (bz0 + bz2) - 12*B2*st2*bz1 + z2*(A + B*ct)**2*(3*bz1 + bz3))*z/(4*(z6+eta))
          d3_AAB = (-4.d0 * A * B * st2 * z * (A + B*ct)     *  (bz0 + bz2) + 2*B2*st2*z*(A*ct + B)*(bz0 + bz2) + 4*B*st2*(2*A2 + A*B*ct - B2)*bz1 + z2*(A + B*ct)**2*(A*ct + B)*(3*bz1 + bz3))*z/(4*(z6+eta))
          d3_AAT = - B * st * (A3 * z2 * bz1 - A3 * z * bz0 + 2.d0 * A3 * bz1 + 2.d0 * A2 * B * ct * z2 * bz1 + A * B2 * ct2 * z2 * bz1 + A * B2 * ct2 * z * bz0 - A * B2 * ct2 * bz1 + A * B2 * st2 * z * bz0 - A * B2 * st2 * bz1 + 2.d0 * A * B2 * z * bz0 - 5.d0 * A * B2 * bz1 + 2.d0 * B3 * ct * z * bz0 - 4.d0 * B3 * ct * bz1) * z /(z6+eta)
          d3_ABB = (-4.d0 * A * B * st2 * z * (A*ct + B) * (bz0 + bz2) + 4.d0 * A * st2 * (-A2 + A*B*ct + 2.d0*B2)*bz1 + z * (A + B * ct)*(2.d0 * A2 * st2 * (bz0 + bz2) + z * (A*ct + B)*(A*ct + B) * (3.d0 * bz1 + bz3))) * z / ( 4.d0 * (z6+eta))
          d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
          d3_ATT = B * (-A4 * ct * z * bz0 + A4 * ct * bz1 - 2.d0 * A3 * B * ct2 * z * bz0 + A3 * B * st2 * z2 * bz1 - A3 * B * st2 * z * bz0 + 2.d0 * A3 * B * st2 * bz1 - A3 * B * z * bz0 + 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - A * B3 * ct2 * z * bz0 + 2.d0 * A * B3 * st2 * z * bz0 - 2.d0 * A * B3 * st2 * bz1 - 2.d0 * A * B3 * bz1 - B4 * ct * bz1) * z/(z6+eta)
          d3_BBB = ( A * ct + B) * (6.d0 * A2 * st2 * z * (bz0 + bz2) - 12.d0 * A2 * st2 * bz1 + z2 * (A * ct + B)* (A * ct + B) * (3.d0 * bz1 + bz3)) * z /(4.d0*(z6+eta))
          d3_BBT = -A * st * (2.d0 * A3 * ct * z * bz0 - 4.d0 * A3 * ct * bz1 + A2 * B * ct2 * z2 * bz1 + A2 * B * ct2 * z * bz0 - A2 * B * ct2 * bz1 + A2 * B * st2 * z * bz0 - A2 * B * st2 * bz1 + 2.d0 * A2 * B * z * bz0 - 5.d0 * A2 * B * bz1 + 2.d0 * A * B2 * ct * z2 * bz1 + B3 * z2 * bz1 - B3 * z * bz0 + 2.d0 * B3 * bz1) * z /(z6+eta)
          d3_BTT =  A * (-A4 * ct * bz1 - A3 * B * ct2 * z * bz0 + 2.d0 * A3 * B * st2 * z * bz0 - 2.d0 * A3 * B * st2 * bz1 - 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - 2.d0 * A * B3 * ct2 * z * bz0 + A * B3 * st2 * z2 * bz1 - A * B3 * st2 * z * bz0 + 2.d0 * A * B3 * st2 * bz1 - A * B3 * z * bz0 + 2.d0 * A * B3 * bz1 - B4 * ct * z * bz0 + B4 * ct * bz1) * z /(z6+eta)
          d3_TTT = A * B * st * (-A2 * B2 * st2 * z2 * ( 3.d0 * bz1 + bz3) + 6.d0 * A * B * z * (bz0 + bz2) * (A2 * ct + A * B * ct2 + A * B + B2 * ct) + 4.d0 * ( A4 + A3 * B * ct - A2 * B2 * st2 + A * B3 * ct + B4 ) * bz1) * z/(4.d0 * (z6+eta))

          case (4)

          ! Handle four px orbitals

          bz0 = bessi_scaled(0,z)
          bz1 = bessi_scaled(1,z)
          bz2 = bessi_scaled(2,z)
          bz3 = bessi_scaled(3,z)
          bz4 = bessi_scaled(4,z)

          d1_A = (A + B * ct) * bz1 * z / (z2+eta)
          d1_B = (A * ct + B) * bz1 * z / (z2+eta)
          d1_T = -A * B * st  * bz1 * z / (z2+eta)

          d2_AA =   ( 2.d0 * B2  *     st2      * bz1 + z * (A + B*ct) * (A + B*ct) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_BB =   ( 2.d0 * A2  *     st2      * bz1 + z * (A*ct + B) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AB = - ( 2.d0 * A*B *     st2      * bz1 - z * (A + B*ct) * (A*ct + B) * (bz0 + bz2)) * z / (2.d0*(z4+eta))
          d2_AT =    B    * st   * ( - A2 * z   * bz0 + A2                          *  bz1 - A*B*ct*z  * bz0 - B2 *      bz1) * z / (z4+eta)
          d2_BT =    A    * st   * ( - A2       * bz1 - A*B*ct*z                    *  bz0 - B2*z      * bz0 + B2 *      bz1) * z / (z4+eta)
          d2_TT =    A    * B    * ( - A2 * ct  * bz1 + A*B*st2*z                   *  bz0 - 2.d0*A*B  * bz1 - B2 * ct * bz1) * z / (z4+eta)

          d3_AAA = (A + B * ct ) * ( 6.d0   * B2 * st2 * z   *  (bz0 + bz2) - 12*B2*st2*bz1 + z2*(A + B*ct)**2*(3*bz1 + bz3))*z/(4*(z6+eta))
          d3_AAB = (-4.d0 * A * B * st2 * z * (A + B*ct)     *  (bz0 + bz2) + 2*B2*st2*z*(A*ct + B)*(bz0 + bz2) + 4*B*st2*(2*A2 + A*B*ct - B2)*bz1 + z2*(A + B*ct)**2*(A*ct + B)*(3*bz1 + bz3))*z/(4*(z6+eta))
          d3_AAT = - B * st * (A3 * z2 * bz1 - A3 * z * bz0 + 2.d0 * A3 * bz1 + 2.d0 * A2 * B * ct * z2 * bz1 + A * B2 * ct2 * z2 * bz1 + A * B2 * ct2 * z * bz0 - A * B2 * ct2 * bz1 + A * B2 * st2 * z * bz0 - A * B2 * st2 * bz1 + 2.d0 * A * B2 * z * bz0 - 5.d0 * A * B2 * bz1 + 2.d0 * B3 * ct * z * bz0 - 4.d0 * B3 * ct * bz1) * z /(z6+eta)
          d3_ABB = (-4.d0 * A * B * st2 * z * (A*ct + B) * (bz0 + bz2) + 4.d0 * A * st2 * (-A2 + A*B*ct + 2.d0*B2)*bz1 + z * (A + B * ct)*(2.d0 * A2 * st2 * (bz0 + bz2) + z * (A*ct + B)*(A*ct + B) * (3.d0 * bz1 + bz3))) * z / ( 4.d0 * (z6+eta))
          d3_ABT = -st * (A4 * z * bz0 - A4 * bz1 + A3 * B * ct * z2 * bz1 + A3 * B * ct * z * bz0 + 2.d0 * A3 * B * ct * bz1 + A2 * B2 * ct2 * z2 * bz1 + A2 * B2 * ct2 * z * bz0 - A2 * B2 * st2 * z * bz0 + A2 * B2 * z2 * bz1 - A2 * B2 * z * bz0 + 6.d0 * A2 * B2 * bz1 + A * B3 * ct * z2 * bz1 + A * B3 * ct * z * bz0 + 2.d0 * A * B3 * ct * bz1 + B4 * z * bz0 - B4 * bz1) * z /(z6+eta)
          d3_ATT = B * (-A4 * ct * z * bz0 + A4 * ct * bz1 - 2.d0 * A3 * B * ct2 * z * bz0 + A3 * B * st2 * z2 * bz1 - A3 * B * st2 * z * bz0 + 2.d0 * A3 * B * st2 * bz1 - A3 * B * z * bz0 + 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - A * B3 * ct2 * z * bz0 + 2.d0 * A * B3 * st2 * z * bz0 - 2.d0 * A * B3 * st2 * bz1 - 2.d0 * A * B3 * bz1 - B4 * ct * bz1) * z/(z6+eta)
          d3_BBB = ( A * ct + B) * (6.d0 * A2 * st2 * z * (bz0 + bz2) - 12.d0 * A2 * st2 * bz1 + z2 * (A * ct + B)* (A * ct + B) * (3.d0 * bz1 + bz3)) * z /(4.d0*(z6+eta))
          d3_BBT = -A * st * (2.d0 * A3 * ct * z * bz0 - 4.d0 * A3 * ct * bz1 + A2 * B * ct2 * z2 * bz1 + A2 * B * ct2 * z * bz0 - A2 * B * ct2 * bz1 + A2 * B * st2 * z * bz0 - A2 * B * st2 * bz1 + 2.d0 * A2 * B * z * bz0 - 5.d0 * A2 * B * bz1 + 2.d0 * A * B2 * ct * z2 * bz1 + B3 * z2 * bz1 - B3 * z * bz0 + 2.d0 * B3 * bz1) * z /(z6+eta)
          d3_BTT =  A * (-A4 * ct * bz1 - A3 * B * ct2 * z * bz0 + 2.d0 * A3 * B * st2 * z * bz0 - 2.d0 * A3 * B * st2 * bz1 - 2.d0 * A3 * B * bz1 - A2 * B2 * ct3 * z * bz0 + A2 * B2 * ct * st2 * z2 * bz1 + A2 * B2 * ct * st2 * z * bz0 - 2.d0 * A2 * B2 * ct * z * bz0 - 2.d0 * A * B3 * ct2 * z * bz0 + A * B3 * st2 * z2 * bz1 - A * B3 * st2 * z * bz0 + 2.d0 * A * B3 * st2 * bz1 - A * B3 * z * bz0 + 2.d0 * A * B3 * bz1 - B4 * ct * z * bz0 + B4 * ct * bz1) * z /(z6+eta)
          d3_TTT = A * B * st * (-A2 * B2 * st2 * z2 * ( 3.d0 * bz1 + bz3) + 6.d0 * A * B * z * (bz0 + bz2) * (A2 * ct + A * B * ct2 + A * B + B2 * ct) + 4.d0 * ( A4 + A3 * B * ct - A2 * B2 * st2 + A * B3 * ct + B4 ) * bz1) * z/(4.d0 * (z6+eta))

          d4_ABTT = ( - A6 * ct * z * bz0 + A6 * ct * bz1 - A5 * B * ct2 * z2 * bz1 - A5 * B * ct2 * z * bz0 - 4.d0 * A5 * B * ct2 * bz1 + 2.d0 * A5 * B * st2 * z2 * bz1 - 2.d0 * A5 * B * st2 * z * bz0 + 4.d0 * A5 * B * st2 * bz1 - 2.d0 * A5 * B * z * bz0 + 4.d0 * A5 * B * bz1 - 2.d0 * A4 * B2 * ct3 * z2 * bz1 - 2.d0 * A4 * B2 * ct3 * z * bz0 - A4 * B2 * ct3 * bz1 + A4 * B2 * ct * st2 * z3 * bz0 + 2.d0 * A4 * B2 * ct * st2 * z2 * bz1 + 6.d0 * A4 * B2 * ct * st2 * z * bz0 - 5.d0 * A4 * B2 * ct * st2 * bz1 - 2.d0 * A4 * B2 * ct * z2 * bz1 - A4 * B2 * ct * z * bz0 - 8.d0 * A4 * B **2*ct*bz1 - A3 * B3 * ct4 * z2 * bz1 - A3 * B3 * ct4 * z * bz0 - A3 * B3 * ct4 * bz1 + A3 * B3 * ct2 * st2 * z3 * bz0 + 2.d0 * A3 * B3 * ct2 * st2 * z2 * bz1 + 2.d0 * A3 * B3 * ct2 * st2 * z * bz0 - A3 * B3 * ct2 * st2 * bz1 - 4.d0 * A3 * B3 * ct2 * z2 * bz1 - 4.d0 * A3 * B3 * ct2 * z * bz0 - A3 * B3 * st4 * z2 * bz1 - A3 * B3 * st4 * z * bz0 + A3 * B3 * st2 * z3 * bz0 - 2.d0 * A3 * B3 * st2 * z2 * bz1 + 14.d0 * A3 * B3 * st*2 * z * bz0 - 17.d0 * A3 * B3 * st2 * bz1 - A3 * B3 * z2 * bz1 + 3.d0 * A3 * B3 * z * bz0 - 15.d0 * A3 * B3 * bz1 - 2.d0 * A2 * B4 * ct3 * z2 * bz1 - 2.d0 * A2 * B4 * ct3 * z * bz0 - A2 * B4 * ct3 * bz1 + A2 * B4 * ct * st2 * z3 * bz0 + 2.d0 * A2 * B4 * ct * st2 * z2 * bz1 + 6.d0 * A2 * B4 * ct * st2 * z * bz0 - 5.d0 * A2 * B4 * ct * st2 * bz1 - 2.d0 * A2 * B4 * ct * z2 * bz1 - A2 * B4 * ct * z * bz0 - 8.d0 * A2 * B4 * ct * bz1 - A * B5 * ct2 * z2 * bz1 - A * B5 * ct2 * z * bz0 - 4.d0 * A * B5 * ct2 * bz1 + 2.d0 * A * B5 * st2 * z2 * bz1 - 2.d0 * A * B5 * st2 * z * bz0 + 4.d0 * A * B5 * st2 * bz1 - 2.d0 * A * B5 * z * bz0 + 4.d0 * A * B5 * bz1 - B6 * ct * z * bz0 + B6 * ct * bz1)*z/(z8+eta)
          d4_AABT = -st * (A5 * z2 * bz1 - A5 * z * bz0 + 2.d0 * A5 * bz1 + A4 * B * ct * z3 * bz0 + A4 * B * ct * z2 * bz1 + 4.d0 * A4 * B * ct * z * bz0 - 8.d0 * A4 * B * ct * bz1 + 2.d0 * A3 * B2 * ct2 * z3 * bz0 + A3 * B2 * ct2 * z2 * bz1 + 3.d0 * A3 * B2 * ct2 * z * bz0 - 5.d0 * A3 * B2 * ct2 * bz1 - 2 * A3 * B2 * st2 * z2 * bz1 + 3.d0 * A3 * B2 * st2 * z * bz0 - 5.d0 * A3 * B2 * st2 * bz1 + A3 * B2 * z3 * bz0 - 2.d0 * A3 * B2 * z2 * bz1 + 11.d0 * A3 * B2 * z * bz0 - 23.d0 * A3 * B2 * bz1 + A2 * B3 * ct3 * z3 * bz0 + A2 * B3 * ct3 * z2 * bz1 - A2 * B3 * ct3 * z * bz0 - A2 * B3 * ct * st2 * z2 * bz1 - A2 * B3 * ct * st2 * z * bz0 + 2.d0 * A2 * B3 * ct * z3 * bz0 + 9.d0 * A2 * B3 * ct * z * bz0 - 16.d0 * A2 * B3 * ct * bz1 + A * B4 * ct2 * z3 * bz0 + 2.d0 * A * B4 * ct2 * z2 * bz1 - 3.d0 * A * B4 * ct2 * bz1 + A * B4 * st2 * z2 * bz1 - 4.d0 * A * B4 * st2 * z * bz0 + 5.d0 * A * B4 * st2 * bz1 + 2.d0 * A * B4 * z2 * bz1 - 5.d0 * A * B4 * z * bz0 + 13.d0 * A * B4 * bz1 + 2.d0 * B5 * ct * z2 * bz1 - 4.d0 * B5 * ct * z * bz0 + 8.d0 * B5 * ct * bz1)*z/(z8+eta)
          d4_ABBT = st * ( -2.d0 * A7 * ct * bz1 - A6 * B * ct2 * z * bz0 - 5.d0 * A6 * B * ct2 * bz1 - 3.d0 * A6 * B * bz1 - 3.d0 * A5 * B2 * ct3 * z * bz0 - 4.d0 * A5 * B2 * ct3 * bz1 - 2.d0 * A5 * B2 * ct * z * bz0 - 7.d0 * A5 * B2 * ct * bz1 + 4.d0 * A5 * ct * z * bz0 - 8.d0 * A5 * ct * bz1 - 2.d0 * A4 * B3 * ct4 * z * bz0 - 4.d0 * A4 * B3 * ct4 * bz1 - 7.d0 * A4 * B3 * ct2 * z * bz0 - 2.d0 * A4 * B3 * ct2 * bz1 - A4 * B3 * z * bz0 + A4 * B3 * bz1 - 4.d0 * A4 * B * ct2 * z * bz0 + 8.d0 * A4 * B * ct2 * bz1 + 9.d0 * A4 * B * z * bz0 - 18.d0 * A4 * B * bz1 - 5.d0 * A3 * B4 * ct3 * z * bz0 - 8.d0 * A3 * B4 * ct3 * bz1 - 5.d0 * A3 * B4 * ct * z * bz0 + 8.d0 * A3 * B4 * ct * bz1 - 8.d0 * A3 * B2 * ct * z * bz0 + 16.d0 * A3 * B2 * ct * bz1 - 4.d0 * A2 * B5 * ct2 * z * bz0 - 5.d0 * A2 * B5 * ct2 * bz1 - A2 * B5 * z * bz0 + 3.d0 * A2 * B5 * bz1 - 14.d0 * A2 * B3 * z * bz0 + 28.d0 * A2 * B3 * bz1 - A * B6 * ct * z * bz0 - 3.d0 * A * B6 * ct * bz1 - 4.d0 * A * B4 * ct * z * bz0 + 8.d0  * A * B4 * ct * bz1 - B7 * bz1 + B5 * z * bz0 - 2.d0 * B5 * bz1)*z/(z8+eta)
          d4_AABB = (12.d0 * A2 * B2 * st4 * z * (bz0 + bz2) - 8.d0  * A * B * st2 * z2 * ( A + B * ct) * (A * ct + B) * ( 3.d0 * bz1 + bz3) + 8.d0 * A * st2 * z *(A + B*ct) * (bz0 + bz2) * (- A2 + A*B*ct + 2.d0 * B2) + 2.d0 * B2 * st2 * z2 * (A * ct + B)**2 * (3.d0 * bz1 + bz3) + 8.d0 * B * st2 * z * (A * ct + B) * (bz0 + bz2) * (2.d0 * A2 + A * B * ct - B2) + 8.d0 * st2 *(2.d0 * A4 - 4.d0 * A3 * B * ct + A2 * B2 * st2 - 12.d0 * A2 * B2 - 4.d0 * A * B3 * ct + 2.d0 * B4)*bz1 + z2 * (A + B * ct)**2 *(2.d0 * A2 * st2 * (bz1 + bz3) + 4.d0 * A2 * st2 * bz1 + 2.d0 * z * (A*ct + B)**2*(bz0 + bz2) + z * (A*ct + B)**2 * (bz0 + 2.d0*bz2 + bz4))) * z /( 8.d0 * (z8+eta))


      end select

      ! ----------------------- End Derivative ------------------------ !  
                      
      select case (pattern_id)
        
         case (0000) ! | s   s   s   s    ( 1 ) 
           const_t          = int_xxxx
           integral_t       = bessi_scaled(0, z)
           f                = const_t * integral_t

        case (0001) ! | s   s   s   px   ( 2 ) 
          const_t          = int_xxxx
          integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
          f                = const_t * integral_t
      
      case (0002) ! | s   s   s   py   ( 3 ) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0003) ! | s   s   s   pz   ( 4 ) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
        case (0010) ! | s   s   px  s    ( 5 ) 
         const_t          = int_xxxx
         integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
         f                = const_t * integral_t
      
        case (0011) ! | s   s   px  px   ( 6 ) 
         const_t          = int_xxxx
         integral_t       = (ccd*bessi_scaled(0, z) - cqcd*d2_BB - d1_T*sqcd*B/(B2*B+eta) + d2_BT*sqcd*B/(B*B+eta))/ax2 
         f                = const_t * integral_t 

      case (0012) ! | s   s   px  py   ( 7 ) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0013) ! | s   s   px  pz   ( 8 ) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0020) ! | s   s   py  s    ( 9 ) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0021) ! | s   s   py  px   ( 10) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0022) ! | s   s   py  py   ( 11) 
       const_t          = int_xxyy
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0023) ! | s   s   py  pz   ( 12) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0030) ! | s   s   pz  s    ( 13) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0031) ! | s   s   pz  px   ( 14) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0032) ! | s   s   pz  py   ( 15) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0033) ! | s   s   pz  pz   ( 16) 
        const_t          = int_xxyy
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
        case (0100) ! | s   px  s   s    ( 17) 
          const_t          = int_xxxx
          integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
          f                = const_t * integral_t
      
        case (0101) ! | s   px  s   px   ( 18) 
          const_t          = int_xxxx
          integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
          f                = const_t * integral_t
      
      case (0102) ! | s   px  s   py   ( 19) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (0103) ! | s   px  s   pz   ( 20) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
        case (0110) ! | s   px  px  s    ( 21) 
          const_t          = int_xxxx
          integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
          f                = const_t * integral_t
      
        case (0111) ! | s   px  px  px   ( 22) 
          const_t          = int_xxxx
          integral_t       = (A*B*B2*spb*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spb*sqcd + 1.0*A*B2*d3_ABT*spb*sqcd + 1.0*B*B2*cpb*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpb*d2_TT*sqcd + 1.0*B2*cpb*d3_BTT*sqcd)*(A*B*B2)/((A*B*B2*A*B*B2+eta)*ax*ax2)
          f                = const_t * integral_t
      
      case (0112) ! | s   px  px  py   ( 23) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (0113) ! | s   px  px  pz   ( 24) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (0120) ! | s   px  py  s    ( 25) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (0121) ! | s   px  py  px   ( 26) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
       case (0122) ! | s   px  py  py   ( 27) 
         const_t          = int_xxyy
         integral_t       = (A*d1_A*spb + cpb*d1_T)/(A*ax)
         f                = const_t * integral_t
      
      case (0123) ! | s   px  py  pz   ( 28) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (0130) ! | s   px  pz  s    ( 29) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (0131) ! | s   px  pz  px   ( 30) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (0132) ! | s   px  pz  py   ( 31) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (0133) ! | s   px  pz  pz   ( 32) 
        const_t          = int_xxyy
        integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (0200) ! | s   py  s   s    ( 33) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0201) ! | s   py  s   px   ( 34) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
       case (0202) ! | s   py  s   py   ( 35) 
         const_t          = int_yxyx
         integral_t       = bessi_scaled(0, z)
         f                = const_t * integral_t
      
      case (0203) ! | s   py  s   pz   ( 36) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0210) ! | s   py  px  s    ( 37) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0211) ! | s   py  px  px   ( 38) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
      case (0212) ! | s   py  px  py   ( 39) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (0213) ! | s   py  px  pz   ( 40) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0220) ! | s   py  py  s    ( 41) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (0221) ! | s   py  py  px   ( 42) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (0222) ! | s   py  py  py   ( 43) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0223) ! | s   py  py  pz   ( 44) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0230) ! | s   py  pz  s    ( 45) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0231) ! | s   py  pz  px   ( 46) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0232) ! | s   py  pz  py   ( 47) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0233) ! | s   py  pz  pz   ( 48) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0300) ! | s   pz  s   s    ( 49) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0301) ! | s   pz  s   px   ( 50) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0302) ! | s   pz  s   py   ( 51) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
       case (0303) ! | s   pz  s   pz   ( 52) 
         const_t          = int_yxyx
         integral_t       = bessi_scaled(0, z)
         f                = const_t * integral_t
      
      case (0310) ! | s   pz  px  s    ( 53) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0311) ! | s   pz  px  px   ( 54) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
      case (0312) ! | s   pz  px  py   ( 55) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0313) ! | s   pz  px  pz   ( 56) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (0320) ! | s   pz  py  s    ( 57) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0321) ! | s   pz  py  px   ( 58) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (0322) ! | s   pz  py  py   ( 59) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0323) ! | s   pz  py  pz   ( 60) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0330) ! | s   pz  pz  s    ( 61) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (0331) ! | s   pz  pz  px   ( 62) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (0332) ! | s   pz  pz  py   ( 63) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (0333) ! | s   pz  pz  pz   ( 64) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
        case (1000) ! | px  s   s   s    ( 65) 
          const_t          = int_xxxx
          integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
          f                = const_t * integral_t
      
        case (1001) ! | px  s   s   px   ( 66) 
          const_t          = int_xxxx
          integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
          f                = const_t * integral_t
      
      case (1002) ! | px  s   s   py   ( 67) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1003) ! | px  s   s   pz   ( 68) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
        case (1010) ! | px  s   px  s    ( 69) 
          const_t          = int_xxxx
          integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
          f                = const_t * integral_t
      
       case (1011) ! | px  s   px  px   ( 70) 
         const_t          = int_xxxx
         integral_t       = (A*B*B2*spa*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spa*sqcd + 1.0*A*B2*d3_ABT*spa*sqcd + 1.0*B*B2*cpa*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpa*d2_TT*sqcd + 1.0*B2*cpa*d3_BTT*sqcd)/(A*B*B2*ax*ax2)
         f                = const_t * integral_t
      
      case (1012) ! | px  s   px  py   ( 71) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1013) ! | px  s   px  pz   ( 72) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1020) ! | px  s   py  s    ( 73) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1021) ! | px  s   py  px   ( 74) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
       case (1022) ! | px  s   py  py   ( 75) 
         const_t          = int_xxyy
         integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
         f                = const_t * integral_t
      
      case (1023) ! | px  s   py  pz   ( 76) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1030) ! | px  s   pz  s    ( 77) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1031) ! | px  s   pz  px   ( 78) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1032) ! | px  s   pz  py   ( 79) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1033) ! | px  s   pz  pz   ( 80) 
        const_t          = int_xxyy
        integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
        f                = const_t * integral_t
      
        case (1100) ! | px  px  s   s    ( 81) 
          const_t          = int_xxxx
          integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - d1_T*spab*A/((A2*A+eta)*ax2) + d2_AT*spab*A/((A*A+eta)*ax2)
          f                = const_t * integral_t
      
       case (1101) ! | px  px  s   px   ( 82) 
         const_t          = int_xxxx
         integral_t       = (A*A2*B*sqd*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqd*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqd + 1.0*A*cqd*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqd - 1.0*A2*cqd*d3_ATT*spab)/(A*A2*B*ax*ax2)
         f                = const_t * integral_t
      
      case (1102) ! | px  px  s   py   ( 83) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
      case (1103) ! | px  px  s   pz   ( 84) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
       case (1110) ! | px  px  px  s    ( 85) 
         const_t          = int_xxxx
         integral_t       = (A*A2*B*sqc*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqc*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqc + 1.0*A*cqc*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqc - 1.0*A2*cqc*d3_ATT*spab)/(A*A2*B*ax*ax2)
         f                = const_t * integral_t
      
      case (1111) ! | px  px  px  px   ( 86) 
        const_t          = int_xxxx
        integral_t       = (A*A2*B*B2*(cab*ccd*bessi_scaled(0, z) - cab*cqcd*d2_BB - ccd*cpab*d2_AA + cpab*cqcd*d4_AABB) + 1.0*A*A2*B*sqcd*(-cab*d1_T + cpab*d3_AAT) + 1.0*A*A2*B2*sqcd*(cab*d2_BT - cpab*d4_AABT) + 1.0*A*B*B2*spab*(-ccd*d1_T + cqcd*d3_BBT) + 1.0*A*B*d2_TT*spab*sqcd - 1.0*A*B2*d3_BTT*spab*sqcd + 1.0*A2*B*B2*spab*(ccd*d2_AT - cqcd*d4_ABBT) - 1.0*A2*B*d3_ATT*spab*sqcd + 1.0*A2*B2*d4_ABTT*spab*sqcd)/(A*A2*B*B2*ax2**2)
        f                = const_t * integral_t

      case (1112) ! | px  px  px  py   ( 87) 
       const_t          = 0
       integral_t       = (A*A2*B*sqc*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqc*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqc + 1.0*A*cqc*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqc - 1.0*A2*cqc*d3_ATT*spab)/(A*A2*B*ax*ax2)
       f                = const_t * integral_t
      
      case (1113) ! | px  px  px  pz   ( 88) 
       const_t          = 0
       integral_t       = (A*A2*B*sqc*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqc*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqc + 1.0*A*cqc*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqc - 1.0*A2*cqc*d3_ATT*spab)/(A*A2*B*ax*ax2)
       f                = const_t * integral_t
      
      case (1120) ! | px  px  py  s    ( 89) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
      case (1121) ! | px  px  py  px   ( 90) 
       const_t          = 0
       integral_t       = (A*A2*B*sqd*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqd*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqd + 1.0*A*cqd*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqd - 1.0*A2*cqd*d3_ATT*spab)/(A*A2*B*ax*ax2)
       f                = const_t * integral_t
      
      case (1122) ! | px  px  py  py   ( 91) 
        const_t          = int_xxyy
        integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab*A/((A2*A+eta)*ax2) + 1.0*d2_AT*spab*A/((A*A+eta)*ax2)
        f                = const_t * integral_t
      
      case (1123) ! | px  px  py  pz   ( 92) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
      case (1130) ! | px  px  pz  s    ( 93) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
      case (1131) ! | px  px  pz  px   ( 94) 
       const_t          = 0
       integral_t       = (A*A2*B*sqd*(cab*d1_B - cpab*d3_AAB) + 1.0*A*A2*cqd*(-cab*d1_T + cpab*d3_AAT) - 1.0*A*B*d2_BT*spab*sqd + 1.0*A*cqd*d2_TT*spab + 1.0*A2*B*d3_ABT*spab*sqd - 1.0*A2*cqd*d3_ATT*spab)/(A*A2*B*ax*ax2)
       f                = const_t * integral_t
      
      case (1132) ! | px  px  pz  py   ( 95) 
       const_t          = 0
       integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab/(A2*ax2) + 1.0*d2_AT*spab/(A*ax2)
       f                = const_t * integral_t
      
      case (1133) ! | px  px  pz  pz   ( 96) 
        const_t          = int_xxyy
        integral_t       = cab*bessi_scaled(0, z)/ax2 - cpab*d2_AA/ax2 - 1.0*d1_T*spab*A/((A2*A+eta)*ax2) + 1.0*d2_AT*spab*A/((A*A+eta)*ax2)
        f                = const_t * integral_t
      
      case (1200) ! | px  py  s   s    ( 97) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1201) ! | px  py  s   px   ( 98) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1202) ! | px  py  s   py   ( 99) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (1203) ! | px  py  s   pz   ( 100) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1210) ! | px  py  px  s    ( 101) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1211) ! | px  py  px  px   ( 102) 
       const_t          = 0
       integral_t       = (A*B*B2*spa*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spa*sqcd + 1.0*A*B2*d3_ABT*spa*sqcd + 1.0*B*B2*cpa*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpa*d2_TT*sqcd + 1.0*B2*cpa*d3_BTT*sqcd)/(A*B*B2*ax*ax2)
       f                = const_t * integral_t
      
      case (1212) ! | px  py  px  py   ( 103) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (1213) ! | px  py  px  pz   ( 104) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1220) ! | px  py  py  s    ( 105) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (1221) ! | px  py  py  px   ( 106) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (1222) ! | px  py  py  py   ( 107) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1223) ! | px  py  py  pz   ( 108) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1230) ! | px  py  pz  s    ( 109) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1231) ! | px  py  pz  px   ( 110) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1232) ! | px  py  pz  py   ( 111) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1233) ! | px  py  pz  pz   ( 112) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1300) ! | px  pz  s   s    ( 113) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1301) ! | px  pz  s   px   ( 114) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1302) ! | px  pz  s   py   ( 115) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1303) ! | px  pz  s   pz   ( 116) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (1310) ! | px  pz  px  s    ( 117) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1311) ! | px  pz  px  px   ( 118) 
       const_t          = 0
       integral_t       = (A*B*B2*spa*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spa*sqcd + 1.0*A*B2*d3_ABT*spa*sqcd + 1.0*B*B2*cpa*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpa*d2_TT*sqcd + 1.0*B2*cpa*d3_BTT*sqcd)/(A*B*B2*ax*ax2)
       f                = const_t * integral_t
      
      case (1312) ! | px  pz  px  py   ( 119) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1313) ! | px  pz  px  pz   ( 120) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spa*sqc - 1.0*A*cqc*d2_AT*spa + 1.0*B*cpa*d2_BT*sqc - 1.0*cpa*cqc*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (1320) ! | px  pz  py  s    ( 121) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1321) ! | px  pz  py  px   ( 122) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (1322) ! | px  pz  py  py   ( 123) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1323) ! | px  pz  py  pz   ( 124) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1330) ! | px  pz  pz  s    ( 125) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (1331) ! | px  pz  pz  px   ( 126) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spa*sqd - 1.0*A*cqd*d2_AT*spa + 1.0*B*cpa*d2_BT*sqd - 1.0*cpa*cqd*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (1332) ! | px  pz  pz  py   ( 127) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (1333) ! | px  pz  pz  pz   ( 128) 
       const_t          = 0
       integral_t       = (A*d1_A*spa + 1.0*cpa*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2000) ! | py  s   s   s    ( 129) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2001) ! | py  s   s   px   ( 130) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2002) ! | py  s   s   py   ( 131) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2003) ! | py  s   s   pz   ( 132) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2010) ! | py  s   px  s    ( 133) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2011) ! | py  s   px  px   ( 134) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
       case (2012) ! | py  s   px  py   ( 135) 
         const_t          = int_yxyx
         integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
         f                = const_t * integral_t
      
      case (2013) ! | py  s   px  pz   ( 136) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2020) ! | py  s   py  s    ( 137) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2021) ! | py  s   py  px   ( 138) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (2022) ! | py  s   py  py   ( 139) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2023) ! | py  s   py  pz   ( 140) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2030) ! | py  s   pz  s    ( 141) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2031) ! | py  s   pz  px   ( 142) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2032) ! | py  s   pz  py   ( 143) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2033) ! | py  s   pz  pz   ( 144) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2100) ! | py  px  s   s    ( 145) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2101) ! | py  px  s   px   ( 146) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (2102) ! | py  px  s   py   ( 147) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (2103) ! | py  px  s   pz   ( 148) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2110) ! | py  px  px  s    ( 149) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (2111) ! | py  px  px  px   ( 150) 
       const_t          = 0
       integral_t       = (A*B*B2*spb*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spb*sqcd + 1.0*A*B2*d3_ABT*spb*sqcd + 1.0*B*B2*cpb*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpb*d2_TT*sqcd + 1.0*B2*cpb*d3_BTT*sqcd)/(A*B*B2*ax*ax2)
       f                = const_t * integral_t
      
      case (2112) ! | py  px  px  py   ( 151) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (2113) ! | py  px  px  pz   ( 152) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (2120) ! | py  px  py  s    ( 153) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (2121) ! | py  px  py  px   ( 154) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (2122) ! | py  px  py  py   ( 155) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2123) ! | py  px  py  pz   ( 156) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2130) ! | py  px  pz  s    ( 157) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2131) ! | py  px  pz  px   ( 158) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (2132) ! | py  px  pz  py   ( 159) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2133) ! | py  px  pz  pz   ( 160) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (2200) ! | py  py  s   s    ( 161) 
        const_t          = int_yyxx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2201) ! | py  py  s   px   ( 162) 
        const_t          = int_yyxx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (2202) ! | py  py  s   py   ( 163) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2203) ! | py  py  s   pz   ( 164) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2210) ! | py  py  px  s    ( 165) 
        const_t          = int_yyxx
        integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (2211) ! | py  py  px  px   ( 166) 
        const_t          = int_yyxx
        integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd*B/(B2*(B+eta)*ax2) + 1.0*d2_BT*sqcd*B/(B*(B+eta)*ax2)
        f                = const_t * integral_t
      
      case (2212) ! | py  py  px  py   ( 167) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2213) ! | py  py  px  pz   ( 168) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2220) ! | py  py  py  s    ( 169) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2221) ! | py  py  py  px   ( 170) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2222) ! | py  py  py  py   ( 171) 
        const_t          = int_yyyy
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2223) ! | py  py  py  pz   ( 172) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2230) ! | py  py  pz  s    ( 173) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2231) ! | py  py  pz  px   ( 174) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2232) ! | py  py  pz  py   ( 175) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2233) ! | py  py  pz  pz   ( 176) 
        const_t          = int_yyzz
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2300) ! | py  pz  s   s    ( 177) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2301) ! | py  pz  s   px   ( 178) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2302) ! | py  pz  s   py   ( 179) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2303) ! | py  pz  s   pz   ( 180) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2310) ! | py  pz  px  s    ( 181) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2311) ! | py  pz  px  px   ( 182) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
      case (2312) ! | py  pz  px  py   ( 183) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2313) ! | py  pz  px  pz   ( 184) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2320) ! | py  pz  py  s    ( 185) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2321) ! | py  pz  py  px   ( 186) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2322) ! | py  pz  py  py   ( 187) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
       case (2323) ! | py  pz  py  pz   ( 188) 
         const_t          = int_yzyz
         integral_t       = bessi_scaled(0, z)
         f                = const_t * integral_t
      
      case (2330) ! | py  pz  pz  s    ( 189) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (2331) ! | py  pz  pz  px   ( 190) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (2332) ! | py  pz  pz  py   ( 191) 
        const_t          = int_yzyz
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (2333) ! | py  pz  pz  pz   ( 192) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3000) ! | pz  s   s   s    ( 193) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3001) ! | pz  s   s   px   ( 194) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3002) ! | pz  s   s   py   ( 195) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3003) ! | pz  s   s   pz   ( 196) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3010) ! | pz  s   px  s    ( 197) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3011) ! | pz  s   px  px   ( 198) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
      case (3012) ! | pz  s   px  py   ( 199) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3013) ! | pz  s   px  pz   ( 200) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (3020) ! | pz  s   py  s    ( 201) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3021) ! | pz  s   py  px   ( 202) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3022) ! | pz  s   py  py   ( 203) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3023) ! | pz  s   py  pz   ( 204) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3030) ! | pz  s   pz  s    ( 205) 
        const_t          = int_yxyx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3031) ! | pz  s   pz  px   ( 206) 
        const_t          = int_yxyx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (3032) ! | pz  s   pz  py   ( 207) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3033) ! | pz  s   pz  pz   ( 208) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3100) ! | pz  px  s   s    ( 209) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3101) ! | pz  px  s   px   ( 210) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (3102) ! | pz  px  s   py   ( 211) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3103) ! | pz  px  s   pz   ( 212) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (3110) ! | pz  px  px  s    ( 213) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (3111) ! | pz  px  px  px   ( 214) 
       const_t          = 0
       integral_t       = (A*B*B2*spb*(ccd*d1_A - cqcd*d3_ABB) - 1.0*A*B*d2_AT*spb*sqcd + 1.0*A*B2*d3_ABT*spb*sqcd + 1.0*B*B2*cpb*(ccd*d1_T - cqcd*d3_BBT) - 1.0*B*cpb*d2_TT*sqcd + 1.0*B2*cpb*d3_BTT*sqcd)/(A*B*B2*ax*ax2)
       f                = const_t * integral_t
      
      case (3112) ! | pz  px  px  py   ( 215) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (3113) ! | pz  px  px  pz   ( 216) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spb*sqc - 1.0*A*cqc*d2_AT*spb + 1.0*B*cpb*d2_BT*sqc - 1.0*cpb*cqc*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (3120) ! | pz  px  py  s    ( 217) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3121) ! | pz  px  py  px   ( 218) 
       const_t          = 0
       integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
       f                = const_t * integral_t
      
      case (3122) ! | pz  px  py  py   ( 219) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3123) ! | pz  px  py  pz   ( 220) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3130) ! | pz  px  pz  s    ( 221) 
        const_t          = int_yxyx
        integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
        f                = const_t * integral_t
      
      case (3131) ! | pz  px  pz  px   ( 222) 
        const_t          = int_yxyx
        integral_t       = (A*B*d2_AB*spb*sqd - 1.0*A*cqd*d2_AT*spb + 1.0*B*cpb*d2_BT*sqd - 1.0*cpb*cqd*d2_TT)/(A*B*ax**2)
        f                = const_t * integral_t
      
      case (3132) ! | pz  px  pz  py   ( 223) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3133) ! | pz  px  pz  pz   ( 224) 
       const_t          = 0
       integral_t       = (A*d1_A*spb + 1.0*cpb*d1_T)/(A*ax)
       f                = const_t * integral_t
      
      case (3200) ! | pz  py  s   s    ( 225) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3201) ! | pz  py  s   px   ( 226) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3202) ! | pz  py  s   py   ( 227) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3203) ! | pz  py  s   pz   ( 228) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3210) ! | pz  py  px  s    ( 229) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3211) ! | pz  py  px  px   ( 230) 
       const_t          = 0
       integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
       f                = const_t * integral_t
      
      case (3212) ! | pz  py  px  py   ( 231) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3213) ! | pz  py  px  pz   ( 232) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3220) ! | pz  py  py  s    ( 233) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3221) ! | pz  py  py  px   ( 234) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3222) ! | pz  py  py  py   ( 235) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3223) ! | pz  py  py  pz   ( 236) 
        const_t          = int_yzyz
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3230) ! | pz  py  pz  s    ( 237) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3231) ! | pz  py  pz  px   ( 238) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3232) ! | pz  py  pz  py   ( 239) 
        const_t          = int_yzyz
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3233) ! | pz  py  pz  pz   ( 240) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3300) ! | pz  pz  s   s    ( 241) 
        const_t          = int_yyxx
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3301) ! | pz  pz  s   px   ( 242) 
        const_t          = int_yyxx
        integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (3302) ! | pz  pz  s   py   ( 243) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3303) ! | pz  pz  s   pz   ( 244) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3310) ! | pz  pz  px  s    ( 245) 
        const_t          = int_yyxx
        integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
        f                = const_t * integral_t
      
      case (3311) ! | pz  pz  px  px   ( 246) 
        const_t          = int_yyxx
        integral_t       = ccd*bessi_scaled(0, z)/ax2 - cqcd*d2_BB/ax2 - 1.0*d1_T*sqcd/(B2*ax2) + 1.0*d2_BT*sqcd/(B*ax2)
        f                = const_t * integral_t
      
      case (3312) ! | pz  pz  px  py   ( 247) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3313) ! | pz  pz  px  pz   ( 248) 
       const_t          = 0
       integral_t       = (B*d1_B*sqc - 1.0*cqc*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3320) ! | pz  pz  py  s    ( 249) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3321) ! | pz  pz  py  px   ( 250) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3322) ! | pz  pz  py  py   ( 251) 
        const_t          = int_yyzz 
        integral_t       = bessi_scaled(0, z)
        f                = const_t * integral_t
      
      case (3323) ! | pz  pz  py  pz   ( 252) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3330) ! | pz  pz  pz  s    ( 253) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
      case (3331) ! | pz  pz  pz  px   ( 254) 
       const_t          = 0
       integral_t       = (B*d1_B*sqd - 1.0*cqd*d1_T)/(B*ax)
       f                = const_t * integral_t
      
      case (3332) ! | pz  pz  pz  py   ( 255) 
       const_t          = 0
       integral_t       = bessi_scaled(0, z)
       f                = const_t * integral_t
      
       case (3333) ! | pz  pz  pz  pz   ( 256) 
         const_t          = int_yyyy
         integral_t       = bessi_scaled(0, z)
         f                = const_t * integral_t

       case default
         f = 0.0d0
       end select

       end function f

end subroutine integrate_ERI_integral