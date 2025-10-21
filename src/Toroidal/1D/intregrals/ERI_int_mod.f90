subroutine integrate_ERI_integral_mod(pattern_id,p,q,p_x,q_x,phi,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,result)

      use quadpack , only : dqag , dqags, dqagp
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod
      use tools 
      use, intrinsic :: ieee_arithmetic
      use constants_module

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
      
      double precision                   :: A, A2, A3, A4, A5, A6, A7 
      double precision                   :: B, B2, B3, B4, B5, B6, B7 
      double precision                   :: z, z2, z3, z4, z5, z6, z7
      double precision                   :: bz0, bz1, bz2, bz3, bz4 
      double precision                   :: bz5, bz6, bz7, bz8
      double precision                   :: ax2, ax3, ax4
      double precision                   :: term
      double precision                   :: const

      double precision                   :: spa   , cpa  , spb   , cpb
      double precision                   :: sqc   , cqc  , sqd   , cqd  
      double precision                   :: spab  , cpab , sqcd  , cqcd
      double precision                   :: cab   , ccd
      double precision                   :: psi 

      ! ------------------- integral variables ------------------------ !

      double precision                   :: der_t , integral_t 
      double precision                   :: int_xxxx , int_yyxx
      double precision                   :: int_yxyx , int_xxyy
      double precision                   :: int_yyyy , int_yyzz
      double precision                   :: int_yzyz , int_xyxy
      double precision                   :: beta, beta2, p2, q2, p3, q3 
      double precision                   :: term1 , term2
      double precision                   :: inv   , Iinv 
      double precision                   :: ppq , ppq2, ppq3 , ppq4 

      ! ------------------- derivative variables ---------------------- !

      double precision                   :: d1_T    , d2_TT
      double precision                   :: d1_A    , d2_AA 
      double precision                   :: d1_B    , d2_BB
      double precision                   :: d2_AB   , d2_AT  , d2_BT 
      double precision                   :: d3_AAB  , d3_ABB , d3_AAT , d3_ATT
      double precision                   :: d3_BBT  , d3_BTT , d3_ABT

      double precision                   :: d4_ABTT    , d4_ABBT    , d4_AABT 
      double precision                   :: d4_AABB 

      double precision                   ::  st ,  ct 
      double precision                   :: st2 , ct2 
      double precision                   :: st3 , ct3 
      double precision                   :: st4 , ct4 

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !

      ax2 = ax  * ax
      ax3 = ax2 * ax
      ax4 = ax3 * ax

      sf = 15.d0 

      select case (pattern_id)
        
      ! --------------------- only s function ------------------------- ! 

        case (0000) ! | s   s   s   s    ( 1 ) 

        call dqag(f0000, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                  abserr, neval, ier, limit, lenw, last, &
                  iwork, work)


         if (ier > 2) then
           write(*,'(A,I4)') 'Error code from the case ', pattern_id
           stop 
         end if
        
      ! --------------------- one  p function ------------------------- !

        case (0001) ! | s   s   s   px   ( 2 ) 

          call dqag(f0001, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                  abserr, neval, ier, limit, lenw, last, &
                  iwork, work)


          if (ier > 2) then
            write(*,'(A,I4)') 'Error code from the case ', pattern_id
            stop 
          end if

        case (0010) ! | s   s   px  s    ( 5 ) 

          call dqag(f0010, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

          if (ier > 2) then
            write(*,'(A,I4)') 'Error code from the case ', pattern_id
            stop 
          end if

        case (0100) ! | s   px  s   s    ( 17)  

          call dqag(f0100, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


          if (ier > 2) then
            write(*,'(A,I4)') 'Error code from the case ', pattern_id
            stop 
          end if
         
        case (1000) ! | px  s   s   s    ( 65) 

          call dqag(f1000, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

          if (ier > 2) then
            write(*,'(A,I4)') 'Error code from the case ', pattern_id
            stop 
          end if

      ! --------------------- two  p function ------------------------- ! 

        case (0011) ! | s   s   px  px   ( 6 ) 

          call dqag(f0011, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0101) ! | s   px  s   px   ( 18) 

          call dqag(f0101, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)
           

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0110) ! | s   px  px  s    ( 21)

          call dqag(f0110, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1001) ! | px  s   s   px   ( 66)

          call dqag(f1001, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1010) ! | px  s   px  s    ( 69)

          call dqag(f1010, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1100) ! | px  px  s   s    ( 81)

          call dqag(f1100, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

      ! ---------------------- three  p function ---------------------- !
        
        case (0111) ! | s   px  px  px   ( 22) 

           call dqag(f0111, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                     abserr, neval, ier, limit, lenw, last, &
                     iwork, work)


             if (ier > 2) then
               write(*,'(A,I4)') 'Error code from the case ', pattern_id
               stop 
             end if

        case (1011) ! | px  s   px  px   ( 70) 

            call dqag(f1011, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)

              if (ier > 2) then
                write(*,'(A,I4)') 'Error code from the case ', pattern_id
                stop 
              end if

        case (1101) ! | px  px  s   px   ( 82) 

            call dqag(f1101, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                       abserr, neval, ier, limit, lenw, last, &
                       iwork, work)
              

              if (ier > 2) then
                write(*,'(A,I4)') 'Error code from the case ', pattern_id
                stop 
              end if

        case (1110) ! | px  px  px  s    ( 85) 

            call dqag(f1110, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                       abserr, neval, ier, limit, lenw, last, &
                       iwork, work)


              if (ier > 2) then
                write(*,'(A,I4)') 'Error code from the case ', pattern_id
                stop 
              end if
          
      ! ---------------------- four   p function ---------------------- !

        case (1111) ! | px  px  px  px   ( 86)

          call dqag(f1111, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

      ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!! Y Y !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      
        case (0022,0033) ! | s   s   py  py   ( 11)

          call dqag(f0022, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)
            
            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0122,0133) ! | s   px  py  py   ( 27) 

          call dqag(f0122, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0202,0303) ! | s   py  s   py   ( 35) 

          call dqag(f0202, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0212,0313) ! | s   py  px  py   ( 39)
        
          call dqag(f0212, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0220,0330) ! | s   py  py  s    ( 41)

          call dqag(f0220, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

             

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (0221,0331) ! | s   py  py  px   ( 42)

          call dqag(f0221, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1022,1033) ! | px  s   py  py   ( 75)

          call dqag(f1022, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1122,1133) ! | px  px  py  py   ( 91)

          call dqag(f1122, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1202,1303) ! | px  py  s   py   ( 99)
          call dqag(f1202, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1212,1313) ! | px  py  px  py   ( 103)

          call dqag(f1212, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1220,1330) ! | px  py  py  s    ( 105)

          call dqag(f1220, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (1221,1331) ! | px  py  py  px   ( 106)
          
          call dqag(f1221, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

         case (2002,3003) ! | py  s   s   py   ( 131 , 196 ) 
           
          call dqag(f2002, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                     abserr, neval, ier, limit, lenw, last, &
                     iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2012,3013) ! | py  s   px  py   ( 135,200)
          
          call dqag(f2012, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2020,3030) ! | py  s   py  s    ( 137,205)

          call dqag(f2020, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2021,3031) ! | py  s   py  px   ( 138)
          
          call dqag(f2021, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if


        case (2102,3103) ! | py  px  s   py   ( 147) 
            
          call dqag(f2102, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                  abserr, neval, ier, limit, lenw, last, &
                  iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2112,3113) ! | py  px  px  py   ( 151) 

          call dqag(f2112, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2120,3130) ! | py  px  py  s    ( 153)

          call dqag(f2120, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)


            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2121,3131) ! | py  px  py  px   ( 154)

          call dqag(f2121, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2200,3300) ! | py  py  s   s    ( 161)
          
          call dqag(f2200, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if
      
        case (2201,3301) ! | py  py  s   px   ( 162)

          call dqag(f2201, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                     abserr, neval, ier, limit, lenw, last, &
                     iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2210,3310) ! | py  py  px  s    ( 165) 

         call dqag(f2210, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

        case (2211,3311) ! | py  py  px  px   ( 166)

          call dqag(f2211, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                    abserr, neval, ier, limit, lenw, last, &
                    iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

          case (2222,3333) ! | py  py  py  py   ( 171)

          call dqag(f2222, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                     abserr, neval, ier, limit, lenw, last, &
                     iwork, work)

            if (ier > 2) then
              write(*,'(A,I4)') 'Error code from the case ', pattern_id
              stop 
            end if

          case (2233,3322) ! | py  py  pz  pz   ( 176)

            call dqag(f2233, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)

              if (ier > 2) then
                write(*,'(A,I4)') 'Error code from the case ', pattern_id
                stop 
              end if
            
          case (2323,3232) ! | py  pz  py  pz   ( 188)

            call dqag(f2323, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                      abserr, neval, ier, limit, lenw, last, &
                      iwork, work)

              if (ier > 2) then
                write(*,'(A,I4)') 'Error code from the case ', pattern_id
                stop 
              end if

          case (2332,3223) ! | py  pz  pz  py   ( 191)

              call dqag(f2332, 0.d0, 2.d0*pi, epsabs, epsrel, key, result, &
                        abserr, neval, ier, limit, lenw, last, &
                        iwork, work)

                if (ier > 2) then
                  write(*,'(A,I4)') 'Error code from the case ', pattern_id
                  stop 
                end if

        case default 
          result = 0.d0 

      end select 

      result = result * dexp(-sf)


          
















      contains 

      double precision function f0000(theta) Result(f)

      double precision,intent(in)  :: theta

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      bz0   = bessi_scaled(0,z)
      bz2   = bessi_scaled(2,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      sqd  = dsin(xqd)
      cqd  = dcos(xqd)

      d1_B  =  0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = -0.5d0 *  A * B * st  * (bz0-bz2)


      integral_t  = int_xxxx

      if (dabs(B) < 1.d-15) then 
        der_t       = (d1_B*sqd)/ax
      else 
        der_t       = (d1_B*sqd - cqd*d1_T/B)/ax
      end if 

      f           =  der_t * integral_t

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz2   = bessi_scaled(2,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      sqc  = dsin(xqc)
      cqc  = dcos(xqc)

      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_xxxx

      if (dabs(B) < 1.d-15) then 
        der_t       = (d1_B*sqc)/ax
      else 
        der_t       = (d1_B*sqc - cqc*d1_T/B)/ax
      end if 

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz2   = bessi_scaled(2,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      spb  = dsin(xpb)
      cpb  = dcos(xpb)


      d1_A  =   0.5d0 * (A + ct * B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)
      
      integral_t  = int_xxxx

      if (dabs(A) < 1.d-15) then 
        der_t = (d1_A*spb)/ax
      else 
        der_t = (d1_A*spb + cpb*d1_T/A)/ax
      end if 
    
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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz2   = bessi_scaled(2,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      spa  = dsin(xpa)
      cpa  = dcos(xpa)

      d1_A  =   0.5d0 * (A + ct * B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)
      
      integral_t  = int_xxxx

      if (dabs(A) < 1.d-15) then 
        der_t       = (d1_A*spa)/ax
      else 
        der_t       = (d1_A*spa + cpa*d1_T/A)/ax
      end if 

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      
      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      d1_T  = - 0.5d0 * A * B * st  * (bz0-bz2)
      d2_BB = (3*A**2*bz0*ct**2 - 4*A**2*bz2*ct**2 + A**2*bz4*ct**2 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24

      integral_t  = int_xxxx
    
      if (dabs(B) < 1.d-15) then
        der_t = (ccd*bz0 - cqcd*d2_BB)/ ax2 
      else 
        der_t = (ccd*bz0 - cqcd*d2_BB + d1_T*sqcd/B2 - d2_BT*sqcd/B)/ ax2
      end if

        f     = der_t * integral_t

      end function f0011

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqd = dsin(xqd)
      cqd = dcos(xqd)

      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
          der_t = d2_AB*spb*sqd/ax2
        else if (dabs(A) < 1.d-15) then
          der_t = d2_AB*spb*sqd/ax2 - cqd*d2_AT*spb/(B*ax2)
        else if (dabs(B) < 1.d-15) then
          der_t = d2_AB*spb*sqd/ax2 + cpb*d2_BT*sqd/(A*ax2)
      endif
      
      else
        der_t = d2_AB*spb*sqd/ax2 - cqd*d2_AT*spb/(B*ax2) + cpb*d2_BT*sqd/(A*ax2) - cpb*cqd*d2_TT/(A*B*ax2)
      endif

      f           = der_t * integral_t

      end function f0101

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      
      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqc = dsin(xqc)
      cqc = dcos(xqc)

      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_xxxx

       if (dabs(A*B) < 1.d-30) then 
         if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
           der_t    = d2_AB*spb*sqc/ax2
         else if (dabs(A) < 1.d-15) then
           der_t    = (d2_AB*spb*sqc - cqc*d2_AT*spb/B)/ax2
         else if (dabs(B) < 1.d-15) then
           der_t    = (d2_AB*spb*sqc + cpb*d2_BT*sqc/A)/ax2
         end if
       else 
       der_t       = (d2_AB*spb*sqc - cqc*d2_AT*spb/B + cpb*d2_BT*sqc/A - cpb*cqc*d2_TT/(A*B))/ax2
       end if 

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqc = dsin(xqc)
      cqc = dcos(xqc)
      
      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
            der_t   = d2_AB*spa*sqc/ax2
        else if (dabs(A) < 1.d-15) then
            der_t   = d2_AB*spa*sqc/ax2 - cqc*d2_AT*spa/(B*ax2)
        else if (dabs(B) < 1.d-15) then
            der_t   = d2_AB*spa*sqc/ax2 + cpa*d2_BT*sqc/(A*ax2)
        endif
      else
        der_t       = d2_AB*spa*sqc/ax2 - cqc*d2_AT*spa/(B*ax2) + cpa*d2_BT*sqc/(A*ax2) - cpa*cqc*d2_TT/(A*B*ax2)
      end if

      f           = der_t * integral_t

      end function f1010

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqd = dsin(xqd)
      cqd = dcos(xqd)

      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
          der_t   = d2_AB*spa*sqd/ax2
        else if (dabs(A) < 1.d-15) then
          der_t   = d2_AB*spa*sqd/ax2 - cqd*d2_AT*spa/(B*ax2)
        else if (dabs(B) < 1.d-15) then
          der_t   = d2_AB*spa*sqd/ax2 + cpa*d2_BT*sqd/(A*ax2)
        endif
      else
        der_t     = d2_AB*spa*sqd/ax2 - cqd*d2_AT*spa/(B*ax2) + cpa*d2_BT*sqd/(A*ax2) - cpa*cqd*d2_TT/(A*B*ax2)
      end if 

      f    = der_t * integral_t

      end function f1001

      
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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)

      st  = dsin(theta)       ;    ct = dcos(theta)
      
      cab  = dcos(xpa) *  dcos(xpb)
      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))


      d1_T  = - 0.5d0 * A * B * st  * (bz0-bz2)
      d2_AA = (3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0*ct**2 - 4*B**2*bz2*ct**2 + B**2*bz4*ct**2 + 12*bz0 - 12*bz2)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24

      integral_t  = int_xxxx
      
      if (dabs(A) < 1.d-15) then
        der_t = (cab*bz0 - cpab*d2_AA) / ax2
      else
        der_t = (cab*bz0 - cpab*d2_AA - d1_T*spab/A2 + d2_AT*spab/A) / ax2
      endif

      f           = der_t * integral_t

      end function f1100

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
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
      bz4   = bessi_scaled(4,z)
      bz5   = bessi_scaled(5,z)
      bz6   = bessi_scaled(6,z)


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


      d1_A   =  0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T   = -0.5d0 *  A * B * st  * (bz0-bz2)
      d2_AT  = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_TT  =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24
      d3_ABB =  (10*A**3*bz0*ct**2 - 15*A**3*bz2*ct**2 + 6*A**3*bz4*ct**2 - A**3*bz6*ct**2 + 10*A**2*B*bz0*ct**3 + 20*A**2*B*bz0*ct - 15*A**2*B*bz2*ct**3 - 30*A**2*B*bz2*ct + 6*A**2*B*bz4*ct**3 + 12*A**2*B*bz4*ct - A**2*B*bz6*ct**3 - 2*A**2*B*bz6*ct + 20*A*B**2*bz0*ct**2 + 10*A*B**2*bz0 - 30*A*B**2*bz2*ct**2 - 15*A*B**2*bz2 + 12*A*B**2*bz4*ct**2 + 6*A*B**2*bz4 - 2*A*B**2*bz6*ct**2 - A*B**2*bz6 + 120*A*bz0*ct**2 + 60*A*bz0 - 160*A*bz2*ct**2 - 80*A*bz2 + 40*A*bz4*ct**2 + 20*A*bz4 + 10*B**3*bz0*ct - 15*B**3*bz2*ct + 6*B**3*bz4*ct - B**3*bz6*ct + 180*B*bz0*ct - 240*B*bz2*ct + 60*B*bz4*ct)/480
      d3_ABT =  -st*(10*A**3*B*bz0*ct - 15*A**3*B*bz2*ct + 6*A**3*B*bz4*ct - A**3*B*bz6*ct + 10*A**2*B**2*bz0*ct**2 + 10*A**2*B**2*bz0 - 15*A**2*B**2*bz2*ct**2 - 15*A**2*B**2*bz2 + 6*A**2*B**2*bz4*ct**2 + 6*A**2*B**2*bz4 - A**2*B**2*bz6*ct**2 - A**2*B**2*bz6 + 60*A**2*bz0 - 80*A**2*bz2 + 20*A**2*bz4 + 10*A*B**3*bz0*ct - 15*A*B**3*bz2*ct + 6*A*B**3*bz4*ct - A*B**3*bz6*ct + 180*A*B*bz0*ct - 240*A*B*bz2*ct + 60*A*B*bz4*ct + 60*B**2*bz0 - 80*B**2*bz2 + 20*B**2*bz4 + 240*bz0 - 240*bz2)/480
      d3_BBT =  -A*st*(10*A**2*B*bz0*ct**2 - 15*A**2*B*bz2*ct**2 + 6*A**2*B*bz4*ct**2 - A**2*B*bz6*ct**2 + 20*A*B**2*bz0*ct - 30*A*B**2*bz2*ct + 12*A*B**2*bz4*ct - 2*A*B**2*bz6*ct + 120*A*bz0*ct - 160*A*bz2*ct + 40*A*bz4*ct + 10*B**3*bz0 - 15*B**3*bz2 + 6*B**3*bz4 - B**3*bz6 + 180*B*bz0 - 240*B*bz2 + 60*B*bz4)/480
      d3_BTT =  A*(10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 + 10*A*B**3*bz0*st**2 - 15*A*B**3*bz2*st**2 + 6*A*B**3*bz4*st**2 - A*B**3*bz6*st**2 - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 60*B**2*bz0*ct + 80*B**2*bz2*ct - 20*B**2*bz4*ct - 240*bz0*ct + 240*bz2*ct)/480

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30 ) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
        der_t       = spb *  (ccd*d1_A - cqcd*d3_ABB)     /(ax3)   
        else if (dabs(A) < 1.d-15) then 
        der_t       = spb *  (ccd*d1_A - cqcd*d3_ABB)     /(ax3)      &
      &             + spb *  d2_AT*sqcd                   /(B2*ax3)   & 
      &             - spb *  d3_ABT*sqcd                  /(B*ax3)    
        else if (dabs(B) < 1.d-15) then 
        der_t       = spb *  (ccd*d1_A - cqcd*d3_ABB)     /(ax3)      &
      &             + cpb *  (ccd*d1_T - cqcd*d3_BBT)     /(A*ax3) 
        end if 
      else 
      der_t       =   spb *  (ccd*d1_A - cqcd*d3_ABB)     /(ax3)      &
      &             + spb *  d2_AT*sqcd                   /(B2*ax3)   & 
      &             - spb *  d3_ABT*sqcd                  /(B*ax3)    &
      &             + cpb *  (ccd*d1_T - cqcd*d3_BBT)     /(A*ax3)    & 
      &             + cpb *  d2_TT*sqcd                   /(A*B2*ax3) &
      &             - cpb *  d3_BTT*sqcd                  /(A*B*ax3)
      end if 

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
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
      bz4   = bessi_scaled(4,z)
      bz5   = bessi_scaled(5,z)
      bz6   = bessi_scaled(6,z)

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

      d1_A   =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T   = - 0.5d0 *  A * B * st  * (bz0-bz2)

      d2_AT  = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_TT  =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24 
      d3_ABB =  (10*A**3*bz0*ct**2 - 15*A**3*bz2*ct**2 + 6*A**3*bz4*ct**2 - A**3*bz6*ct**2 + 10*A**2*B*bz0*ct**3 + 20*A**2*B*bz0*ct - 15*A**2*B*bz2*ct**3 - 30*A**2*B*bz2*ct + 6*A**2*B*bz4*ct**3 + 12*A**2*B*bz4*ct - A**2*B*bz6*ct**3 - 2*A**2*B*bz6*ct + 20*A*B**2*bz0*ct**2 + 10*A*B**2*bz0 - 30*A*B**2*bz2*ct**2 - 15*A*B**2*bz2 + 12*A*B**2*bz4*ct**2 + 6*A*B**2*bz4 - 2*A*B**2*bz6*ct**2 - A*B**2*bz6 + 120*A*bz0*ct**2 + 60*A*bz0 - 160*A*bz2*ct**2 - 80*A*bz2 + 40*A*bz4*ct**2 + 20*A*bz4 + 10*B**3*bz0*ct - 15*B**3*bz2*ct + 6*B**3*bz4*ct - B**3*bz6*ct + 180*B*bz0*ct - 240*B*bz2*ct + 60*B*bz4*ct)/480
      d3_ABT =  -st*(10*A**3*B*bz0*ct - 15*A**3*B*bz2*ct + 6*A**3*B*bz4*ct - A**3*B*bz6*ct + 10*A**2*B**2*bz0*ct**2 + 10*A**2*B**2*bz0 - 15*A**2*B**2*bz2*ct**2 - 15*A**2*B**2*bz2 + 6*A**2*B**2*bz4*ct**2 + 6*A**2*B**2*bz4 - A**2*B**2*bz6*ct**2 - A**2*B**2*bz6 + 60*A**2*bz0 - 80*A**2*bz2 + 20*A**2*bz4 + 10*A*B**3*bz0*ct - 15*A*B**3*bz2*ct + 6*A*B**3*bz4*ct - A*B**3*bz6*ct + 180*A*B*bz0*ct - 240*A*B*bz2*ct + 60*A*B*bz4*ct + 60*B**2*bz0 - 80*B**2*bz2 + 20*B**2*bz4 + 240*bz0 - 240*bz2)/480
      d3_BBT =  -A*st*(10*A**2*B*bz0*ct**2 - 15*A**2*B*bz2*ct**2 + 6*A**2*B*bz4*ct**2 - A**2*B*bz6*ct**2 + 20*A*B**2*bz0*ct - 30*A*B**2*bz2*ct + 12*A*B**2*bz4*ct - 2*A*B**2*bz6*ct + 120*A*bz0*ct - 160*A*bz2*ct + 40*A*bz4*ct + 10*B**3*bz0 - 15*B**3*bz2 + 6*B**3*bz4 - B**3*bz6 + 180*B*bz0 - 240*B*bz2 + 60*B*bz4)/480
      d3_BTT =  A*(10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 + 10*A*B**3*bz0*st**2 - 15*A*B**3*bz2*st**2 + 6*A*B**3*bz4*st**2 - A*B**3*bz6*st**2 - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 60*B**2*bz0*ct + 80*B**2*bz2*ct - 20*B**2*bz4*ct - 240*bz0*ct + 240*bz2*ct)/480
      

      integral_t  = int_xxxx

       if (dabs(A*B) < 1.d-30) then
         if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
          der_t      =    spa * (ccd*d1_A - cqcd*d3_ABB)    /(ax3)
         else if (dabs(A) < 1.d-15) then
          der_t      =      spa * (ccd*d1_A - cqcd*d3_ABB)  /(ax3)      &
      &                 +   spa * d2_AT*sqcd                /(B2*ax3)   &
      &                 -   spa * d3_ABT*sqcd               /(B*ax3)
         else if (dabs(B) < 1.d-15) then
          der_t      =      spa * (ccd*d1_A - cqcd*d3_ABB)  /(ax3)      &
      &                 +   cpa * (ccd*d1_T - cqcd*d3_BBT)  /(A*ax3)
         endif
       else

        der_t      =      spa * (ccd*d1_A - cqcd*d3_ABB)  /(ax3)        &
      &               +   spa * d2_AT*sqcd                /(B2*ax3)     &
      &               -   spa * d3_ABT*sqcd               /(B*ax3)      &
      &               +   cpa * (ccd*d1_T - cqcd*d3_BBT)  /(A*ax3)      &
      &               +   cpa * d2_TT*sqcd                /(A*B2*ax3)   &
      &               -   cpa * d3_BTT*sqcd               /(A*B*ax3)

      end if 

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
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
      bz4   = bessi_scaled(4,z)
      bz5   = bessi_scaled(5,z)
      bz6   = bessi_scaled(6,z)

      st  = dsin(theta)        ;    ct  = dcos(theta)
      st2 = st  * st           ;    ct2 = ct  * ct  
      st3 = st2 * st           ;    ct3 = ct2 * ct  
      
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


      d1_B   =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T   = - 0.5d0 *  A * B * st  * (bz0-bz2)
      
      d2_BT  = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT  =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24
      d3_AAB =  (10*A**3*bz0*ct - 15*A**3*bz2*ct + 6*A**3*bz4*ct - A**3*bz6*ct + 20*A**2*B*bz0*ct**2 + 10*A**2*B*bz0 - 30*A**2*B*bz2*ct**2 - 15*A**2*B*bz2 + 12*A**2*B*bz4*ct**2 + 6*A**2*B*bz4 - 2*A**2*B*bz6*ct**2 - A**2*B*bz6 + 10*A*B**2*bz0*ct**3 + 20*A*B**2*bz0*ct - 15*A*B**2*bz2*ct**3 - 30*A*B**2*bz2*ct + 6*A*B**2*bz4*ct**3 + 12*A*B**2*bz4*ct - A*B**2*bz6*ct**3 - 2*A*B**2*bz6*ct + 180*A*bz0*ct - 240*A*bz2*ct + 60*A*bz4*ct + 10*B**3*bz0*ct**2 - 15*B**3*bz2*ct**2 + 6*B**3*bz4*ct**2 - B**3*bz6*ct**2 + 120*B*bz0*ct**2 + 60*B*bz0 - 160*B*bz2*ct**2 - 80*B*bz2 + 40*B*bz4*ct**2 + 20*B*bz4)/480
      d3_ABT = -st*(10*A**3*B*bz0*ct - 15*A**3*B*bz2*ct + 6*A**3*B*bz4*ct - A**3*B*bz6*ct + 10*A**2*B**2*bz0*ct**2 + 10*A**2*B**2*bz0 - 15*A**2*B**2*bz2*ct**2 - 15*A**2*B**2*bz2 + 6*A**2*B**2*bz4*ct**2 + 6*A**2*B**2*bz4 - A**2*B**2*bz6*ct**2 - A**2*B**2*bz6 + 60*A**2*bz0 - 80*A**2*bz2 + 20*A**2*bz4 + 10*A*B**3*bz0*ct - 15*A*B**3*bz2*ct + 6*A*B**3*bz4*ct - A*B**3*bz6*ct + 180*A*B*bz0*ct - 240*A*B*bz2*ct + 60*A*B*bz4*ct + 60*B**2*bz0 - 80*B**2*bz2 + 20*B**2*bz4 + 240*bz0 - 240*bz2)/480
      d3_AAT = -B*st*(10*A**3*bz0 - 15*A**3*bz2 + 6*A**3*bz4 - A**3*bz6 + 20*A**2*B*bz0*ct - 30*A**2*B*bz2*ct + 12*A**2*B*bz4*ct - 2*A**2*B*bz6*ct + 10*A*B**2*bz0*ct**2 - 15*A*B**2*bz2*ct**2 + 6*A*B**2*bz4*ct**2 - A*B**2*bz6*ct**2 + 180*A*bz0 - 240*A*bz2 + 60*A*bz4 + 120*B*bz0*ct - 160*B*bz2*ct + 40*B*bz4*ct)/480
      d3_ATT = B*(10*A**3*B*bz0*st**2 - 15*A**3*B*bz2*st**2 + 6*A**3*B*bz4*st**2 - A**3*B*bz6*st**2 + 10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 - 60*A**2*bz0*ct + 80*A**2*bz2*ct - 20*A**2*bz4*ct - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 240*bz0*ct + 240*bz2*ct)/480

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
          der_t        = sqd * (cab*d1_B - cpab*d3_AAB)        /(ax3)
        else if (dabs(A) < 1.d-15) then
          der_t        = sqd * (cab*d1_B - cpab*d3_AAB)        /(ax3)             &
      &                + cqd * (-cab*d1_T + cpab*d3_AAT)       /(B*ax3)
        else if (dabs(B) < 1.d-15) then
          der_t        = sqd * (cab*d1_B - cpab*d3_AAB)        /(ax3)             &
      &                - sqd *  d2_BT*spab                     /(A2*ax3)          &
      &                + sqd * d3_ABT*spab                     /(A*ax3)
        endif
      else
        der_t        =   sqd * (cab*d1_B - cpab*d3_AAB)        /(ax3)             &
      &                + cqd * (-cab*d1_T + cpab*d3_AAT)       /(B*ax3)           &
      &                - sqd *  d2_BT*spab                     /(A2*ax3)          &
      &                + cqd * d2_TT*spab                      /(A2*B*ax3)        &
      &                + sqd * d3_ABT*spab                     /(A*ax3)           &
                       - cqd * d3_ATT*spab                     /(A*B*ax3)
      end if

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
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
      bz4   = bessi_scaled(4,z)
      bz5   = bessi_scaled(5,z)
      bz6   = bessi_scaled(6,z)

      st  = dsin(theta)        ;    ct  = dcos(theta)
      st2 = st  * st           ;    ct2 = ct  * ct  
      st3 = st2 * st           ;    ct3 = ct2 * ct  
      
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

      
      d1_B   =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T   = - 0.5d0 *  A * B * st  * (bz0-bz2)

      d2_BT  = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT  =  A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24
      d3_AAB =  (10*A**3*bz0*ct - 15*A**3*bz2*ct + 6*A**3*bz4*ct - A**3*bz6*ct + 20*A**2*B*bz0*ct**2 + 10*A**2*B*bz0 - 30*A**2*B*bz2*ct**2 - 15*A**2*B*bz2 + 12*A**2*B*bz4*ct**2 + 6*A**2*B*bz4 - 2*A**2*B*bz6*ct**2 - A**2*B*bz6 + 10*A*B**2*bz0*ct**3 + 20*A*B**2*bz0*ct - 15*A*B**2*bz2*ct**3 - 30*A*B**2*bz2*ct + 6*A*B**2*bz4*ct**3 + 12*A*B**2*bz4*ct - A*B**2*bz6*ct**3 - 2*A*B**2*bz6*ct + 180*A*bz0*ct - 240*A*bz2*ct + 60*A*bz4*ct + 10*B**3*bz0*ct**2 - 15*B**3*bz2*ct**2 + 6*B**3*bz4*ct**2 - B**3*bz6*ct**2 + 120*B*bz0*ct**2 + 60*B*bz0 - 160*B*bz2*ct**2 - 80*B*bz2 + 40*B*bz4*ct**2 + 20*B*bz4)/480
      d3_ABT = -st*(10*A**3*B*bz0*ct - 15*A**3*B*bz2*ct + 6*A**3*B*bz4*ct - A**3*B*bz6*ct + 10*A**2*B**2*bz0*ct**2 + 10*A**2*B**2*bz0 - 15*A**2*B**2*bz2*ct**2 - 15*A**2*B**2*bz2 + 6*A**2*B**2*bz4*ct**2 + 6*A**2*B**2*bz4 - A**2*B**2*bz6*ct**2 - A**2*B**2*bz6 + 60*A**2*bz0 - 80*A**2*bz2 + 20*A**2*bz4 + 10*A*B**3*bz0*ct - 15*A*B**3*bz2*ct + 6*A*B**3*bz4*ct - A*B**3*bz6*ct + 180*A*B*bz0*ct - 240*A*B*bz2*ct + 60*A*B*bz4*ct + 60*B**2*bz0 - 80*B**2*bz2 + 20*B**2*bz4 + 240*bz0 - 240*bz2)/480
      d3_AAT = -B*st*(10*A**3*bz0 - 15*A**3*bz2 + 6*A**3*bz4 - A**3*bz6 + 20*A**2*B*bz0*ct - 30*A**2*B*bz2*ct + 12*A**2*B*bz4*ct - 2*A**2*B*bz6*ct + 10*A*B**2*bz0*ct**2 - 15*A*B**2*bz2*ct**2 + 6*A*B**2*bz4*ct**2 - A*B**2*bz6*ct**2 + 180*A*bz0 - 240*A*bz2 + 60*A*bz4 + 120*B*bz0*ct - 160*B*bz2*ct + 40*B*bz4*ct)/480
      d3_ATT = B*(10*A**3*B*bz0*st**2 - 15*A**3*B*bz2*st**2 + 6*A**3*B*bz4*st**2 - A**3*B*bz6*st**2 + 10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 - 60*A**2*bz0*ct + 80*A**2*bz2*ct - 20*A**2*bz4*ct - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 240*bz0*ct + 240*bz2*ct)/480

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15) then
          der_t        =  sqc*(cab*d1_B - cpab*d3_AAB) /(ax3)
        else if (dabs(A) < 1.d-15) then
          der_t     =   sqc*(cab*d1_B - cpab*d3_AAB)  /(ax3)            &
      &               + cqc*(-cab*d1_T + cpab*d3_AAT) /(B*ax3)
        else if (dabs(B) < 1.d-15) then
          der_t      =  sqc*(cab*d1_B - cpab*d3_AAB) /(ax3)            &
      &               - d2_BT*spab*sqc                 /(A2*ax3)       &
      &               + d3_ABT*spab*sqc                /(A*ax3)
        endif
      else
        der_t        =  sqc*(cab*d1_B - cpab*d3_AAB) /(ax3)            &
      &               + cqc*(-cab*d1_T + cpab*d3_AAT)  /(B*ax3)        &
      &               - d2_BT*spab*sqc                 /(A2*ax3)       &
      &               + cqc*d2_TT*spab                 /(A2*B*ax3)     &
      &               + d3_ABT*spab*sqc                /(A*ax3)        &
      &               - cqc*d3_ATT*spab                /(A*B*ax3)
      end if 

      f           = der_t * integral_t

      end function f1110

      double precision function f1111(theta) Result(f)

      double precision,intent(in) :: theta
      double precision            :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12

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
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)
      int_xxxx    = const * 0.5d0 * pi * dsqrt(inv) * erfcx(term) / (p+q)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z
      z5 = z4 * z 
      z6 = z5 * z
      z7 = z6 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)
      bz5   = bessi_scaled(5,z)
      bz6   = bessi_scaled(6,z)
      bz7   = bessi_scaled(7,z)
      bz8   = bessi_scaled(8,z)

      st  = dsin(theta)        ;    ct  = dcos(theta)
      st2 = st  * st           ;    ct2 = ct  * ct
      st3 = st2 * st           ;    ct3 = ct2 * ct
      st4 = st3 * st           ;    ct4 = ct3 * ct
      

      cab =  dcos(xpa) * dcos(xpb)
      ccd  = dcos(xqc) * dcos(xqd)

      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))

      d1_T = - 0.5d0 * A * B * st  * (bz0-bz2)

      d2_AA = (3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0*ct**2 - 4*B**2*bz2*ct**2 + B**2*bz4*ct**2 + 12*bz0 - 12*bz2)/24
      d2_BB = (3*A**2*bz0*ct**2 - 4*A**2*bz2*ct**2 + A**2*bz4*ct**2 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT = A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24
      d2_AT = -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24


      d3_AAT = -B*st*(10*A**3*bz0 - 15*A**3*bz2 + 6*A**3*bz4 - A**3*bz6 + 20*A**2*B*bz0*ct - 30*A**2*B*bz2*ct + 12*A**2*B*bz4*ct - 2*A**2*B*bz6*ct + 10*A*B**2*bz0*ct**2 - 15*A*B**2*bz2*ct**2 + 6*A*B**2*bz4*ct**2 - A*B**2*bz6*ct**2 + 180*A*bz0 - 240*A*bz2 + 60*A*bz4 + 120*B*bz0*ct - 160*B*bz2*ct + 40*B*bz4*ct)/480
      d3_ATT =  B*(10*A**3*B*bz0*st**2 - 15*A**3*B*bz2*st**2 + 6*A**3*B*bz4*st**2 - A**3*B*bz6*st**2 + 10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 - 60*A**2*bz0*ct + 80*A**2*bz2*ct - 20*A**2*bz4*ct - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 240*bz0*ct + 240*bz2*ct)/480
      d3_BBT = -A*st*(10*A**2*B*bz0*ct**2 - 15*A**2*B*bz2*ct**2 + 6*A**2*B*bz4*ct**2 - A**2*B*bz6*ct**2 + 20*A*B**2*bz0*ct - 30*A*B**2*bz2*ct + 12*A*B**2*bz4*ct - 2*A*B**2*bz6*ct + 120*A*bz0*ct - 160*A*bz2*ct + 40*A*bz4*ct + 10*B**3*bz0 - 15*B**3*bz2 + 6*B**3*bz4 - B**3*bz6 + 180*B*bz0 - 240*B*bz2 + 60*B*bz4)/480
      d3_BTT = A*(10*A**2*B**2*bz0*ct*st**2 - 15*A**2*B**2*bz2*ct*st**2 + 6*A**2*B**2*bz4*ct*st**2 - A**2*B**2*bz6*ct*st**2 + 10*A*B**3*bz0*st**2 - 15*A*B**3*bz2*st**2 + 6*A*B**3*bz4*st**2 - A*B**3*bz6*st**2 - 60*A*B*bz0*ct**2 + 120*A*B*bz0*st**2 + 80*A*B*bz2*ct**2 - 160*A*B*bz2*st**2 - 20*A*B*bz4*ct**2 + 40*A*B*bz4*st**2 - 60*B**2*bz0*ct + 80*B**2*bz2*ct - 20*B**2*bz4*ct - 240*bz0*ct + 240*bz2*ct)/480


      d4_AABT = -st*(35*A**4*B*bz0*ct - 56*A**4*B*bz2*ct + 28*A**4*B*bz4*ct - 8*A**4*B*bz6*ct + A**4*B*bz8*ct + 70*A**3*B**2*bz0*ct**2 + 35*A**3*B**2*bz0 - 112*A**3*B**2*bz2*ct**2 - 56*A**3*B**2*bz2 + 56*A**3*B**2*bz4*ct**2 + 28*A**3*B**2*bz4 - 16*A**3*B**2*bz6*ct**2 - 8*A**3*B**2*bz6 + 2*A**3*B**2*bz8*ct**2 + A**3*B**2*bz8 + 280*A**3*bz0 - 420*A**3*bz2 + 168*A**3*bz4 - 28*A**3*bz6 + 35*A**2*B**3*bz0*ct**3 + 70*A**2*B**3*bz0*ct - 56*A**2*B**3*bz2*ct**3 - 112*A**2*B**3*bz2*ct + 28*A**2*B**3*bz4*ct**3 + 56*A**2*B**3*bz4*ct - 8*A**2*B**3*bz6*ct**3 - 16*A**2*B**3*bz6*ct + A**2*B**3*bz8*ct**3 + 2*A**2*B**3*bz8*ct + 1960*A**2*B*bz0*ct - 2940*A**2*B*bz2*ct + 1176*A**2*B*bz4*ct - 196*A**2*B*bz6*ct + 35*A*B**4*bz0*ct**2 - 56*A*B**4*bz2*ct**2 + 28*A*B**4*bz4*ct**2 - 8*A*B**4*bz6*ct**2 + A*B**4*bz8*ct**2 + 1400*A*B**2*bz0*ct**2 + 840*A*B**2*bz0 - 2100*A*B**2*bz2*ct**2 - 1260*A*B**2*bz2 + 840*A*B**2*bz4*ct**2 + 504*A*B**2*bz4 - 140*A*B**2*bz6*ct**2 - 84*A*B**2*bz6 + 5040*A*bz0 - 6720*A*bz2 + 1680*A*bz4 + 560*B**3*bz0*ct - 840*B**3*bz2*ct + 336*B**3*bz4*ct - 56*B**3*bz6*ct + 6720*B*bz0*ct - 8960*B*bz2*ct + 2240*B*bz4*ct)/13440
      d4_ABBT = -st*(35*A**4*B*bz0*ct**2 - 56*A**4*B*bz2*ct**2 + 28*A**4*B*bz4*ct**2 - 8*A**4*B*bz6*ct**2 + A**4*B*bz8*ct**2 + 35*A**3*B**2*bz0*ct**3 + 70*A**3*B**2*bz0*ct - 56*A**3*B**2*bz2*ct**3 - 112*A**3*B**2*bz2*ct + 28*A**3*B**2*bz4*ct**3 + 56*A**3*B**2*bz4*ct - 8*A**3*B**2*bz6*ct**3 - 16*A**3*B**2*bz6*ct + A**3*B**2*bz8*ct**3 + 2*A**3*B**2*bz8*ct + 560*A**3*bz0*ct - 840*A**3*bz2*ct + 336*A**3*bz4*ct - 56*A**3*bz6*ct + 70*A**2*B**3*bz0*ct**2 + 35*A**2*B**3*bz0 - 112*A**2*B**3*bz2*ct**2 - 56*A**2*B**3*bz2 + 56*A**2*B**3*bz4*ct**2 + 28*A**2*B**3*bz4 - 16*A**2*B**3*bz6*ct**2 - 8*A**2*B**3*bz6 + 2*A**2*B**3*bz8*ct**2 + A**2*B**3*bz8 + 1400*A**2*B*bz0*ct**2 + 840*A**2*B*bz0 - 2100*A**2*B*bz2*ct**2 - 1260*A**2*B*bz2 + 840*A**2*B*bz4*ct**2 + 504*A**2*B*bz4 - 140*A**2*B*bz6*ct**2 - 84*A**2*B*bz6 + 35*A*B**4*bz0*ct - 56*A*B**4*bz2*ct + 28*A*B**4*bz4*ct - 8*A*B**4*bz6*ct + A*B**4*bz8*ct + 1960*A*B**2*bz0*ct - 2940*A*B**2*bz2*ct + 1176*A*B**2*bz4*ct - 196*A*B**2*bz6*ct + 6720*A*bz0*ct - 8960*A*bz2*ct + 2240*A*bz4*ct + 280*B**3*bz0 - 420*B**3*bz2 + 168*B**3*bz4 - 28*B**3*bz6 + 5040*B*bz0 - 6720*B*bz2 + 1680*B*bz4)/13440
      d4_AABB =     (35*A**4*bz0*ct**2 - 56*A**4*bz2*ct**2 + 28*A**4*bz4*ct**2 - 8*A**4*bz6*ct**2 + A**4*bz8*ct**2 + 70*A**3*B*bz0*ct**3 + 70*A**3*B*bz0*ct - 112*A**3*B*bz2*ct**3 - 112*A**3*B*bz2*ct + 56*A**3*B*bz4*ct**3 + 56*A**3*B*bz4*ct - 16*A**3*B*bz6*ct**3 - 16*A**3*B*bz6*ct + 2*A**3*B*bz8*ct**3 + 2*A**3*B*bz8*ct + 35*A**2*B**2*bz0*ct**4 + 140*A**2*B**2*bz0*ct**2 + 35*A**2*B**2*bz0 - 56*A**2*B**2*bz2*ct**4 - 224*A**2*B**2*bz2*ct**2 - 56*A**2*B**2*bz2 + 28*A**2*B**2*bz4*ct**4 + 112*A**2*B**2*bz4*ct**2 + 28*A**2*B**2*bz4 - 8*A**2*B**2*bz6*ct**4 - 32*A**2*B**2*bz6*ct**2 - 8*A**2*B**2*bz6 + A**2*B**2*bz8*ct**4 + 4*A**2*B**2*bz8*ct**2 + A**2*B**2*bz8 + 1400*A**2*bz0*ct**2 + 280*A**2*bz0 - 2100*A**2*bz2*ct**2 - 420*A**2*bz2 + 840*A**2*bz4*ct**2 + 168*A**2*bz4 - 140*A**2*bz6*ct**2 - 28*A**2*bz6 + 70*A*B**3*bz0*ct**3 + 70*A*B**3*bz0*ct - 112*A*B**3*bz2*ct**3 - 112*A*B**3*bz2*ct + 56*A*B**3*bz4*ct**3 + 56*A*B**3*bz4*ct - 16*A*B**3*bz6*ct**3 - 16*A*B**3*bz6*ct + 2*A*B**3*bz8*ct**3 + 2*A*B**3*bz8*ct + 1120*A*B*bz0*ct**3 + 2240*A*B*bz0*ct - 1680*A*B*bz2*ct**3 - 3360*A*B*bz2*ct + 672*A*B*bz4*ct**3 + 1344*A*B*bz4*ct - 112*A*B*bz6*ct**3 - 224*A*B*bz6*ct + 35*B**4*bz0*ct**2 - 56*B**4*bz2*ct**2 + 28*B**4*bz4*ct**2 - 8*B**4*bz6*ct**2 + B**4*bz8*ct**2 + 1400*B**2*bz0*ct**2 + 280*B**2*bz0 - 2100*B**2*bz2*ct**2 - 420*B**2*bz2 + 840*B**2*bz4*ct**2 + 168*B**2*bz4 - 140*B**2*bz6*ct**2 - 28*B**2*bz6 + 3360*bz0*ct**2 + 1680*bz0 - 4480*bz2*ct**2 - 2240*bz2 + 1120*bz4*ct**2 + 560*bz4)/13440
      d4_ABTT =     (35*A**4*B**2*bz0*ct*st**2 - 56*A**4*B**2*bz2*ct*st**2 + 28*A**4*B**2*bz4*ct*st**2 - 8*A**4*B**2*bz6*ct*st**2 + A**4*B**2*bz8*ct*st**2 + 35*A**3*B**3*bz0*ct**2*st**2 + 35*A**3*B**3*bz0*st**2 - 56*A**3*B**3*bz2*ct**2*st**2 - 56*A**3*B**3*bz2*st**2 + 28*A**3*B**3*bz4*ct**2*st**2 + 28*A**3*B**3*bz4*st**2 - 8*A**3*B**3*bz6*ct**2*st**2 - 8*A**3*B**3*bz6*st**2 + A**3*B**3*bz8*ct**2*st**2 + A**3*B**3*bz8*st**2 - 280*A**3*B*bz0*ct**2 + 560*A**3*B*bz0*st**2 + 420*A**3*B*bz2*ct**2 - 840*A**3*B*bz2*st**2 - 168*A**3*B*bz4*ct**2 + 336*A**3*B*bz4*st**2 + 28*A**3*B*bz6*ct**2 - 56*A**3*B*bz6*st**2 + 35*A**2*B**4*bz0*ct*st**2 - 56*A**2*B**4*bz2*ct*st**2 + 28*A**2*B**4*bz4*ct*st**2 - 8*A**2*B**4*bz6*ct*st**2 + A**2*B**4*bz8*ct*st**2 - 280*A**2*B**2*bz0*ct**3 + 1400*A**2*B**2*bz0*ct*st**2 - 280*A**2*B**2*bz0*ct + 420*A**2*B**2*bz2*ct**3 - 2100*A**2*B**2*bz2*ct*st**2 + 420*A**2*B**2*bz2*ct - 168*A**2*B**2*bz4*ct**3 + 840*A**2*B**2*bz4*ct*st**2 - 168*A**2*B**2*bz4*ct + 28*A**2*B**2*bz6*ct**3 - 140*A**2*B**2*bz6*ct*st**2 + 28*A**2*B**2*bz6*ct - 1680*A**2*bz0*ct + 2240*A**2*bz2*ct - 560*A**2*bz4*ct - 280*A*B**3*bz0*ct**2 + 560*A*B**3*bz0*st**2 + 420*A*B**3*bz2*ct**2 - 840*A*B**3*bz2*st**2 - 168*A*B**3*bz4*ct**2 + 336*A*B**3*bz4*st**2 + 28*A*B**3*bz6*ct**2 - 56*A*B**3*bz6*st**2 - 5040*A*B*bz0*ct**2 + 6720*A*B*bz0*st**2 + 6720*A*B*bz2*ct**2 - 8960*A*B*bz2*st**2 - 1680*A*B*bz4*ct**2 + 2240*A*B*bz4*st**2 - 1680*B**2*bz0*ct + 2240*B**2*bz2*ct - 560*B**2*bz4*ct - 6720*bz0*ct + 6720*bz2*ct)/13440

      integral_t  = int_xxxx

      if (dabs(A*B) < 1.d-30) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
          t1     =   cab*ccd*bz0                        /(ax4)
          t2     = - cab*cqcd*d2_BB                     /(ax4)
          t3     = - ccd*cpab*d2_AA                     /(ax4)
          t4     =   cpab*cqcd*d4_AABB                  /(ax4)
          der_t = t1 + t2 + t3 + t4
        else if (dabs(A) < 1.d-15) then 
          t1     =   cab*ccd*bz0                        /(ax4)
          t2     = - cab*cqcd*d2_BB                     /(ax4)
          t3     = - ccd*cpab*d2_AA                     /(ax4)
          t4     =   cpab*cqcd*d4_AABB                  /(ax4)
          t5     =   sqcd*(cab*d1_T - cpab*d3_AAT)     /(B2*ax4)
          t6     =   sqcd*(-cab*d2_BT + cpab*d4_AABT)   /(B*ax4)
          der_t = t1 + t2 + t3 + t4 + t5 + t6
        else if (dabs(B) < 1.d-15) then 
          t1     =   cab*ccd*bz0                        /(ax4)
          t2     = - cab*cqcd*d2_BB                     /(ax4)
          t3     = - ccd*cpab*d2_AA                     /(ax4)
          t4     =   cpab*cqcd*d4_AABB                  /(ax4)
          t7     =   spab*(-ccd*d1_T + cqcd*d3_BBT)     /(A2*ax4)
          t10    =   spab*(ccd*d2_AT - cqcd*d4_ABBT)    /(A*ax4)
          der_t = t1 + t2 + t3 + t4 + t7 + t10
        end if 
      
      else 

        t1     =   cab*ccd*bz0                        /(ax4)
        t2     = - cab*cqcd*d2_BB                     /(ax4)
        t3     = - ccd*cpab*d2_AA                     /(ax4)
        t4     =   cpab*cqcd*d4_AABB                  /(ax4)
        t5     =   sqcd*(cab*d1_T - cpab*d3_AAT)     /(B2*ax4)
        t6     =   sqcd*(-cab*d2_BT + cpab*d4_AABT)   /(B*ax4)
        t7     =   spab*(-ccd*d1_T + cqcd*d3_BBT)     /(A2*ax4)
        t8     = - spab*sqcd*d2_TT                    /(A2*B2*ax4)
        t9     =   spab*sqcd*d3_BTT                   /(A2*B*ax4)
        t10    =   spab*(ccd*d2_AT - cqcd*d4_ABBT)    /(A*ax4)
        t11    =   spab*sqcd*d3_ATT                   /(A*B2*ax4)
        t12    = - spab*sqcd*d4_ABTT                  /(A*B*ax4)
        der_t = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12

      end if 

      f           = der_t * integral_t

      end function f1111

      double precision function f0022(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * ppq * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))
      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / ( q * ppq3 )
      
      ! - derivative part - !

      integral_t  = int_xxyy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0022

      double precision function f0033t1(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * ppq * dsqrt(beta)
      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term1 ) / ( q * ppq3 )
      
      ! - derivative part - !

      integral_t  = int_xxyy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0033t1

      double precision function f0033t2(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))
      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term2 ) / ( q * ppq3 )
      
      ! - derivative part - !

      integral_t  = int_xxyy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0033t2


      double precision function f0122(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * ppq * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))
      int_xxyy    = const * 0.125d0 * dsqrt(pi) * ( term1 + term2  ) / ( q * ppq3 )

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      spb  = dsin(xpb)
      cpb  = dcos(xpb)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_xxyy

      if (A < 1.d-15) then 
        der_t       = (d1_A*spb)/(ax)
      else 
        der_t       = (d1_A*spb + cpb*d1_T/A)/(ax)
      end if 

      f           = der_t * integral_t

      end function f0122

      
      double precision function f0202(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_xyxy = const * 0.125d0 * ( term1 + term2 ) / (ppq3)
      
      ! - derivative part - !

      integral_t  = int_xyxy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0202

      double precision function f0212(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_xyxy = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqc  = dsin(xqc)
      cqc  = dcos(xqc)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_xyxy

      if (B < 1.d-15) then 
        der_t       = (d1_B*sqc)/ax
      else 
        der_t       = (d1_B*sqc - cqc*d1_T/B)/ax
      end if 

      f           = der_t * integral_t

      end function f0212

      
      double precision function f0220(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 
      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_xyxy = const * 0.125d0 * ( term1 + term2 ) / (ppq3)
      
      ! - derivative part - !

      integral_t  = int_xyxy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f0220

      double precision function f0221(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_xyxy = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqd  = dsin(xqd)
      cqd  = dcos(xqd)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_xyxy

      if (B < 1.d-15) then 
        der_t       = (d1_B*sqd)/ax
      else
        der_t       = (d1_B*sqd - cqd*d1_T/B)/ax
      end if 

      f           = der_t * integral_t

      end function f0221

      double precision function f1022(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = 2.d0 * p * ppq * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))
      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / ( q * ppq3 ) 

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)

      st = dsin(theta)       ;    ct = dcos(theta)

      spa  = dsin(xpa)
      cpa  = dcos(xpa)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_xxyy
      if (A < 1.d-15) then 
        der_t       = d1_A*spa/ax
      else 
        der_t       = (d1_A*spa + cpa*d1_T/A)/ax
      end if 
      f           = der_t * integral_t

      end function f1022

      double precision function f1122(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * ppq * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi) * dsqrt(inv) * (3.d0*p*q+2.d0*q2+p2*(1.d0-2.d0*q*beta))
      int_xxyy = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / ( q * ppq3 )

      ! - derivative part - !

      z2 = z  * z
      z3 = z2 * z
      z4 = z3 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      cab  = dcos(xpa) *  dcos(xpb)
      spab = dsin(ax*(2.d0*xp-xa-xb))
      cpab = dcos(ax*(2.d0*xp-xa-xb))
      
      d1_T  = - 0.5d0 * A * B * st    * (bz0-bz2)
      d2_AA =  (3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0*ct**2 - 4*B**2*bz2*ct**2 + B**2*bz4*ct**2 + 12*bz0 - 12*bz2)/24
      d2_AT =  -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24

      integral_t  = int_xxyy

      if (A < 1.d-15) then 
        der_t       = (cab*bz0 - cpab*d2_AA) / ax2
      else 
        der_t       = (cab*bz0 - cpab*d2_AA - d1_T*spab/A2 + d2_AT*spab/A) / ax2
      end if 
      
      f           = der_t * integral_t

      end function f1122


      double precision function f1202(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      spa  = dsin(xpa)
      cpa  = dcos(xpa)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yxyx

      if (A < 1.d-15) then 
        der_t       = (d1_A*spa)/ax
      else 
        der_t       = (d1_A*spa + cpa*d1_T/A)/(ax)
      end if 

      f           = der_t * integral_t

      end function f1202

      double precision function f1212(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqc = dsin(xqc)
      cqc = dcos(xqc)
      
      
      d2_AB =            (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT =  -B * st * (3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT =  -A * st * (3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =   A * B  * (3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_yxyx
      if (dabs(A*B) < 1.d-30 ) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
          der_t     = (d2_AB*spa*sqc)/(ax2)
        else if (dabs(A) < 1.d-15) then 
          der_t     = (d2_AB*spa*sqc - cqc*d2_AT*spa/B)/(ax2)
        else if (dabs(B) < 1.d-15) then 
          der_t     = (d2_AB*spa*sqc + cpa*d2_BT*sqc/A)/(ax2)
        end if 
      else 
        der_t     = (d2_AB*spa*sqc - cqc*d2_AT*spa/B + cpa*d2_BT*sqc/A - cpa*cqc*d2_TT/(A*B))/(ax2)
      end if 
      f           = der_t * integral_t

      end function f1212

      double precision function f1220(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      spa  = dsin(xpa)
      cpa  = dcos(xpa)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yxyx
      if (A < 1.d-15) then 
        der_t       = d1_A*spa/ax
      else 
        der_t       = (d1_A*spa + cpa*d1_T/A)/ax
      end if 

      f             = der_t * integral_t

      end function f1220

      double precision function f1221(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spa = dsin(xpa)
      cpa = dcos(xpa)
      sqd = dsin(xqd)
      cqd = dcos(xqd)
      
      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT =  -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT =  -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =   A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_yxyx

      if (dabs(A*B) < 1.d-30 ) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
          der_t       = (d2_AB*spa*sqd)/(ax2)
        else if (dabs(A) < 1.d-15) then 
          der_t       = (d2_AB*spa*sqd - cqd*d2_AT*spa/B)/(ax2)
        else if (dabs(B) < 1.d-15) then 
          der_t       = (d2_AB*spa*sqd + cpa*d2_BT*sqd/A)/(ax2)
        end if 
      else 
        der_t       = (d2_AB*spa*sqd - cqd*d2_AT*spa/B + cpa*d2_BT*sqd/A - cpa*cqd*d2_TT/(A*B))/(ax2)
      end if 
        
      f           = der_t * integral_t

      end function f1221

      double precision function f2002(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)
      
      ! - derivative part - !

      integral_t  = int_yxyx
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2002


      double precision function f2012(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqc  = dsin(xqc)
      cqc  = dcos(xqc)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yxyx
      if (B < 1.d-15) then 
        der_t       = d1_B*sqc/ax
      else 
        der_t       = (d1_B*sqc - cqc*d1_T/B)/ax
      end if 
      
      f           = der_t * integral_t

      end function f2012

      double precision function f2020(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)
      
      ! - derivative part - !

      integral_t  = int_yxyx
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2020

      double precision function f2021(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqd  = dsin(xqd)
      cqd  = dcos(xqd)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)
      

      integral_t  = int_yxyx
      if (B < 1.d-15) then 
        der_t       = d1_B*sqd/ax
      else 
        der_t       = (d1_B*sqd - cqd*d1_T/B)/ax
      end if 
      
      f           = der_t * integral_t

      end function f2021

      double precision function f2102(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      spb  = dsin(xpb)
      cpb  = dcos(xpb)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yxyx

      if (A < 1.d-15) then 
        der_t       = d1_A*spb/ax
      else 
        der_t       = (d1_A*spb + cpb*d1_T/A)/ax
      end if 

      f           = der_t * integral_t

      end function f2102



      double precision function f2112(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      

      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqc = dsin(xqc)
      cqc = dcos(xqc)

      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))
      
      
      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT =  -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT =  -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =   A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_yxyx
      
      if (dabs(A*B) < 1.d-30 ) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
          der_t       = (d2_AB*spb*sqc)/ax2
        else if (dabs(A) < 1.d-15) then 
          der_t       = (d2_AB*spb*sqc - cqc*d2_AT*spb/B)/ax2
        else if (dabs(B) < 1.d-15) then 
          der_t       = (d2_AB*spb*sqc - cpb*d2_BT*sqc/A)/ax2
        end if 
      else 
        der_t       = (d2_AB*spb*sqc - cqc*d2_AT*spb/B + cpb*d2_BT*sqc/A - cpb*cqc*d2_TT/(A*B))/ax2
      end if

      f           = der_t * integral_t

      end function f2112

      double precision function f2120(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      spb  = dsin(xpb)
      cpb  = dcos(xpb)
      
      d1_A  =   0.5d0 * (A + B * ct) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yxyx

      if (A < 1.d-15) then 
        der_t       = (d1_A*spb)/ax
      else 
        der_t       = (d1_A*spb + cpb*d1_T/A)/ax
      end if 

      
      f           = der_t * integral_t

      end function f2120

      double precision function f2121(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = -2.d0 * dsqrt(pi*beta) * ppq
      term2    = erfcx(term) * pi * dsqrt(inv) * ( ppq + 2.d0 * p * q * beta )
      int_yxyx = const * 0.125d0 * ( term1 + term2 ) / (ppq3)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z  

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      spb = dsin(xpb)
      cpb = dcos(xpb)
      sqd = dsin(xqd)
      cqd = dcos(xqd)
      
      d2_AB = (3*A**2*bz0*ct - 4*A**2*bz2*ct + A**2*bz4*ct + 3*A*B*bz0*ct**2 + 3*A*B*bz0 - 4*A*B*bz2*ct**2 - 4*A*B*bz2 + A*B*bz4*ct**2 + A*B*bz4 + 3*B**2*bz0*ct - 4*B**2*bz2*ct + B**2*bz4*ct + 12*bz0*ct - 12*bz2*ct)/24
      d2_AT =  -B*st*(3*A**2*bz0 - 4*A**2*bz2 + A**2*bz4 + 3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 12*bz0 - 12*bz2)/24
      d2_BT =  -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_TT =   A*B*(3*A*B*bz0*st**2 - 4*A*B*bz2*st**2 + A*B*bz4*st**2 - 12*bz0*ct + 12*bz2*ct)/24

      integral_t  = int_yxyx


      if (dabs(A*B) < 1.d-30 ) then 
        if (dabs(A) < 1.d-15 .and. dabs(B) < 1.d-15 ) then 
          der_t       = (d2_AB*spb*sqd)/ax2
        else if (dabs(A) < 1.d-15) then 
          der_t       = (d2_AB*spb*sqd - cqd*d2_AT*spb/B)/ax2
        else if (dabs(B) < 1.d-15) then 
          der_t       = (d2_AB*spb*sqd + cpb*d2_BT*sqd/A)/ax2
        end if 
      else 
        der_t       = (d2_AB*spb*sqd - cqd*d2_AT*spb/B + cpb*d2_BT*sqd/A - cpb*cqd*d2_TT/(A*B))/ax2
      end if

      f           = der_t * integral_t

      end function f2121

      double precision function f2200(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q
      ppq2        = ppq * ppq
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * q2 * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi)  * (2.d0*p2 + q2 + p * q * (3.d0-2.d0*q*beta)) / dsqrt(inv)
      int_yyxx = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / (p2*q*ppq2)
      
      ! - derivative part - !

      integral_t  = int_yyxx
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2200

      double precision function f2201(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq2        = ppq * ppq 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = 2.d0 * p * q2 * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi)  * (2.d0*p2 + q2 + p * q * (3.d0-2.d0*q*beta)) / dsqrt(inv)
      int_yyxx = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / (p2*q*ppq2)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqd  = dsin(xqd)
      cqd  = dcos(xqd)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yyxx
      if (B < 1.d-15) then 
        der_t       = (d1_B*sqd)/ax
      else 
        der_t       = (d1_B*sqd - cqd*d1_T/B)/ax
      end if 

      f           = der_t * integral_t

      end function f2201

      double precision function f2210(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq2        = ppq * ppq 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1    = 2.d0 * p * q2 * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi)  * (2.d0*p2 + q2 + p * q * (3.d0-2.d0*q*beta)) / dsqrt(inv)
      int_yyxx = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / (p2*q*ppq2)

      ! - derivative part - !

      z2 = z * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)


      st = dsin(theta)       ;    ct = dcos(theta)

      sqc  = dsin(xqc)
      cqc  = dcos(xqc)
      
      d1_B  =   0.5d0 * (A * ct + B) * (bz0-bz2)
      d1_T  = - 0.5d0 *  A * B * st  * (bz0-bz2)

      integral_t  = int_yyxx
      if (B < 1.d-15) then 
        der_t       = (d1_B*sqc)/ax
      else 
        der_t       = (d1_B*sqc - cqc*d1_T/B)/ax
      end if 
      
      f           = der_t * integral_t

      end function f2210

      double precision function f2211(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      ppq         = p+q 
      ppq2        = ppq * ppq 
      ppq3        = ppq * ppq * ppq

      p2          = p * p 
      q2          = q * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = 2.d0 * p * q2 * dsqrt(beta)
      term2    = erfcx(term) * dsqrt(pi)  * (2.d0*p2 + q2 + p * q * (3.d0-2.d0*q*beta)) / dsqrt(inv)
      int_yyxx = const * 0.125d0 * dsqrt(pi) * ( term1 + term2 ) / (p2*q*ppq2)

      ! - derivative part - !

      z2 = z  * z 
      z3 = z2 * z
      z4 = z3 * z 

      bz0   = bessi_scaled(0,z)
      bz1   = bessi_scaled(1,z)
      bz2   = bessi_scaled(2,z)
      bz3   = bessi_scaled(3,z)
      bz4   = bessi_scaled(4,z)


      st  = dsin(theta)       ;    ct = dcos(theta)
      st2 = st * st 
      
      ccd  = dcos(xqc) *  dcos(xqd)
      sqcd = dsin(ax*(2.d0*xq-xc-xd))
      cqcd = dcos(ax*(2.d0*xq-xc-xd))
      
      d1_T  = - 0.5d0 * A * B * st  * (bz0-bz2)
      d2_BB =   (3*A**2*bz0*ct**2 - 4*A**2*bz2*ct**2 + A**2*bz4*ct**2 + 6*A*B*bz0*ct - 8*A*B*bz2*ct + 2*A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24
      d2_BT = -A*st*(3*A*B*bz0*ct - 4*A*B*bz2*ct + A*B*bz4*ct + 3*B**2*bz0 - 4*B**2*bz2 + B**2*bz4 + 12*bz0 - 12*bz2)/24

      integral_t  = int_yyxx

      if (dabs(B) < 1.d-15 ) then 
        der_t   = (ccd*bz0 - cqcd*d2_BB)/ax2 
      else 
        der_t   = (ccd*bz0 - cqcd*d2_BB + d1_T*sqcd/B2 - d2_BT*sqcd/B )/ax2 
      end if

      f           = der_t * integral_t

      end function f2211

      double precision function f2222(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      beta2       = beta * beta
      ppq         = p+q 
      ppq2        = ppq  * ppq 
      ppq3        = ppq2 * ppq
      ppq4        = ppq3 * ppq

      p2          = p  * p 
      p3          = p2 * p 
      q2          = q  * q
      q3          = q2 * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    = erfcx(term) * pi * ( ppq2 * (4.d0 * p + q ) * ( p + 4.d0 * q ) - 4.d0 * p * ( p - 2.d0 * q) * ( 2.d0 * p - q ) * q * ppq * beta  + 12.d0 * p3 * q3 * beta2 )
      term2    = 2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * dsqrt(beta) * (4.d0 * p3 - 3.d0 * p * q2 + 4.d0 * q3 - 3.d0 * p2 * q * (1.d0+2.d0*q*beta) ) 
      int_yyyy = const * 0.015625d0 * ( term1 + term2 ) / ( p2 * dsqrt(inv) * q2 * ppq4 )
      
      ! - derivative part - !

      integral_t  = int_yyyy
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2222

      double precision function f2233(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      beta2       = beta * beta 
      ppq         = p+q 
      ppq2        = ppq  * ppq 
      ppq3        = ppq2 * ppq
      ppq4        = ppq3 * ppq

      p2          = p  * p 
      p3          = p2 * p 
      q2          = q  * q
      q3          = q2 * q

      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)


      term1 = erfcx(term) * pi * ( ppq2 * ( 4.d0 * p2 + 11.d0 * p * q + 4.d0 * q2 ) - 4.d0 * p * q * ppq * ( 2.d0 * p2 + p * q + 2.d0 * q2 ) * beta + 4.d0 * p3 * q3 * beta2 )
      term2 = 2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * dsqrt(beta) * ( 4.d0 * p3 + 7.d0 * p * q2 + 4.d0 * q3 + p2 * q * (7.d0 - 2.d0 * q * beta) )
      int_yyzz = const * 0.015625d0 * (term1 + term2 ) / (p2 * dsqrt(inv) * q2 * ppq4)
      
      ! - derivative part - !

      integral_t  = int_yyzz
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2233
      
      double precision function f2323(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      beta2       = beta * beta 
      ppq         = p+q 
      ppq2        = ppq  * ppq 
      ppq3        = ppq2 * ppq
      ppq4        = ppq3 * ppq

      p2          = p  * p 
      p3          = p2 * p 
      q2          = q  * q
      q3          = q2 * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    =   2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * sqrt(beta) * ( 5.d0 * ppq + 2.d0 * p * q * beta )
      term2    =  - erfcx(term) * pi * (3.d0 * ppq2 + 12.d0 * p * q * ppq * beta + 4.d0 * p2 * q2 * beta2)
      int_yzyz =  - const * 0.015625d0 * (term1 + term2) / (p * dsqrt(inv) * q * ppq4)
      
      ! - derivative part - !

      integral_t  = int_yzyz
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2323

      double precision function f2332(theta) Result(f)

      double precision,intent(in) :: theta

      ! - integral part - !

      A           = 2.d0 * p_x / (ax2)
      B           = 2.d0 * q_x / (ax2)
      A2          = A*A
      B2          = B*B
      inv         = (p+q)/(p*q)
      Iinv        = 1.d0/inv
      psi         = phi - theta 

      beta        = 2.d0 * (1.d0-dcos(psi)) / ax2
      beta2       = beta * beta 
      ppq         = p+q 
      ppq2        = ppq  * ppq 
      ppq3        = ppq2 * ppq
      ppq4        = ppq3 * ppq

      p2          = p  * p 
      p3          = p2 * p 
      q2          = q  * q
      q3          = q2 * q


      z           = dsqrt(  dabs( A2 + B2 + 2.d0 * A * B * dcos(theta) ) )
      term        = dsqrt(  dabs( 2.d0 * Iinv * ( 1.d0 - dcos(psi) ) )   ) / ax
      const       = 0.5d0 * pi * dexp(z-2.d0*(p+q)/ax2+sf)

      term1    =   2.d0 * p * dsqrt(pi) * dsqrt(inv) * q * sqrt(beta) * ( 5.d0 * ppq + 2.d0 * p * q * beta )
      term2    =  - erfcx(term) * pi * (3.d0 * ppq2 + 12.d0 * p * q * ppq * beta + 4.d0 * p2 * q2 * beta2)
      int_yzyz =  - const * 0.015625d0 * (term1 + term2) / (p * dsqrt(inv) * q * ppq4)
      
      ! - derivative part - !

      integral_t  = int_yzyz
      der_t       = bessi_scaled(0, z)
      f           = der_t * integral_t

      end function f2332

end subroutine integrate_ERI_integral_mod