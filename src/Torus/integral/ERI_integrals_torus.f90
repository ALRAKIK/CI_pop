subroutine ERI_integral_4_function_torus(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 
      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l 
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision,parameter      :: pi     = 3.14159265358979323846D00
      double precision,parameter      :: R2PI52 = 2.0d0*pi**(5.0d0/2.0d0)
      
      double precision                :: xa , ya , za 
      double precision                :: xb , yb , zb
      double precision                :: xc , yc , zc
      double precision                :: xd , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: zAB , zCD
      double precision                ::  DAB ,  DCD
      double precision                :: D2AB , D2CD 
      double precision                :: xp , xq 
      double precision                :: yp , yq
      double precision                :: zp , zq
      double precision                :: xpa , ypa , zpa
      double precision                :: xpb , ypb , zpb 
      double precision                :: xqc , yqc , zqc 
      double precision                :: xqd , yqd , zqd 
      double precision                :: c1 , c2 , c3 , c4 
      double precision                :: mu1 , mu2 
      double precision                :: E_AB , E_CD 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: p , q  , ip , iq 
      double precision                :: p_plus_q  , ip_plus_q
      double precision                :: p_times_q , ip_times_q
      double precision                :: popq , qopq , pq 
      double precision                :: tu   , tv 
      double precision                :: c_boys
      double precision                :: const  
      double precision                :: F_in 
      double precision                :: Boys_func
      double precision                :: F0 , F1 , F2 , F3 , F4 
      double precision                :: xPQ , yPQ , zPQ
      double precision                :: R2PQ
      double precision                :: value_s 

      !-----------------------------------------------------------------!

      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      call PBC(xa,xb,xAB)
      call PBC(xc,xd,xCD)

      D2AB = (xAB*xAB + yAB*yAB  + zAB*zAB)
      D2CD = (xCD*xCD + yCD*yCD  + zCD*zCD)

      DAB  = dsqrt(D2AB)
      DCD  = dsqrt(D2CD)
     
      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital

          p   = alpha + beta 
          ip  = 1.0d0/p

          mu1 = alpha*beta*ip

          E_AB = dexp(-mu1*D2AB)

          xp  = ( xa * alpha + xb * beta ) * ip
          call bary_center(alpha,xa,beta,xb,ip,xp)
          yp  = ( ya * alpha + yb * beta ) * ip
          zp  = ( za * alpha + zb * beta ) * ip

          xPA = xp - xa 
          yPA = yp - ya
          zPA = zp - za

          xPB = xp - xb 
          yPB = yp - yb
          zPB = zp - zb

          call SSD(xp,xa,xPA)
          call SSD(xp,xb,xPB)

          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital

              q   = gamma + delta 
              iq  = 1.0d0/q
              mu2 = gamma*delta*iq
    
              E_CD = dexp(-mu2*D2CD)

              xq = ( xc * gamma + xd * delta )* iq
              call bary_center(gamma,xc,delta,xd,iq,xq)
              yq = ( yc * gamma + yd * delta )* iq
              zq = ( zc * gamma + zd * delta )* iq

              xQC = xq - xc 
              yQC = yq - yc
              zQC = zq - zc
    
              xQD = xq - xd 
              yQD = yq - yd
              zQD = zq - zd

              xPQ = xp - xq
              yPQ = yp - yq
              zPQ = zp - zq

              call SSD(xq,xc,xQC)
              call SSD(xq,xd,xQD)

              call euc(xp,xq,xPQ)
               
              R2PQ = (xPQ*xPQ + yPQ*yPQ + zPQ*zPQ)

               p_plus_q  = p+q 
              ip_plus_q  = 1.0d0/p_plus_q 

                   popq  = p * ip_plus_q     
                   qopq  = q * ip_plus_q

                      tu = 0.50d0*qopq*ip
                      tv = 0.50d0*popq*iq

                      pq = 0.50d0 * ip_plus_q
                   
               p_times_q = p*q
              ip_times_q = 1.0d0/p_times_q 

              c_boys = p_times_q*ip_plus_q

              const  = (c1*c2*c3*c4) * R2PI52 * dsqrt(ip_plus_q) * ip_times_q * E_AB * E_CD

              F_in   = c_boys*R2PQ

              F0     = Boys_func(0,F_in)
              F1     = Boys_func(1,F_in)
              F2     = Boys_func(2,F_in)
              F3     = Boys_func(3,F_in)
              F4     = Boys_func(4,F_in)

              call ERI_value(o1,o2,o3,o4&
                            &,F0,F1,F2,F3,F4&
                            &,popq,qopq,tu,tv,p,q,pq&
                            &,xPA,yPA,zPA,xPB,yPB,zPB&
                            &,xQC,yQC,zQC,xQD,yQD,zQD&
                            &,xPQ,yPQ,zPQ,value_s)

              value = value + const * value_s
              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_torus