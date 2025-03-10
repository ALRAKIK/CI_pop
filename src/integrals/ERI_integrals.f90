subroutine ERI_integral_4_function(one,two,three,four,value)
      
      use torus_init
      use classification_ERI

      implicit none 
      
      type(ERI_function),intent(in)   :: one , two , three , four 

      integer                         ::  i , j   , k  , l 

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


      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision                :: value_s 

      double precision,intent(out)    :: value 


      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      if (torus) call PBC(xa,xb,xAB)
      if (torus) call PBC(xc,xd,xCD)

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
          if (torus) call bary_center(alpha,xa,beta,xb,ip,xp)
          yp  = ( ya * alpha + yb * beta ) * ip
          zp  = ( za * alpha + zb * beta ) * ip

          xPA = xp - xa 
          yPA = yp - ya
          zPA = zp - za

          xPB = xp - xb 
          yPB = yp - yb
          zPB = zp - zb

          if (torus) call SSD(xp,xa,xPA)
          if (torus) call SSD(xp,xb,xPB)

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
              if (torus) call bary_center(gamma,xc,delta,xd,iq,xq)
              yq = ( yc * gamma + yd * delta )* iq
              zq = ( zc * gamma + zd * delta )* iq

              xQC = xq - xc 
              yQC = yq - yc
              zQC = zq - zc
    
              xQD = xq - xd 
              yQD = yq - yd
              zQD = zq - zd

              if (torus) call SSD(xq,xc,xQC)
              if (torus) call SSD(xq,xd,xQD)

              xPQ = xp - xq
              if (torus) call euc(xp,xq,xPQ)
              yPQ = yp - yq
              zPQ = zp - zq
               
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

end subroutine


subroutine ERI_value(o1,o2,o3,o4&
                    &,F0,F1,F2,F3,F4&
                    &,u,v,tu,tv,p,q,pq&
                    &,xPA,yPA,zPA,xPB,yPB,zPB&
                    &,xQC,yQC,zQC,xQD,yQD,zQD&
                    &,xPQ,yPQ,zPQ,value_s)


      implicit none 

      character(LEN=2),intent(in)  :: o1 , o2 , o3 , o4
      double precision,intent(in)  :: F0 , F1 , F2 , F3 , F4 
      

      double precision,intent(in)  :: u , v , tu , tv , p , q , pq 
      double precision,intent(in)  :: xPA,yPA,zPA,xPB,yPB,zPB
      double precision,intent(in)  :: xQC,yQC,zQC,xQD,yQD,zQD

      double precision,intent(in)  :: xPQ , yPQ , zPQ

      double precision,intent(out) :: value_s

      ! ////////////////////////// s s s s //////////////////////////// !

      if (o1 == "s" .and. o2=="s" .and. o3=="s" .and. o4=="s") then 
          value_s = F0 
      end if 

      ! ////////////////////////// p s s s //////////////////////////// !
   
      if ((o1(:1) == "p" .and.     o2=="s" .and.     o3=="s" .and.     o4=="s").or.&
        & (o1     == "s" .and. o2(:1)=="p" .and.     o3=="s" .and.     o4=="s").or.&
        & (o1     == "s" .and.     o2=="s" .and. o3(:1)=="p" .and.     o4=="s").or.&
        & (o1     == "s" .and.     o2=="s" .and.     o3=="s" .and. o4(:1)=="p")     ) then
          
          call ERI_psss(o1,o2,o3,o4&
                      &,F0,F1,F2,F3&
                      &,u,v,tu,tv,p,q,pq&
                      &,xPA,yPA,zPA,xPB,yPB,zPB&
                      &,xQC,yQC,zQC,xQD,yQD,zQD&
                      &,xPQ,yPQ,zPQ,value_s)

      end if
      
      ! ////////////////////////// p p s s //////////////////////////// !
   
      if ((o1(:1) == "p" .and. o2(:1)=="p" .and.      o3=="s" .and.      o4=="s").or.&
        & (o1(:1) == "p" .and.     o2=="s" .and.  o3(:1)=="p" .and.      o4=="s").or.&
        & (o1(:1) == "p" .and.     o2=="s" .and.      o3=="s" .and.  o4(:1)=="p").or.&
        & (o1     == "s" .and. o2(:1)=="p" .and.  o3(:1)=="p" .and.      o4=="s").or.&
        & (o1     == "s" .and. o2(:1)=="p" .and.      o3=="s" .and.  o4(:1)=="p").or.&
        & (o1     == "s" .and.     o2=="s" .and.  o3(:1)=="p" .and.  o4(:1)=="p")) then
          
          call ERI_ppss(o1,o2,o3,o4&
                      &,F0,F1,F2,F3&
                      &,u,v,tu,tv,p,q,pq&
                      &,xPA,yPA,zPA,xPB,yPB,zPB&
                      &,xQC,yQC,zQC,xQD,yQD,zQD&
                      &,xPQ,yPQ,zPQ,value_s)

      end if

      ! ////////////////////////// p p p s //////////////////////////// !
   
      if ((o1(:1) == "p" .and. o2(:1)=="p" .and.   o3(:1)=="p" .and.      o4=="s").or.&
        & (o1(:1) == "p" .and. o2(:1)=="p" .and.       o3=="s" .and.  o4(:1)=="p").or.&
        & (o1(:1) == "p" .and.     o2=="s" .and.   o3(:1)=="p" .and.  o4(:1)=="p").or.&
        & (o1     == "s" .and. o2(:1)=="p" .and.   o3(:1)=="p" .and.  o4(:1)=="p")) then
          
          call ERI_ppps(o1,o2,o3,o4&
                      &,F0,F1,F2,F3&
                      &,u,v,tu,tv,p,q,pq&
                      &,xPA,yPA,zPA,xPB,yPB,zPB&
                      &,xQC,yQC,zQC,xQD,yQD,zQD&
                      &,xPQ,yPQ,zPQ,value_s)

      end if

            ! ////////////////////////// p p p s //////////////////////////// !
   
      if (o1(:1) == "p" .and. o2(:1)=="p" .and.   o3(:1)=="p" .and.  o4(:1)=="p") then
          
          call ERI_pppp(o1,o2,o3,o4&
                      &,F0,F1,F2,F3,F4&
                      &,u,v,tu,tv,p,q,pq&
                      &,xPA,yPA,zPA,xPB,yPB,zPB&
                      &,xQC,yQC,zQC,xQD,yQD,zQD&
                      &,xPQ,yPQ,zPQ,value_s)

      end if






end subroutine


subroutine ERI_psss(o1,o2,o3,o4&
                  &,F0,F1,F2,F3&
                  &,u,v,tu,tv,p,q,pq&
                  &,xPA,yPA,zPA,xPB,yPB,zPB&
                  &,xQC,yQC,zQC,xQD,yQD,zQD&
                  &,xPQ,yPQ,zPQ,value)


      implicit none 


      character(LEN=2),intent(in)  :: o1 , o2 , o3 , o4
      double precision,intent(in)  :: F0 , F1 , F2 , F3 
      

      double precision,intent(in)  :: xPA,yPA,zPA,xPB,yPB,zPB
      double precision,intent(in)  :: xQC,yQC,zQC,xQD,yQD,zQD

      double precision,intent(in)  :: xPQ , yPQ , zPQ

      double precision,intent(in)  :: u , v , tu , tv , p ,q  , pq 
      
      double precision,intent(out) :: value

      ! ////////////////////////// p s s s //////////////////////////// !

      !                            p_x s s s                            !

      if (o1 == "px" .and. o2=="s" .and. o3=="s" .and. o4=="s") then 
        value = -F1*v*xPQ + xPA*F0  
      end if

      !                            p_y s s s                            !

      if (o1 == "py" .and. o2=="s" .and. o3=="s" .and. o4=="s") then 
        value = -F1*v*yPQ + yPA*F0 
      end if

      !                            p_z s s s                            !

      if (o1 == "pz" .and. o2=="s" .and. o3=="s" .and. o4=="s") then 
        value = -F1*v*zPQ + zPA * F0 
      end if

      ! ////////////////////////// s p s s //////////////////////////// !

      !                            s p_x s s                            !

      if (o1 == "s" .and. o2=="px" .and. o3=="s" .and. o4=="s") then 
        value =  -F1*v*xPQ + xPB * F0 
      end if

      !                            s p_y s s                            !

      if (o1 == "s" .and. o2=="py" .and. o3=="s" .and. o4=="s") then 
        value = -F1*v*yPQ + yPB * F0 
      end if

      !                            s p_z s s                            !

      if (o1 == "s" .and. o2=="pz" .and. o3=="s" .and. o4=="s") then 
        value = -F1*v*zPQ + zPB*F0 
      end if

      ! ////////////////////////// s s p s //////////////////////////// !

      !                            s s p_x s                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="px" .and. o4=="s") then 
          value = F1*u*xPQ + xQC * F0 
      end if

      !                            s s p_y s                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="py" .and. o4=="s") then 
          value = F1*u*yPQ + yQC * F0 
      end if

      !                            s s p_z s                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="pz" .and. o4=="s") then 
        value = F1*u*zPQ + zQC * F0 
      end if


      ! ////////////////////////// s s s p //////////////////////////// !

      !                            s s s p_x                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="s" .and. o4=="px") then 
          value = F1*u*xPQ + xQD * F0 
      end if

      !                            s s s p_y                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="s" .and. o4=="py") then 
        value = F1*u*yPQ + yQD * F0 
      end if

      !                            s s s p_z                            !

      if (o1 == "s" .and. o2=="s" .and. o3=="s" .and. o4=="pz") then 
        value = F1*u*zPQ + zQD * F0 
      end if


end subroutine


subroutine ERI_ppss(o1,o2,o3,o4&
             &,F0,F1,F2,F3&
             &,u,v,tu,tv,p,q,pq&
             &,xPA,yPA,zPA,xPB,yPB,zPB&
             &,xQC,yQC,zQC,xQD,yQD,zQD&
             &,xPQ,yPQ,zPQ,value)


      implicit none 


      character(LEN=2),intent(in)  :: o1 , o2 , o3 , o4
      double precision,intent(in)  :: F0,F1,F2,F3 
      

      double precision,intent(in)  :: xPA,yPA,zPA,xPB,yPB,zPB
      double precision,intent(in)  :: xQC,yQC,zQC,xQD,yQD,zQD

      double precision,intent(in)  :: xPQ , yPQ , zPQ

      double precision,intent(in)  :: u , v , tu , tv , p , q ,pq 
      
      double precision,intent(out) :: value

      ! Combination 1: p_x p_x s s
      if (o1 == "px" .and. o2=="px" .and. o3=="s" .and. o4=="s") then 
        value = F1*(-tu + v*(-xPA*xPQ - xPB*xPQ)) + F2*v**2*xPQ**2 + (xPA*xPB + 1/(2*p)) * F0 
      end if

      ! Combination 2: p_x s p_x s
      if (o1 == "px" .and. o2 == "s" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq + u*xPA*xPQ - v*xPQ*xQC) - F2*u*v*xPQ**2 + xPA*xQC * F0 
      end if

      ! Combination 3: p_x s s p_x
      if (o1 == "px" .and. o2 == "s" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq + u*xPA*xPQ - v*xPQ*xQD) - F2*u*v*xPQ**2 + xPA*xQD * F0
      end if

      ! Combination 4: p_x p_y s s
      if (o1 == "px" .and. o2 == "py" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-xPA*yPQ - xPQ*yPB) + F2*v**2*xPQ*yPQ + xPA*yPB * F0
      end if
    
      ! Combination 5: p_x s p_y s
      if (o1 == "px" .and. o2 == "s" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*xPA*yPQ - v*xPQ*yQC) - F2*u*v*xPQ*yPQ + xPA*yQC * F0
      end if
    
      ! Combination 6: p_x s s p_y
      if (o1 == "px" .and. o2 == "s" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*xPA*yPQ - v*xPQ*yQD) - F2*u*v*xPQ*yPQ + xPA*yQD * F0
      end if
    
      ! Combination 7: p_x p_z s s
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-xPA*zPQ - xPQ*zPB) + F2*v**2*xPQ*zPQ + xPA*zPB * F0
      end if
    
      ! Combination 8: p_x s p_z s
      if (o1 == "px" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*xPA*zPQ - v*xPQ*zQC) - F2*u*v*xPQ*zPQ + xPA*zQC * F0
      end if

      ! Combination 9: p_x s s p_z
      if (o1 == "px" .and. o2 == "s" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*xPA*zPQ - v*xPQ*zQD) - F2*u*v*xPQ*zPQ + xPA*zQD * F0
      end if
    
      ! Combination 10: s p_x p_x s
      if (o1 == "s" .and. o2 == "px" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq + u*xPB*xPQ - v*xPQ*xQC) - F2*u*v*xPQ**2 + xPB*xQC * F0
      end if
    
      ! Combination 11: s p_x s p_x
      if (o1 == "s" .and. o2 == "px" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq + u*xPB*xPQ - v*xPQ*xQD) - F2*u*v*xPQ**2 + xPB*xQD * F0
      end if
    
      ! Combination 12: p_y p_x s s
      if (o1 == "py" .and. o2 == "px" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-xPB*yPQ - xPQ*yPA) + F2*v**2*xPQ*yPQ + xPB*yPA * F0
      end if
    
      ! Combination 13: s p_x p_y s
      if (o1 == "s" .and. o2 == "px" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*xPB*yPQ - v*xPQ*yQC) - F2*u*v*xPQ*yPQ + xPB*yQC * F0
          end if
        
      ! Combination 14: s p_x s p_y
      if (o1 == "s" .and. o2 == "px" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*xPB*yPQ - v*xPQ*yQD) - F2*u*v*xPQ*yPQ + xPB*yQD * F0
      end if
        
          ! Combination 15: p_z p_x s s
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-xPB*zPQ - xPQ*zPA) + F2*v**2*xPQ*zPQ + xPB*zPA * F0
      end if
    
      ! Combination 16: s p_x p_z s
      if (o1 == "s" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*xPB*zPQ - v*xPQ*zQC) - F2*u*v*xPQ*zPQ + xPB*zQC * F0
      end if
    
      ! Combination 17: s p_x s p_z
      if (o1 == "s" .and. o2 == "px" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*xPB*zPQ - v*xPQ*zQD) - F2*u*v*xPQ*zPQ + xPB*zQD * F0
      end if
    
      ! Combination 18: s s p_x p_x
      if (o1 == "s" .and. o2 == "s" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv + u*(xPQ*xQC + xPQ*xQD)) + F2*u**2*xPQ**2 + (xQC*xQD + 1/(2*q)) * F0
      end if
    
      ! Combination 19: p_y s p_x s
      if (o1 == "py" .and. o2 == "s" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*yPA - v*xQC*yPQ) - F2*u*v*xPQ*yPQ + xQC*yPA * F0
      end if
    
      ! Combination 20: s p_y p_x s
      if (o1 == "s" .and. o2 == "py" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*yPB - v*xQC*yPQ) - F2*u*v*xPQ*yPQ + xQC*yPB * F0
      end if
    
      ! Combination 21: s s p_x p_y
      if (o1 == "s" .and. o2 == "s" .and. o3 == "px" .and. o4 == "py") then
        value = F1*u*(xPQ*yQD + xQC*yPQ) + F2*u**2*xPQ*yPQ + xQC*yQD * F0
      end if
    
      ! Combination 22: p_z s p_x s
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*zPA - v*xQC*zPQ) - F2*u*v*xPQ*zPQ + xQC*zPA * F0
      end if
    
      ! Combination 23: s p_z p_x s
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*zPB - v*xQC*zPQ) - F2*u*v*xPQ*zPQ + xQC*zPB * F0
      end if
    
      ! Combination 24: s s p_x p_z
      if (o1 == "s" .and. o2 == "s" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*u*(xPQ*zQD + xQC*zPQ) + F2*u**2*xPQ*zPQ + xQC*zQD * F0
      end if
    
      ! Combination 25: p_y s s p_x
      if (o1 == "py" .and. o2 == "s" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*yPA - v*xQD*yPQ) - F2*u*v*xPQ*yPQ + xQD*yPA * F0
      end if
    
      ! Combination 26: s p_y s p_x
      if (o1 == "s" .and. o2 == "py" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*yPB - v*xQD*yPQ) - F2*u*v*xPQ*yPQ + xQD*yPB * F0
      end if
    
      ! Combination 27: s s p_y p_x
      if (o1 == "s" .and. o2 == "s" .and. o3 == "py" .and. o4 == "px") then
        value = F1*u*(xPQ*yQC + xQD*yPQ) + F2*u**2*xPQ*yPQ + xQD*yQC * F0
      end if
    
      ! Combination 28: p_z s s p_x
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*zPA - v*xQD*zPQ) - F2*u*v*xPQ*zPQ + xQD*zPA * F0
      end if
    
      ! Combination 29: s p_z s p_x
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*zPB - v*xQD*zPQ) - F2*u*v*xPQ*zPQ + xQD*zPB * F0
      end if
    
      ! Combination 30: s s p_z p_x
      if (o1 == "s" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*u*(xPQ*zQC + xQD*zPQ) + F2*u**2*xPQ*zPQ + xQD*zQC * F0
      end if
    
      ! Combination 31: p_y p_y s s
      if (o1 == "py" .and. o2 == "py" .and. o3 == "s" .and. o4 == "s") then
        value = F1*(-tu + v*(-yPA*yPQ - yPB*yPQ)) + F2*v**2*yPQ**2 + (yPA*yPB + 1/(2*p)) * F0
      end if
    
      ! Combination 32: p_y s p_y s
      if (o1 == "py" .and. o2 == "s" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq + u*yPA*yPQ - v*yPQ*yQC) - F2*u*v*yPQ**2 + yPA*yQC * F0
      end if
    
      ! Combination 33: p_y s s p_y
      if (o1 == "py" .and. o2 == "s" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq + u*yPA*yPQ - v*yPQ*yQD) - F2*u*v*yPQ**2 + yPA*yQD * F0
      end if
    
      ! Combination 34: p_y p_z s s
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-yPA*zPQ - yPQ*zPB) + F2*v**2*yPQ*zPQ + yPA*zPB * F0
      end if
    
      ! Combination 35: p_y s p_z s
      if (o1 == "py" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*yPA*zPQ - v*yPQ*zQC) - F2*u*v*yPQ*zPQ + yPA*zQC * F0
      end if
    
      ! Combination 36: p_y s s p_z
      if (o1 == "py" .and. o2 == "s" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*yPA*zPQ - v*yPQ*zQD) - F2*u*v*yPQ*zPQ + yPA*zQD * F0
      end if
    
      ! Combination 37: s p_y p_y s
      if (o1 == "s" .and. o2 == "py" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq + u*yPB*yPQ - v*yPQ*yQC) - F2*u*v*yPQ**2 + yPB*yQC * F0
      end if
    
      ! Combination 38: s p_y s p_y
      if (o1 == "s" .and. o2 == "py" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq + u*yPB*yPQ - v*yPQ*yQD) - F2*u*v*yPQ**2 + yPB*yQD * F0
      end if
    
      ! Combination 39: p_z p_y s s
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "s" .and. o4 == "s") then
        value = F1*v*(-yPB*zPQ - yPQ*zPA) + F2*v**2*yPQ*zPQ + yPB*zPA * F0
      end if
    
      ! Combination 40: s p_y p_z s
      if (o1 == "s" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*yPB*zPQ - v*yPQ*zQC) - F2*u*v*yPQ*zPQ + yPB*zQC * F0
      end if
    
      ! Combination 41: s p_y s p_z
      if (o1 == "s" .and. o2 == "py" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*yPB*zPQ - v*yPQ*zQD) - F2*u*v*yPQ*zPQ + yPB*zQD * F0
      end if
    
      ! Combination 42: s s p_y p_y
      if (o1 == "s" .and. o2 == "s" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv + u*(yPQ*yQC + yPQ*yQD)) + F2*u**2*yPQ**2 + (yQC*yQD + 1/(2*q)) * F0
      end if
    
      ! Combination 43: p_z s p_y s
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*yPQ*zPA - v*yQC*zPQ) - F2*u*v*yPQ*zPQ + yQC*zPA * F0
      end if
    
      ! Combination 44: s p_z p_y s
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*yPQ*zPB - v*yQC*zPQ) - F2*u*v*yPQ*zPQ + yQC*zPB * F0
      end if
    
      ! Combination 45: s s p_y p_z
      if (o1 == "s" .and. o2 == "s" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*u*(yPQ*zQD + yQC*zPQ) + F2*u**2*yPQ*zPQ + yQC*zQD * F0
      end if
    
      ! Combination 46: p_z s s p_y
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*yPQ*zPA - v*yQD*zPQ) - F2*u*v*yPQ*zPQ + yQD*zPA * F0
      end if
    
      ! Combination 47: s p_z s p_y
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*yPQ*zPB - v*yQD*zPQ) - F2*u*v*yPQ*zPQ + yQD*zPB * F0
      end if
    
      ! Combination 48: s s p_z p_y
      if (o1 == "s" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*u*(yPQ*zQC + yQD*zPQ) + F2*u**2*yPQ*zPQ + yQD*zQC * F0
      end if
    
      ! Combination 49: p_z p_z s s
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "s") then
        value = F1*(-tu + v*(-zPA*zPQ - zPB*zPQ)) + F2*v**2*zPQ**2 + (zPA*zPB + 1/(2*p)) * F0
      end if
    
      ! Combination 50: p_z s p_z s
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq + u*zPA*zPQ - v*zPQ*zQC) - F2*u*v*zPQ**2 + zPA*zQC * F0
      end if
    
      ! Combination 51: p_z s s p_z
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq + u*zPA*zPQ - v*zPQ*zQD) - F2*u*v*zPQ**2 + zPA*zQD * F0
      end if
    
      ! Combination 52: s p_z p_z s
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq + u*zPB*zPQ - v*zPQ*zQC) - F2*u*v*zPQ**2 + zPB*zQC * F0
      end if
    
      ! Combination 53: s p_z s p_z
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq + u*zPB*zPQ - v*zPQ*zQD) - F2*u*v*zPQ**2 + zPB*zQD * F0
      end if
    
      ! Combination 54: s s p_z p_z
      if (o1 == "s" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv + u*(zPQ*zQC + zPQ*zQD)) + F2*u**2*zPQ**2 + (zQC*zQD + 1/(2*q)) * F0
      end if

end subroutine


subroutine ERI_ppps(o1,o2,o3,o4&
                  &,F0,F1,F2,F3&
                  &,u,v,tu,tv,p,q,pq&
                  &,xPA,yPA,zPA,xPB,yPB,zPB&
                  &,xQC,yQC,zQC,xQD,yQD,zQD&
                  &,xPQ,yPQ,zPQ,value)


      implicit none 


      character(LEN=2),intent(in)  :: o1 , o2 , o3 , o4
      double precision,intent(in)  :: F0,F1,F2,F3
      

      double precision,intent(in)  :: xPA,yPA,zPA,xPB,yPB,zPB
      double precision,intent(in)  :: xQC,yQC,zQC,xQD,yQD,zQD

      double precision,intent(in)  :: xPQ , yPQ , zPQ

      double precision,intent(in)  :: u , v , tu , tv , p , q , pq
      
      double precision,intent(out) :: value
      
      ! ////////////////////////// p p p s //////////////////////////// !

      ! Combination 1: p_x p_x p_x s
      if (o1 == "px" .and. o2 == "px" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq*(xPA + xPB) - tu*xQC + u*(xPA*xPB*xPQ + xPQ/(2*p)) + v*(-xPA*xPQ*xQC - xPB*xPQ*xQC)) + F2*(-2*pq*v*xPQ + u*(-tu*xPQ + v*(-xPA*xPQ**2 - xPB*xPQ**2)) + v**2*xPQ**2*xQC) + F3*u*v**2*xPQ**3 + (xPA*xPB*xQC + xQC/(2*p)) * F0
      end if

      ! Combination 2: p_x p_x s p_x
      if (o1 == "px" .and. o2 == "px" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq*(xPA + xPB) - tu*xQD + u*(xPA*xPB*xPQ + xPQ/(2*p)) + v*(-xPA*xPQ*xQD - xPB*xPQ*xQD)) + F2*(-2*pq*v*xPQ + u*(-tu*xPQ + v*(-xPA*xPQ**2 - xPB*xPQ**2)) + v**2*xPQ**2*xQD) + F3*u*v**2*xPQ**3 + (xPA*xPB*xQD + xQD/(2*p)) * F0
      end if
    
      ! Combination 3: p_x p_x p_y s
      if (o1 == "px" .and. o2 == "px" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(-tu*yQC + u*(xPA*xPB*yPQ + yPQ/(2*p)) + v*(-xPA*xPQ*yQC - xPB*xPQ*yQC)) + F2*(u*(-tu*yPQ + v*(-xPA*xPQ*yPQ - xPB*xPQ*yPQ)) + v**2*xPQ**2*yQC) + F3*u*v**2*xPQ**2*yPQ + (xPA*xPB*yQC + yQC/(2*p)) * F0
      end if
    
      ! Combination 4: p_x p_x s p_y
      if (o1 == "px" .and. o2 == "px" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(-tu*yQD + u*(xPA*xPB*yPQ + yPQ/(2*p)) + v*(-xPA*xPQ*yQD - xPB*xPQ*yQD)) + F2*(u*(-tu*yPQ + v*(-xPA*xPQ*yPQ - xPB*xPQ*yPQ)) + v**2*xPQ**2*yQD) + F3*u*v**2*xPQ**2*yPQ + (xPA*xPB*yQD + yQD/(2*p)) * F0
      end if
    
      ! Combination 5: p_x p_x p_z s
      if (o1 == "px" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(-tu*zQC + u*(xPA*xPB*zPQ + zPQ/(2*p)) + v*(-xPA*xPQ*zQC - xPB*xPQ*zQC)) + F2*(u*(-tu*zPQ + v*(-xPA*xPQ*zPQ - xPB*xPQ*zPQ)) + v**2*xPQ**2*zQC) + F3*u*v**2*xPQ**2*zPQ + (xPA*xPB*zQC + zQC/(2*p)) * F0
      end if
    
      ! Combination 6: p_x p_x s p_z
      if (o1 == "px" .and. o2 == "px" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(-tu*zQD + u*(xPA*xPB*zPQ + zPQ/(2*p)) + v*(-xPA*xPQ*zQD - xPB*xPQ*zQD)) + F2*(u*(-tu*zPQ + v*(-xPA*xPQ*zPQ - xPB*xPQ*zPQ)) + v**2*xPQ**2*zQD) + F3*u*v**2*xPQ**2*zPQ + (xPA*xPB*zQD + zQD/(2*p)) * F0
      end if
    
      ! Combination 7: p_x s p_x p_x
      if (o1 == "px" .and. o2 == "s" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC + xQD) - tv*xPA + u*(xPA*xPQ*xQC + xPA*xPQ*xQD) + v*(-xPQ*xQC*xQD - xPQ/(2*q))) + F2*(tv*v*xPQ + u**2*xPA*xPQ**2 + u*(2*pq*xPQ + v*(-xPQ**2*xQC - xPQ**2*xQD))) - F3*u**2*v*xPQ**3 + (xPA*xQC*xQD + xPA/(2*q)) * F0
      end if
    
      ! Combination 8: p_x p_y p_x s
      if (o1 == "px" .and. o2 == "py" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq*yPB + u*xPA*xPQ*yPB + v*(-xPA*xQC*yPQ - xPQ*xQC*yPB)) + F2*(v**2*xPQ*xQC*yPQ + v*(-pq*yPQ + u*(-xPA*xPQ*yPQ - xPQ**2*yPB))) + F3*u*v**2*xPQ**2*yPQ + xPA*xQC*yPB * F0
      end if
    
      ! Combination 9: p_x s p_x p_y
      if (o1 == "px" .and. o2 == "s" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*yQD + u*(xPA*xPQ*yQD + xPA*xQC*yPQ) - v*xPQ*xQC*yQD) + F2*(u**2*xPA*xPQ*yPQ + u*(pq*yPQ + v*(-xPQ**2*yQD - xPQ*xQC*yPQ))) - F3*u**2*v*xPQ**2*yPQ + xPA*xQC*yQD * F0
      end if
    
      ! Combination 10: p_x p_z p_x s
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq*zPB + u*xPA*xPQ*zPB + v*(-xPA*xQC*zPQ - xPQ*xQC*zPB)) + F2*(v**2*xPQ*xQC*zPQ + v*(-pq*zPQ + u*(-xPA*xPQ*zPQ - xPQ**2*zPB))) + F3*u*v**2*xPQ**2*zPQ + xPA*xQC*zPB * F0
      end if
    
      ! Combination 11: p_x s p_x p_z
      if (o1 == "px" .and. o2 == "s" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*zQD + u*(xPA*xPQ*zQD + xPA*xQC*zPQ) - v*xPQ*xQC*zQD) + F2*(u**2*xPA*xPQ*zPQ + u*(pq*zPQ + v*(-xPQ**2*zQD - xPQ*xQC*zPQ))) - F3*u**2*v*xPQ**2*zPQ + xPA*xQC*zQD * F0
      end if
    
      ! Combination 12: p_x p_y s p_x
      if (o1 == "px" .and. o2 == "py" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq*yPB + u*xPA*xPQ*yPB + v*(-xPA*xQD*yPQ - xPQ*xQD*yPB)) + F2*(v**2*xPQ*xQD*yPQ + v*(-pq*yPQ + u*(-xPA*xPQ*yPQ - xPQ**2*yPB))) + F3*u*v**2*xPQ**2*yPQ + xPA*xQD*yPB * F0
      end if
    
      ! Combination 13: p_x s p_y p_x
      if (o1 == "px" .and. o2 == "s" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*yQC + u*(xPA*xPQ*yQC + xPA*xQD*yPQ) - v*xPQ*xQD*yQC) + F2*(u**2*xPA*xPQ*yPQ + u*(pq*yPQ + v*(-xPQ**2*yQC - xPQ*xQD*yPQ))) - F3*u**2*v*xPQ**2*yPQ + xPA*xQD*yQC * F0
      end if
    
      ! Combination 14: p_x p_z s p_x
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq*zPB + u*xPA*xPQ*zPB + v*(-xPA*xQD*zPQ - xPQ*xQD*zPB)) + F2*(v**2*xPQ*xQD*zPQ + v*(-pq*zPQ + u*(-xPA*xPQ*zPQ - xPQ**2*zPB))) + F3*u*v**2*xPQ**2*zPQ + xPA*xQD*zPB * F0
      end if
    
      ! Combination 15: p_x s p_z p_x
      if (o1 == "px" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*zQC + u*(xPA*xPQ*zQC + xPA*xQD*zPQ) - v*xPQ*xQD*zQC) + F2*(u**2*xPA*xPQ*zPQ + u*(pq*zPQ + v*(-xPQ**2*zQC - xPQ*xQD*zPQ))) - F3*u**2*v*xPQ**2*zPQ + xPA*xQD*zQC * F0
      end if
    
      ! Combination 16: p_x p_y p_y s
      if (o1 == "px" .and. o2 == "py" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq*xPA + u*xPA*yPB*yPQ + v*(-xPA*yPQ*yQC - xPQ*yPB*yQC)) + F2*(v**2*xPQ*yPQ*yQC + v*(-pq*xPQ + u*(-xPA*yPQ**2 - xPQ*yPB*yPQ))) + F3*u*v**2*xPQ*yPQ**2 + xPA*yPB*yQC * F0
      end if
    
      ! Combination 17: p_x p_y s p_y
      if (o1 == "px" .and. o2 == "py" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq*xPA + u*xPA*yPB*yPQ + v*(-xPA*yPQ*yQD - xPQ*yPB*yQD)) + F2*(v**2*xPQ*yPQ*yQD + v*(-pq*xPQ + u*(-xPA*yPQ**2 - xPQ*yPB*yPQ))) + F3*u*v**2*xPQ*yPQ**2 + xPA*yPB*yQD * F0
      end if
    
      ! Combination 18: p_x p_y p_z s
      if (o1 == "px" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*xPA*yPB*zPQ + v*(-xPA*yPQ*zQC - xPQ*yPB*zQC)) + F2*(u*v*(-xPA*yPQ*zPQ - xPQ*yPB*zPQ) + v**2*xPQ*yPQ*zQC) + F3*u*v**2*xPQ*yPQ*zPQ + xPA*yPB*zQC * F0
      end if
    
      ! Combination 19: p_x p_y s p_z
      if (o1 == "px" .and. o2 == "py" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*xPA*yPB*zPQ + v*(-xPA*yPQ*zQD - xPQ*yPB*zQD)) + F2*(u*v*(-xPA*yPQ*zPQ - xPQ*yPB*zPQ) + v**2*xPQ*yPQ*zQD) + F3*u*v**2*xPQ*yPQ*zPQ + xPA*yPB*zQD * F0
      end if
    
      ! Combination 20: p_x s p_y p_y
      if (o1 == "px" .and. o2 == "s" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*xPA + u*(xPA*yPQ*yQC + xPA*yPQ*yQD) + v*(-xPQ*yQC*yQD - xPQ/(2*q))) + F2*(u**2*xPA*yPQ**2 + v*(tv*xPQ + u*(-xPQ*yPQ*yQC - xPQ*yPQ*yQD))) - F3*u**2*v*xPQ*yPQ**2 + (xPA*yQC*yQD + xPA/(2*q)) * F0
      end if
    
      ! Combination 21: p_x p_z p_y s
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*xPA*yPQ*zPB + v*(-xPA*yQC*zPQ - xPQ*yQC*zPB)) + F2*(u*v*(-xPA*yPQ*zPQ - xPQ*yPQ*zPB) + v**2*xPQ*yQC*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xPA*yQC*zPB * F0
      end if
    
      ! Combination 22: p_x s p_y p_z
      if (o1 == "px" .and. o2 == "s" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(u*(xPA*yPQ*zQD + xPA*yQC*zPQ) - v*xPQ*yQC*zQD) + F2*(u**2*xPA*yPQ*zPQ + u*v*(-xPQ*yPQ*zQD - xPQ*yQC*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xPA*yQC*zQD * F0
      end if
    
      ! Combination 23: p_x p_z s p_y
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*xPA*yPQ*zPB + v*(-xPA*yQD*zPQ - xPQ*yQD*zPB)) + F2*(u*v*(-xPA*yPQ*zPQ - xPQ*yPQ*zPB) + v**2*xPQ*yQD*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xPA*yQD*zPB * F0
      end if
    
      ! Combination 24: p_x s p_z p_y
      if (o1 == "px" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(u*(xPA*yPQ*zQC + xPA*yQD*zPQ) - v*xPQ*yQD*zQC) + F2*(u**2*xPA*yPQ*zPQ + u*v*(-xPQ*yPQ*zQC - xPQ*yQD*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xPA*yQD*zQC * F0
      end if
    
      ! Combination 25: p_x p_z p_z s
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq*xPA + u*xPA*zPB*zPQ + v*(-xPA*zPQ*zQC - xPQ*zPB*zQC)) + F2*(v**2*xPQ*zPQ*zQC + v*(-pq*xPQ + u*(-xPA*zPQ**2 - xPQ*zPB*zPQ))) + F3*u*v**2*xPQ*zPQ**2 + xPA*zPB*zQC * F0
      end if
    
      ! Combination 26: p_x p_z s p_z
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq*xPA + u*xPA*zPB*zPQ + v*(-xPA*zPQ*zQD - xPQ*zPB*zQD)) + F2*(v**2*xPQ*zPQ*zQD + v*(-pq*xPQ + u*(-xPA*zPQ**2 - xPQ*zPB*zPQ))) + F3*u*v**2*xPQ*zPQ**2 + xPA*zPB*zQD * F0
      end if
    
      ! Combination 27: p_x s p_z p_z
      if (o1 == "px" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*xPA + u*(xPA*zPQ*zQC + xPA*zPQ*zQD) + v*(-xPQ*zQC*zQD - xPQ/(2*q))) + F2*(u**2*xPA*zPQ**2 + v*(tv*xPQ + u*(-xPQ*zPQ*zQC - xPQ*zPQ*zQD))) - F3*u**2*v*xPQ*zPQ**2 + (xPA*zQC*zQD + xPA/(2*q)) * F0
      end if
    
      ! Combination 28: s p_x p_x p_x
      if (o1 == "s" .and. o2 == "px" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC + xQD) - tv*xPB + u*(xPB*xPQ*xQC + xPB*xPQ*xQD) + v*(-xPQ*xQC*xQD - xPQ/(2*q))) + F2*(tv*v*xPQ + u**2*xPB*xPQ**2 + u*(2*pq*xPQ + v*(-xPQ**2*xQC - xPQ**2*xQD))) - F3*u**2*v*xPQ**3 + (xPB*xQC*xQD + xPB/(2*q)) * F0
      end if
    
      ! Combination 29: p_y p_x p_x s
      if (o1 == "py" .and. o2 == "px" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq*yPA + u*xPB*xPQ*yPA + v*(-xPB*xQC*yPQ - xPQ*xQC*yPA)) + F2*(v**2*xPQ*xQC*yPQ + v*(-pq*yPQ + u*(-xPB*xPQ*yPQ - xPQ**2*yPA))) + F3*u*v**2*xPQ**2*yPQ + xPB*xQC*yPA * F0
      end if
    
      ! Combination 30: s p_x p_x p_y
      if (o1 == "s" .and. o2 == "px" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*yQD + u*(xPB*xPQ*yQD + xPB*xQC*yPQ) - v*xPQ*xQC*yQD) + F2*(u**2*xPB*xPQ*yPQ + u*(pq*yPQ + v*(-xPQ**2*yQD - xPQ*xQC*yPQ))) - F3*u**2*v*xPQ**2*yPQ + xPB*xQC*yQD * F0
      end if
    
      ! Combination 31: p_z p_x p_x s
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(pq*zPA + u*xPB*xPQ*zPA + v*(-xPB*xQC*zPQ - xPQ*xQC*zPA)) + F2*(v**2*xPQ*xQC*zPQ + v*(-pq*zPQ + u*(-xPB*xPQ*zPQ - xPQ**2*zPA))) + F3*u*v**2*xPQ**2*zPQ + xPB*xQC*zPA * F0
      end if
    
      ! Combination 32: s p_x p_x p_z
      if (o1 == "s" .and. o2 == "px" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*zQD + u*(xPB*xPQ*zQD + xPB*xQC*zPQ) - v*xPQ*xQC*zQD) + F2*(u**2*xPB*xPQ*zPQ + u*(pq*zPQ + v*(-xPQ**2*zQD - xPQ*xQC*zPQ))) - F3*u**2*v*xPQ**2*zPQ + xPB*xQC*zQD * F0
      end if
    
      ! Combination 33: p_y p_x s p_x
      if (o1 == "py" .and. o2 == "px" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq*yPA + u*xPB*xPQ*yPA + v*(-xPB*xQD*yPQ - xPQ*xQD*yPA)) + F2*(v**2*xPQ*xQD*yPQ + v*(-pq*yPQ + u*(-xPB*xPQ*yPQ - xPQ**2*yPA))) + F3*u*v**2*xPQ**2*yPQ + xPB*xQD*yPA * F0
      end if
    
      ! Combination 34: s p_x p_y p_x
      if (o1 == "s" .and. o2 == "px" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*yQC + u*(xPB*xPQ*yQC + xPB*xQD*yPQ) - v*xPQ*xQD*yQC) + F2*(u**2*xPB*xPQ*yPQ + u*(pq*yPQ + v*(-xPQ**2*yQC - xPQ*xQD*yPQ))) - F3*u**2*v*xPQ**2*yPQ + xPB*xQD*yQC * F0
      end if
    
      ! Combination 35: p_z p_x s p_x
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(pq*zPA + u*xPB*xPQ*zPA + v*(-xPB*xQD*zPQ - xPQ*xQD*zPA)) + F2*(v**2*xPQ*xQD*zPQ + v*(-pq*zPQ + u*(-xPB*xPQ*zPQ - xPQ**2*zPA))) + F3*u*v**2*xPQ**2*zPQ + xPB*xQD*zPA * F0
      end if
    
      ! Combination 36: s p_x p_z p_x
      if (o1 == "s" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*zQC + u*(xPB*xPQ*zQC + xPB*xQD*zPQ) - v*xPQ*xQD*zQC) + F2*(u**2*xPB*xPQ*zPQ + u*(pq*zPQ + v*(-xPQ**2*zQC - xPQ*xQD*zPQ))) - F3*u**2*v*xPQ**2*zPQ + xPB*xQD*zQC * F0
      end if
    
      ! Combination 37: p_y p_x p_y s
      if (o1 == "py" .and. o2 == "px" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq*xPB + u*xPB*yPA*yPQ + v*(-xPB*yPQ*yQC - xPQ*yPA*yQC)) + F2*(v**2*xPQ*yPQ*yQC + v*(-pq*xPQ + u*(-xPB*yPQ**2 - xPQ*yPA*yPQ))) + F3*u*v**2*xPQ*yPQ**2 + xPB*yPA*yQC * F0
      end if
    
      ! Combination 38: p_y p_x s p_y
      if (o1 == "py" .and. o2 == "px" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq*xPB + u*xPB*yPA*yPQ + v*(-xPB*yPQ*yQD - xPQ*yPA*yQD)) + F2*(v**2*xPQ*yPQ*yQD + v*(-pq*xPQ + u*(-xPB*yPQ**2 - xPQ*yPA*yPQ))) + F3*u*v**2*xPQ*yPQ**2 + xPB*yPA*yQD * F0
      end if
    
      ! Combination 39: p_y p_x p_z s
      if (o1 == "py" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(u*xPB*yPA*zPQ + v*(-xPB*yPQ*zQC - xPQ*yPA*zQC)) + F2*(u*v*(-xPB*yPQ*zPQ - xPQ*yPA*zPQ) + v**2*xPQ*yPQ*zQC) + F3*u*v**2*xPQ*yPQ*zPQ + xPB*yPA*zQC * F0
      end if
    
      ! Combination 40: p_y p_x s p_z
      if (o1 == "py" .and. o2 == "px" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(u*xPB*yPA*zPQ + v*(-xPB*yPQ*zQD - xPQ*yPA*zQD)) + F2*(u*v*(-xPB*yPQ*zPQ - xPQ*yPA*zPQ) + v**2*xPQ*yPQ*zQD) + F3*u*v**2*xPQ*yPQ*zPQ + xPB*yPA*zQD * F0
      end if
    
      ! Combination 41: s p_x p_y p_y
      if (o1 == "s" .and. o2 == "px" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*xPB + u*(xPB*yPQ*yQC + xPB*yPQ*yQD) + v*(-xPQ*yQC*yQD - xPQ/(2*q))) + F2*(u**2*xPB*yPQ**2 + v*(tv*xPQ + u*(-xPQ*yPQ*yQC - xPQ*yPQ*yQD))) - F3*u**2*v*xPQ*yPQ**2 + (xPB*yQC*yQD + xPB/(2*q)) * F0
      end if
    
      ! Combination 42: p_z p_x p_y s
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(u*xPB*yPQ*zPA + v*(-xPB*yQC*zPQ - xPQ*yQC*zPA)) + F2*(u*v*(-xPB*yPQ*zPQ - xPQ*yPQ*zPA) + v**2*xPQ*yQC*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xPB*yQC*zPA * F0
      end if
    
      ! Combination 43: s p_x p_y p_z
      if (o1 == "s" .and. o2 == "px" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(u*(xPB*yPQ*zQD + xPB*yQC*zPQ) - v*xPQ*yQC*zQD) + F2*(u**2*xPB*yPQ*zPQ + u*v*(-xPQ*yPQ*zQD - xPQ*yQC*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xPB*yQC*zQD * F0
      end if
    
      ! Combination 44: p_z p_x s p_y
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(u*xPB*yPQ*zPA + v*(-xPB*yQD*zPQ - xPQ*yQD*zPA)) + F2*(u*v*(-xPB*yPQ*zPQ - xPQ*yPQ*zPA) + v**2*xPQ*yQD*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xPB*yQD*zPA * F0
      end if
    
      ! Combination 45: s p_x p_z p_y
      if (o1 == "s" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(u*(xPB*yPQ*zQC + xPB*yQD*zPQ) - v*xPQ*yQD*zQC) + F2*(u**2*xPB*yPQ*zPQ + u*v*(-xPQ*yPQ*zQC - xPQ*yQD*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xPB*yQD*zQC * F0
      end if
    
      ! Combination 46: p_z p_x p_z s
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq*xPB + u*xPB*zPA*zPQ + v*(-xPB*zPQ*zQC - xPQ*zPA*zQC)) + F2*(v**2*xPQ*zPQ*zQC + v*(-pq*xPQ + u*(-xPB*zPQ**2 - xPQ*zPA*zPQ))) + F3*u*v**2*xPQ*zPQ**2 + xPB*zPA*zQC * F0
      end if
    
      ! Combination 47: p_z p_x s p_z
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq*xPB + u*xPB*zPA*zPQ + v*(-xPB*zPQ*zQD - xPQ*zPA*zQD)) + F2*(v**2*xPQ*zPQ*zQD + v*(-pq*xPQ + u*(-xPB*zPQ**2 - xPQ*zPA*zPQ))) + F3*u*v**2*xPQ*zPQ**2 + xPB*zPA*zQD * F0
      end if
    
      ! Combination 48: s p_x p_z p_z
      if (o1 == "s" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*xPB + u*(xPB*zPQ*zQC + xPB*zPQ*zQD) + v*(-xPQ*zQC*zQD - xPQ/(2*q))) + F2*(u**2*xPB*zPQ**2 + v*(tv*xPQ + u*(-xPQ*zPQ*zQC - xPQ*zPQ*zQD))) - F3*u**2*v*xPQ*zPQ**2 + (xPB*zQC*zQD + xPB/(2*q)) * F0
      end if
    
      ! Combination 49: p_y s p_x p_x
      if (o1 == "py" .and. o2 == "s" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*yPA + u*(xPQ*xQC*yPA + xPQ*xQD*yPA) + v*(-xQC*xQD*yPQ - yPQ/(2*q))) + F2*(u**2*xPQ**2*yPA + v*(tv*yPQ + u*(-xPQ*xQC*yPQ - xPQ*xQD*yPQ))) - F3*u**2*v*xPQ**2*yPQ + (xQC*xQD*yPA + yPA/(2*q)) * F0
      end if
    
      ! Combination 50: s p_y p_x p_x
      if (o1 == "s" .and. o2 == "py" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*yPB + u*(xPQ*xQC*yPB + xPQ*xQD*yPB) + v*(-xQC*xQD*yPQ - yPQ/(2*q))) + F2*(u**2*xPQ**2*yPB + v*(tv*yPQ + u*(-xPQ*xQC*yPQ - xPQ*xQD*yPQ))) - F3*u**2*v*xPQ**2*yPQ + (xQC*xQD*yPB + yPB/(2*q)) * F0
      end if
    
      ! Combination 51: p_z s p_x p_x
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*zPA + u*(xPQ*xQC*zPA + xPQ*xQD*zPA) + v*(-xQC*xQD*zPQ - zPQ/(2*q))) + F2*(u**2*xPQ**2*zPA + v*(tv*zPQ + u*(-xPQ*xQC*zPQ - xPQ*xQD*zPQ))) - F3*u**2*v*xPQ**2*zPQ + (xQC*xQD*zPA + zPA/(2*q)) * F0
      end if
    
      ! Combination 52: s p_z p_x p_x
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*zPB + u*(xPQ*xQC*zPB + xPQ*xQD*zPB) + v*(-xQC*xQD*zPQ - zPQ/(2*q))) + F2*(u**2*xPQ**2*zPB + v*(tv*zPQ + u*(-xPQ*xQC*zPQ - xPQ*xQD*zPQ))) - F3*u**2*v*xPQ**2*zPQ + (xQC*xQD*zPB + zPB/(2*q)) * F0
      end if
    
      ! Combination 53: p_y p_y p_x s
      if (o1 == "py" .and. o2 == "py" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(-tu*xQC + u*(xPQ*yPA*yPB + xPQ/(2*p)) + v*(-xQC*yPA*yPQ - xQC*yPB*yPQ)) + F2*(u*(-tu*xPQ + v*(-xPQ*yPA*yPQ - xPQ*yPB*yPQ)) + v**2*xQC*yPQ**2) + F3*u*v**2*xPQ*yPQ**2 + (xQC*yPA*yPB + xQC/(2*p)) * F0
      end if
    
      ! Combination 54: p_y s p_x p_y
      if (o1 == "py" .and. o2 == "s" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*xQC + u*(xPQ*yPA*yQD + xQC*yPA*yPQ) - v*xQC*yPQ*yQD) + F2*(u**2*xPQ*yPA*yPQ + u*(pq*xPQ + v*(-xPQ*yPQ*yQD - xQC*yPQ**2))) - F3*u**2*v*xPQ*yPQ**2 + xQC*yPA*yQD * F0
      end if
    
      ! Combination 55: p_y p_z p_x s
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*yPA*zPB + v*(-xQC*yPA*zPQ - xQC*yPQ*zPB)) + F2*(u*v*(-xPQ*yPA*zPQ - xPQ*yPQ*zPB) + v**2*xQC*yPQ*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xQC*yPA*zPB * F0
      end if
    
      ! Combination 56: p_y s p_x p_z
      if (o1 == "py" .and. o2 == "s" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(u*(xPQ*yPA*zQD + xQC*yPA*zPQ) - v*xQC*yPQ*zQD) + F2*(u**2*xPQ*yPA*zPQ + u*v*(-xPQ*yPQ*zQD - xQC*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQC*yPA*zQD * F0
      end if
    
      ! Combination 57: s p_y p_x p_y
      if (o1 == "s" .and. o2 == "py" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*xQC + u*(xPQ*yPB*yQD + xQC*yPB*yPQ) - v*xQC*yPQ*yQD) + F2*(u**2*xPQ*yPB*yPQ + u*(pq*xPQ + v*(-xPQ*yPQ*yQD - xQC*yPQ**2))) - F3*u**2*v*xPQ*yPQ**2 + xQC*yPB*yQD * F0
      end if
    
      ! Combination 58: p_z p_y p_x s
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(u*xPQ*yPB*zPA + v*(-xQC*yPB*zPQ - xQC*yPQ*zPA)) + F2*(u*v*(-xPQ*yPB*zPQ - xPQ*yPQ*zPA) + v**2*xQC*yPQ*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xQC*yPB*zPA * F0
      end if
    
      ! Combination 59: s p_y p_x p_z
      if (o1 == "s" .and. o2 == "py" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(u*(xPQ*yPB*zQD + xQC*yPB*zPQ) - v*xQC*yPQ*zQD) + F2*(u**2*xPQ*yPB*zPQ + u*v*(-xPQ*yPQ*zQD - xQC*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQC*yPB*zQD * F0
      end if
    
      ! Combination 60: p_z s p_x p_y
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(u*(xPQ*yQD*zPA + xQC*yPQ*zPA) - v*xQC*yQD*zPQ) + F2*(u**2*xPQ*yPQ*zPA + u*v*(-xPQ*yQD*zPQ - xQC*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQC*yQD*zPA * F0
      end if
    
      ! Combination 61: s p_z p_x p_y
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(u*(xPQ*yQD*zPB + xQC*yPQ*zPB) - v*xQC*yQD*zPQ) + F2*(u**2*xPQ*yPQ*zPB + u*v*(-xPQ*yQD*zPQ - xQC*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQC*yQD*zPB * F0
      end if
    
      ! Combination 62: p_z p_z p_x s
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "s") then
        value = F1*(-tu*xQC + u*(xPQ*zPA*zPB + xPQ/(2*p)) + v*(-xQC*zPA*zPQ - xQC*zPB*zPQ)) + F2*(u*(-tu*xPQ + v*(-xPQ*zPA*zPQ - xPQ*zPB*zPQ)) + v**2*xQC*zPQ**2) + F3*u*v**2*xPQ*zPQ**2 + (xQC*zPA*zPB + xQC/(2*p)) * F0
      end if
    
      ! Combination 63: p_z s p_x p_z
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*xQC + u*(xPQ*zPA*zQD + xQC*zPA*zPQ) - v*xQC*zPQ*zQD) + F2*(u**2*xPQ*zPA*zPQ + u*(pq*xPQ + v*(-xPQ*zPQ*zQD - xQC*zPQ**2))) - F3*u**2*v*xPQ*zPQ**2 + xQC*zPA*zQD * F0
      end if
    
      ! Combination 64: s p_z p_x p_z
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*xQC + u*(xPQ*zPB*zQD + xQC*zPB*zPQ) - v*xQC*zPQ*zQD) + F2*(u**2*xPQ*zPB*zPQ + u*(pq*xPQ + v*(-xPQ*zPQ*zQD - xQC*zPQ**2))) - F3*u**2*v*xPQ*zPQ**2 + xQC*zPB*zQD * F0
      end if
    
      ! Combination 65: p_y p_y s p_x
      if (o1 == "py" .and. o2 == "py" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(-tu*xQD + u*(xPQ*yPA*yPB + xPQ/(2*p)) + v*(-xQD*yPA*yPQ - xQD*yPB*yPQ)) + F2*(u*(-tu*xPQ + v*(-xPQ*yPA*yPQ - xPQ*yPB*yPQ)) + v**2*xQD*yPQ**2) + F3*u*v**2*xPQ*yPQ**2 + (xQD*yPA*yPB + xQD/(2*p)) * F0
      end if
    
      ! Combination 66: p_y s p_y p_x
      if (o1 == "py" .and. o2 == "s" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*xQD + u*(xPQ*yPA*yQC + xQD*yPA*yPQ) - v*xQD*yPQ*yQC) + F2*(u**2*xPQ*yPA*yPQ + u*(pq*xPQ + v*(-xPQ*yPQ*yQC - xQD*yPQ**2))) - F3*u**2*v*xPQ*yPQ**2 + xQD*yPA*yQC * F0
      end if
    
      ! Combination 67: p_y p_z s p_x
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*yPA*zPB + v*(-xQD*yPA*zPQ - xQD*yPQ*zPB)) + F2*(u*v*(-xPQ*yPA*zPQ - xPQ*yPQ*zPB) + v**2*xQD*yPQ*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xQD*yPA*zPB * F0
      end if
    
      ! Combination 68: p_y s p_z p_x
      if (o1 == "py" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(u*(xPQ*yPA*zQC + xQD*yPA*zPQ) - v*xQD*yPQ*zQC) + F2*(u**2*xPQ*yPA*zPQ + u*v*(-xPQ*yPQ*zQC - xQD*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQD*yPA*zQC * F0
      end if
    
      ! Combination 69: s p_y p_y p_x
      if (o1 == "s" .and. o2 == "py" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*xQD + u*(xPQ*yPB*yQC + xQD*yPB*yPQ) - v*xQD*yPQ*yQC) + F2*(u**2*xPQ*yPB*yPQ + u*(pq*xPQ + v*(-xPQ*yPQ*yQC - xQD*yPQ**2))) - F3*u**2*v*xPQ*yPQ**2 + xQD*yPB*yQC * F0
      end if
    
      ! Combination 70: p_z p_y s p_x
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(u*xPQ*yPB*zPA + v*(-xQD*yPB*zPQ - xQD*yPQ*zPA)) + F2*(u*v*(-xPQ*yPB*zPQ - xPQ*yPQ*zPA) + v**2*xQD*yPQ*zPQ) + F3*u*v**2*xPQ*yPQ*zPQ + xQD*yPB*zPA * F0
      end if
    
      ! Combination 71: s p_y p_z p_x
      if (o1 == "s" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(u*(xPQ*yPB*zQC + xQD*yPB*zPQ) - v*xQD*yPQ*zQC) + F2*(u**2*xPQ*yPB*zPQ + u*v*(-xPQ*yPQ*zQC - xQD*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQD*yPB*zQC * F0
      end if
    
      ! Combination 72: p_z s p_y p_x
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(u*(xPQ*yQC*zPA + xQD*yPQ*zPA) - v*xQD*yQC*zPQ) + F2*(u**2*xPQ*yPQ*zPA + u*v*(-xPQ*yQC*zPQ - xQD*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQD*yQC*zPA * F0
      end if
    
      ! Combination 73: s p_z p_y p_x
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(u*(xPQ*yQC*zPB + xQD*yPQ*zPB) - v*xQD*yQC*zPQ) + F2*(u**2*xPQ*yPQ*zPB + u*v*(-xPQ*yQC*zPQ - xQD*yPQ*zPQ)) - F3*u**2*v*xPQ*yPQ*zPQ + xQD*yQC*zPB * F0
      end if
    
      ! Combination 74: p_z p_z s p_x
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "px") then
        value = F1*(-tu*xQD + u*(xPQ*zPA*zPB + xPQ/(2*p)) + v*(-xQD*zPA*zPQ - xQD*zPB*zPQ)) + F2*(u*(-tu*xPQ + v*(-xPQ*zPA*zPQ - xPQ*zPB*zPQ)) + v**2*xQD*zPQ**2) + F3*u*v**2*xPQ*zPQ**2 + (xQD*zPA*zPB + xQD/(2*p)) * F0
      end if
    
      ! Combination 75: p_z s p_z p_x
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*xQD + u*(xPQ*zPA*zQC + xQD*zPA*zPQ) - v*xQD*zPQ*zQC) + F2*(u**2*xPQ*zPA*zPQ + u*(pq*xPQ + v*(-xPQ*zPQ*zQC - xQD*zPQ**2))) - F3*u**2*v*xPQ*zPQ**2 + xQD*zPA*zQC * F0
      end if
    
      ! Combination 76: s p_z p_z p_x
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*xQD + u*(xPQ*zPB*zQC + xQD*zPB*zPQ) - v*xQD*zPQ*zQC) + F2*(u**2*xPQ*zPB*zPQ + u*(pq*xPQ + v*(-xPQ*zPQ*zQC - xQD*zPQ**2))) - F3*u**2*v*xPQ*zPQ**2 + xQD*zPB*zQC * F0
      end if
    
      ! Combination 77: p_y p_y p_y s
      if (o1 == "py" .and. o2 == "py" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq*(yPA + yPB) - tu*yQC + u*(yPA*yPB*yPQ + yPQ/(2*p)) + v*(-yPA*yPQ*yQC - yPB*yPQ*yQC)) + F2*(-2*pq*v*yPQ + u*(-tu*yPQ + v*(-yPA*yPQ**2 - yPB*yPQ**2)) + v**2*yPQ**2*yQC) + F3*u*v**2*yPQ**3 + (yPA*yPB*yQC + yQC/(2*p)) * F0
      end if
    
      ! Combination 78: p_y p_y s p_y
      if (o1 == "py" .and. o2 == "py" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq*(yPA + yPB) - tu*yQD + u*(yPA*yPB*yPQ + yPQ/(2*p)) + v*(-yPA*yPQ*yQD - yPB*yPQ*yQD)) + F2*(-2*pq*v*yPQ + u*(-tu*yPQ + v*(-yPA*yPQ**2 - yPB*yPQ**2)) + v**2*yPQ**2*yQD) + F3*u*v**2*yPQ**3 + (yPA*yPB*yQD + yQD/(2*p)) * F0
      end if
    
      ! Combination 79: p_y p_y p_z s
      if (o1 == "py" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(-tu*zQC + u*(yPA*yPB*zPQ + zPQ/(2*p)) + v*(-yPA*yPQ*zQC - yPB*yPQ*zQC)) + F2*(u*(-tu*zPQ + v*(-yPA*yPQ*zPQ - yPB*yPQ*zPQ)) + v**2*yPQ**2*zQC) + F3*u*v**2*yPQ**2*zPQ + (yPA*yPB*zQC + zQC/(2*p)) * F0
      end if
    
      ! Combination 80: p_y p_y s p_z
      if (o1 == "py" .and. o2 == "py" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(-tu*zQD + u*(yPA*yPB*zPQ + zPQ/(2*p)) + v*(-yPA*yPQ*zQD - yPB*yPQ*zQD)) + F2*(u*(-tu*zPQ + v*(-yPA*yPQ*zPQ - yPB*yPQ*zPQ)) + v**2*yPQ**2*zQD) + F3*u*v**2*yPQ**2*zPQ + (yPA*yPB*zQD + zQD/(2*p)) * F0
      end if
    
      ! Combination 81: p_y s p_y p_y
      if (o1 == "py" .and. o2 == "s" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(yQC + yQD) - tv*yPA + u*(yPA*yPQ*yQC + yPA*yPQ*yQD) + v*(-yPQ*yQC*yQD - yPQ/(2*q))) + F2*(tv*v*yPQ + u**2*yPA*yPQ**2 + u*(2*pq*yPQ + v*(-yPQ**2*yQC - yPQ**2*yQD))) - F3*u**2*v*yPQ**3 + (yPA*yQC*yQD + yPA/(2*q)) * F0
      end if
    
      ! Combination 82: p_y p_z p_y s
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq*zPB + u*yPA*yPQ*zPB + v*(-yPA*yQC*zPQ - yPQ*yQC*zPB)) + F2*(v**2*yPQ*yQC*zPQ + v*(-pq*zPQ + u*(-yPA*yPQ*zPQ - yPQ**2*zPB))) + F3*u*v**2*yPQ**2*zPQ + yPA*yQC*zPB * F0
      end if
    
      ! Combination 83: p_y s p_y p_z
      if (o1 == "py" .and. o2 == "s" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*zQD + u*(yPA*yPQ*zQD + yPA*yQC*zPQ) - v*yPQ*yQC*zQD) + F2*(u**2*yPA*yPQ*zPQ + u*(pq*zPQ + v*(-yPQ**2*zQD - yPQ*yQC*zPQ))) - F3*u**2*v*yPQ**2*zPQ + yPA*yQC*zQD * F0
      end if
    
      ! Combination 84: p_y p_z s p_y
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq*zPB + u*yPA*yPQ*zPB + v*(-yPA*yQD*zPQ - yPQ*yQD*zPB)) + F2*(v**2*yPQ*yQD*zPQ + v*(-pq*zPQ + u*(-yPA*yPQ*zPQ - yPQ**2*zPB))) + F3*u*v**2*yPQ**2*zPQ + yPA*yQD*zPB * F0
      end if
    
      ! Combination 85: p_y s p_z p_y
      if (o1 == "py" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*zQC + u*(yPA*yPQ*zQC + yPA*yQD*zPQ) - v*yPQ*yQD*zQC) + F2*(u**2*yPA*yPQ*zPQ + u*(pq*zPQ + v*(-yPQ**2*zQC - yPQ*yQD*zPQ))) - F3*u**2*v*yPQ**2*zPQ + yPA*yQD*zQC * F0
      end if
    
      ! Combination 86: p_y p_z p_z s
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq*yPA + u*yPA*zPB*zPQ + v*(-yPA*zPQ*zQC - yPQ*zPB*zQC)) + F2*(v**2*yPQ*zPQ*zQC + v*(-pq*yPQ + u*(-yPA*zPQ**2 - yPQ*zPB*zPQ))) + F3*u*v**2*yPQ*zPQ**2 + yPA*zPB*zQC * F0
      end if
    
      ! Combination 87: p_y p_z s p_z
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq*yPA + u*yPA*zPB*zPQ + v*(-yPA*zPQ*zQD - yPQ*zPB*zQD)) + F2*(v**2*yPQ*zPQ*zQD + v*(-pq*yPQ + u*(-yPA*zPQ**2 - yPQ*zPB*zPQ))) + F3*u*v**2*yPQ*zPQ**2 + yPA*zPB*zQD * F0
      end if
    
      ! Combination 88: p_y s p_z p_z
      if (o1 == "py" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*yPA + u*(yPA*zPQ*zQC + yPA*zPQ*zQD) + v*(-yPQ*zQC*zQD - yPQ/(2*q))) + F2*(u**2*yPA*zPQ**2 + v*(tv*yPQ + u*(-yPQ*zPQ*zQC - yPQ*zPQ*zQD))) - F3*u**2*v*yPQ*zPQ**2 + (yPA*zQC*zQD + yPA/(2*q)) * F0
      end if
    
      ! Combination 89: s p_y p_y p_y
      if (o1 == "s" .and. o2 == "py" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(yQC + yQD) - tv*yPB + u*(yPB*yPQ*yQC + yPB*yPQ*yQD) + v*(-yPQ*yQC*yQD - yPQ/(2*q))) + F2*(tv*v*yPQ + u**2*yPB*yPQ**2 + u*(2*pq*yPQ + v*(-yPQ**2*yQC - yPQ**2*yQD))) - F3*u**2*v*yPQ**3 + (yPB*yQC*yQD + yPB/(2*q)) * F0
      end if
    
      ! Combination 90: p_z p_y p_y s
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(pq*zPA + u*yPB*yPQ*zPA + v*(-yPB*yQC*zPQ - yPQ*yQC*zPA)) + F2*(v**2*yPQ*yQC*zPQ + v*(-pq*zPQ + u*(-yPB*yPQ*zPQ - yPQ**2*zPA))) + F3*u*v**2*yPQ**2*zPQ + yPB*yQC*zPA * F0
      end if
    
      ! Combination 91: s p_y p_y p_z
      if (o1 == "s" .and. o2 == "py" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*zQD + u*(yPB*yPQ*zQD + yPB*yQC*zPQ) - v*yPQ*yQC*zQD) + F2*(u**2*yPB*yPQ*zPQ + u*(pq*zPQ + v*(-yPQ**2*zQD - yPQ*yQC*zPQ))) - F3*u**2*v*yPQ**2*zPQ + yPB*yQC*zQD * F0
      end if
    
      ! Combination 92: p_z p_y s p_y
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(pq*zPA + u*yPB*yPQ*zPA + v*(-yPB*yQD*zPQ - yPQ*yQD*zPA)) + F2*(v**2*yPQ*yQD*zPQ + v*(-pq*zPQ + u*(-yPB*yPQ*zPQ - yPQ**2*zPA))) + F3*u*v**2*yPQ**2*zPQ + yPB*yQD*zPA * F0
      end if
    
      ! Combination 93: s p_y p_z p_y
      if (o1 == "s" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*zQC + u*(yPB*yPQ*zQC + yPB*yQD*zPQ) - v*yPQ*yQD*zQC) + F2*(u**2*yPB*yPQ*zPQ + u*(pq*zPQ + v*(-yPQ**2*zQC - yPQ*yQD*zPQ))) - F3*u**2*v*yPQ**2*zPQ + yPB*yQD*zQC * F0
      end if
    
      ! Combination 94: p_z p_y p_z s
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq*yPB + u*yPB*zPA*zPQ + v*(-yPB*zPQ*zQC - yPQ*zPA*zQC)) + F2*(v**2*yPQ*zPQ*zQC + v*(-pq*yPQ + u*(-yPB*zPQ**2 - yPQ*zPA*zPQ))) + F3*u*v**2*yPQ*zPQ**2 + yPB*zPA*zQC * F0
      end if
    
      ! Combination 95: p_z p_y s p_z
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq*yPB + u*yPB*zPA*zPQ + v*(-yPB*zPQ*zQD - yPQ*zPA*zQD)) + F2*(v**2*yPQ*zPQ*zQD + v*(-pq*yPQ + u*(-yPB*zPQ**2 - yPQ*zPA*zPQ))) + F3*u*v**2*yPQ*zPQ**2 + yPB*zPA*zQD * F0
      end if
    
      ! Combination 96: s p_y p_z p_z
      if (o1 == "s" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*yPB + u*(yPB*zPQ*zQC + yPB*zPQ*zQD) + v*(-yPQ*zQC*zQD - yPQ/(2*q))) + F2*(u**2*yPB*zPQ**2 + v*(tv*yPQ + u*(-yPQ*zPQ*zQC - yPQ*zPQ*zQD))) - F3*u**2*v*yPQ*zPQ**2 + (yPB*zQC*zQD + yPB/(2*q)) * F0
      end if
    
      ! Combination 97: p_z s p_y p_y
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*zPA + u*(yPQ*yQC*zPA + yPQ*yQD*zPA) + v*(-yQC*yQD*zPQ - zPQ/(2*q))) + F2*(u**2*yPQ**2*zPA + v*(tv*zPQ + u*(-yPQ*yQC*zPQ - yPQ*yQD*zPQ))) - F3*u**2*v*yPQ**2*zPQ + (yQC*yQD*zPA + zPA/(2*q)) * F0
      end if
    
      ! Combination 98: s p_z p_y p_y
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*zPB + u*(yPQ*yQC*zPB + yPQ*yQD*zPB) + v*(-yQC*yQD*zPQ - zPQ/(2*q))) + F2*(u**2*yPQ**2*zPB + v*(tv*zPQ + u*(-yPQ*yQC*zPQ - yPQ*yQD*zPQ))) - F3*u**2*v*yPQ**2*zPQ + (yQC*yQD*zPB + zPB/(2*q)) * F0
      end if
    
      ! Combination 99: p_z p_z p_y s
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "s") then
        value = F1*(-tu*yQC + u*(yPQ*zPA*zPB + yPQ/(2*p)) + v*(-yQC*zPA*zPQ - yQC*zPB*zPQ)) + F2*(u*(-tu*yPQ + v*(-yPQ*zPA*zPQ - yPQ*zPB*zPQ)) + v**2*yQC*zPQ**2) + F3*u*v**2*yPQ*zPQ**2 + (yQC*zPA*zPB + yQC/(2*p)) * F0
      end if
    
      ! Combination 100: p_z s p_y p_z
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*yQC + u*(yPQ*zPA*zQD + yQC*zPA*zPQ) - v*yQC*zPQ*zQD) + F2*(u**2*yPQ*zPA*zPQ + u*(pq*yPQ + v*(-yPQ*zPQ*zQD - yQC*zPQ**2))) - F3*u**2*v*yPQ*zPQ**2 + yQC*zPA*zQD * F0
      end if
    
      ! Combination 101: s p_z p_y p_z
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*yQC + u*(yPQ*zPB*zQD + yQC*zPB*zPQ) - v*yQC*zPQ*zQD) + F2*(u**2*yPQ*zPB*zPQ + u*(pq*yPQ + v*(-yPQ*zPQ*zQD - yQC*zPQ**2))) - F3*u**2*v*yPQ*zPQ**2 + yQC*zPB*zQD * F0
      end if
    
      ! Combination 102: p_z p_z s p_y
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "py") then
        value = F1*(-tu*yQD + u*(yPQ*zPA*zPB + yPQ/(2*p)) + v*(-yQD*zPA*zPQ - yQD*zPB*zPQ)) + F2*(u*(-tu*yPQ + v*(-yPQ*zPA*zPQ - yPQ*zPB*zPQ)) + v**2*yQD*zPQ**2) + F3*u*v**2*yPQ*zPQ**2 + (yQD*zPA*zPB + yQD/(2*p)) * F0
      end if
    
      ! Combination 103: p_z s p_z p_y
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*yQD + u*(yPQ*zPA*zQC + yQD*zPA*zPQ) - v*yQD*zPQ*zQC) + F2*(u**2*yPQ*zPA*zPQ + u*(pq*yPQ + v*(-yPQ*zPQ*zQC - yQD*zPQ**2))) - F3*u**2*v*yPQ*zPQ**2 + yQD*zPA*zQC * F0
      end if
    
      ! Combination 104: s p_z p_z p_y
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*yQD + u*(yPQ*zPB*zQC + yQD*zPB*zPQ) - v*yQD*zPQ*zQC) + F2*(u**2*yPQ*zPB*zPQ + u*(pq*yPQ + v*(-yPQ*zPQ*zQC - yQD*zPQ**2))) - F3*u**2*v*yPQ*zPQ**2 + yQD*zPB*zQC * F0
      end if
    
      ! Combination 105: p_z p_z p_z s
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "s") then
        value = F1*(pq*(zPA + zPB) - tu*zQC + u*(zPA*zPB*zPQ + zPQ/(2*p)) + v*(-zPA*zPQ*zQC - zPB*zPQ*zQC)) + F2*(-2*pq*v*zPQ + u*(-tu*zPQ + v*(-zPA*zPQ**2 - zPB*zPQ**2)) + v**2*zPQ**2*zQC) + F3*u*v**2*zPQ**3 + (zPA*zPB*zQC + zQC/(2*p)) * F0
      end if
    
      ! Combination 106: p_z p_z s p_z
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "s" .and. o4 == "pz") then
        value = F1*(pq*(zPA + zPB) - tu*zQD + u*(zPA*zPB*zPQ + zPQ/(2*p)) + v*(-zPA*zPQ*zQD - zPB*zPQ*zQD)) + F2*(-2*pq*v*zPQ + u*(-tu*zPQ + v*(-zPA*zPQ**2 - zPB*zPQ**2)) + v**2*zPQ**2*zQD) + F3*u*v**2*zPQ**3 + (zPA*zPB*zQD + zQD/(2*p)) * F0
      end if
    
      ! Combination 107: p_z s p_z p_z
      if (o1 == "pz" .and. o2 == "s" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(zQC + zQD) - tv*zPA + u*(zPA*zPQ*zQC + zPA*zPQ*zQD) + v*(-zPQ*zQC*zQD - zPQ/(2*q))) + F2*(tv*v*zPQ + u**2*zPA*zPQ**2 + u*(2*pq*zPQ + v*(-zPQ**2*zQC - zPQ**2*zQD))) - F3*u**2*v*zPQ**3 + (zPA*zQC*zQD + zPA/(2*q)) * F0
      end if
    
      ! Combination 108: s p_z p_z p_z
      if (o1 == "s" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(zQC + zQD) - tv*zPB + u*(zPB*zPQ*zQC + zPB*zPQ*zQD) + v*(-zPQ*zQC*zQD - zPQ/(2*q))) + F2*(tv*v*zPQ + u**2*zPB*zPQ**2 + u*(2*pq*zPQ + v*(-zPQ**2*zQC - zPQ**2*zQD))) - F3*u**2*v*zPQ**3 + (zPB*zQC*zQD + zPB/(2*q)) * F0
      end if

end subroutine

subroutine ERI_pppp(o1,o2,o3,o4&
                  &,F0,F1,F2,F3,F4&
                  &,u,v,tu,tv,p,q,pq&
                  &,xPA,yPA,zPA,xPB,yPB,zPB&
                  &,xQC,yQC,zQC,xQD,yQD,zQD&
                  &,xPQ,yPQ,zPQ,value)


      implicit none 


      character(LEN=2),intent(in)  :: o1 , o2 , o3 , o4
      double precision,intent(in)  :: F0 , F1 , F2 , F3 , F4 
      

      double precision,intent(in)  :: xPA,yPA,zPA,xPB,yPB,zPB
      double precision,intent(in)  :: xQC,yQC,zQC,xQD,yQD,zQD

      double precision,intent(in)  :: xPQ , yPQ , zPQ

      double precision,intent(in)  :: u , v , tu , tv , p ,q  , pq 
      
      double precision,intent(out) :: value

      ! Combination 1: p_x p_x p_x p_x
      if (o1 == "px" .and. o2 == "px" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xPA*xQC + xPA*xQD + xPB*xQC + xPB*xQD) - tu*xQC*xQD - tv*xPA*xPB + u*(xPA*xPB*xPQ*xQC + xPA*xPB*xPQ*xQD + xPQ*xQC/(2*p) + xPQ*xQD/(2*p)) + v*(-xPA*xPQ*xQC*xQD - xPB*xPQ*xQC*xQD - xPA*xPQ/(2*q) - xPB*xPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(2*pq**2 + tu*tv + u**2*(xPA*xPB*xPQ**2 + xPQ**2/(2*p)) + u*(pq*(2*xPA*xPQ + 2*xPB*xPQ) - tu*xPQ*xQC - tu*xPQ*xQD + v*(-xPA*xPQ**2*xQC - xPA*xPQ**2*xQD - xPB*xPQ**2*xQC - xPB*xPQ**2*xQD)) + v**2*(xPQ**2*xQC*xQD + xPQ**2/(2*q)) + v*(pq*(-2*xPQ*xQC - 2*xPQ*xQD) + tv*xPA*xPQ + tv*xPB*xPQ)) + F3*(-tv*v**2*xPQ**2 + u**2*(-tu*xPQ**2 + v*(-xPA*xPQ**3 - xPB*xPQ**3)) + u*(-4*pq*v*xPQ**2 + v**2*(xPQ**3*xQC + xPQ**3*xQD))) + F4*u**2*v**2*xPQ**4 + (xPA*xPB*xQC*xQD + xPA*xPB/(2*q) + xQC*xQD/(2*p) + 1/(4*p*q)) * F0 
      end if
        
      ! Combination 2: p_x p_x p_x p_y
      if (o1 == "px" .and. o2 == "px" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*(xPA*yQD + xPB*yQD) - tu*xQC*yQD + u*(xPA*xPB*xPQ*yQD + xPA*xPB*xQC*yPQ + xPQ*yQD/(2*p) + xQC*yPQ/(2*p)) + v*(-xPA*xPQ*xQC*yQD - xPB*xPQ*xQC*yQD)) + F2*(-2*pq*v*xPQ*yQD + u**2*(xPA*xPB*xPQ*yPQ + xPQ*yPQ/(2*p)) + u*(pq*(xPA*yPQ + xPB*yPQ) - tu*xPQ*yQD - tu*xQC*yPQ + v*(-xPA*xPQ**2*yQD - xPA*xPQ*xQC*yPQ - xPB*xPQ**2*yQD - xPB*xPQ*xQC*yPQ)) + v**2*xPQ**2*xQC*yQD) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPA*xPQ**2*yPQ - xPB*xPQ**2*yPQ)) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ**3*yQD + xPQ**2*xQC*yPQ))) + F4*u**2*v**2*xPQ**3*yPQ + (xPA*xPB*xQC*yQD + xQC*yQD/(2*p)) * F0 
      end if
    
      ! Combination 3: p_x p_x p_x p_z
      if (o1 == "px" .and. o2 == "px" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*(xPA*zQD + xPB*zQD) - tu*xQC*zQD + u*(xPA*xPB*xPQ*zQD + xPA*xPB*xQC*zPQ + xPQ*zQD/(2*p) + xQC*zPQ/(2*p)) + v*(-xPA*xPQ*xQC*zQD - xPB*xPQ*xQC*zQD)) + F2*(-2*pq*v*xPQ*zQD + u**2*(xPA*xPB*xPQ*zPQ + xPQ*zPQ/(2*p)) + u*(pq*(xPA*zPQ + xPB*zPQ) - tu*xPQ*zQD - tu*xQC*zPQ + v*(-xPA*xPQ**2*zQD - xPA*xPQ*xQC*zPQ - xPB*xPQ**2*zQD - xPB*xPQ*xQC*zPQ)) + v**2*xPQ**2*xQC*zQD) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPA*xPQ**2*zPQ - xPB*xPQ**2*zPQ)) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ**3*zQD + xPQ**2*xQC*zPQ))) + F4*u**2*v**2*xPQ**3*zPQ + (xPA*xPB*xQC*zQD + xQC*zQD/(2*p)) * F0 
      end if
    
      ! Combination 4: p_x p_x p_y p_x
      if (o1 == "px" .and. o2 == "px" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*(xPA*yQC + xPB*yQC) - tu*xQD*yQC + u*(xPA*xPB*xPQ*yQC + xPA*xPB*xQD*yPQ + xPQ*yQC/(2*p) + xQD*yPQ/(2*p)) + v*(-xPA*xPQ*xQD*yQC - xPB*xPQ*xQD*yQC)) + F2*(-2*pq*v*xPQ*yQC + u**2*(xPA*xPB*xPQ*yPQ + xPQ*yPQ/(2*p)) + u*(pq*(xPA*yPQ + xPB*yPQ) - tu*xPQ*yQC - tu*xQD*yPQ + v*(-xPA*xPQ**2*yQC - xPA*xPQ*xQD*yPQ - xPB*xPQ**2*yQC - xPB*xPQ*xQD*yPQ)) + v**2*xPQ**2*xQD*yQC) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPA*xPQ**2*yPQ - xPB*xPQ**2*yPQ)) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ**3*yQC + xPQ**2*xQD*yPQ))) + F4*u**2*v**2*xPQ**3*yPQ + (xPA*xPB*xQD*yQC + xQD*yQC/(2*p)) * F0 
      end if
    
      ! Combination 5: p_x p_x p_z p_x
      if (o1 == "px" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*(xPA*zQC + xPB*zQC) - tu*xQD*zQC + u*(xPA*xPB*xPQ*zQC + xPA*xPB*xQD*zPQ + xPQ*zQC/(2*p) + xQD*zPQ/(2*p)) + v*(-xPA*xPQ*xQD*zQC - xPB*xPQ*xQD*zQC)) + F2*(-2*pq*v*xPQ*zQC + u**2*(xPA*xPB*xPQ*zPQ + xPQ*zPQ/(2*p)) + u*(pq*(xPA*zPQ + xPB*zPQ) - tu*xPQ*zQC - tu*xQD*zPQ + v*(-xPA*xPQ**2*zQC - xPA*xPQ*xQD*zPQ - xPB*xPQ**2*zQC - xPB*xPQ*xQD*zPQ)) + v**2*xPQ**2*xQD*zQC) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPA*xPQ**2*zPQ - xPB*xPQ**2*zPQ)) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ**3*zQC + xPQ**2*xQD*zPQ))) + F4*u**2*v**2*xPQ**3*zPQ + (xPA*xPB*xQD*zQC + xQD*zQC/(2*p)) * F0 
      end if
    
      ! Combination 6: p_x p_x p_y p_y
      if (o1 == "px" .and. o2 == "px" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tu*yQC*yQD - tv*xPA*xPB + u*(xPA*xPB*yPQ*yQC + xPA*xPB*yPQ*yQD + yPQ*yQC/(2*p) + yPQ*yQD/(2*p)) + v*(-xPA*xPQ*yQC*yQD - xPB*xPQ*yQC*yQD - xPA*xPQ/(2*q) - xPB*xPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(xPA*xPB*yPQ**2 + yPQ**2/(2*p)) + u*(-tu*yPQ*yQC - tu*yPQ*yQD + v*(-xPA*xPQ*yPQ*yQC - xPA*xPQ*yPQ*yQD - xPB*xPQ*yPQ*yQC - xPB*xPQ*yPQ*yQD)) + v**2*(xPQ**2*yQC*yQD + xPQ**2/(2*q)) + v*(tv*xPA*xPQ + tv*xPB*xPQ)) + F3*(u**2*(-tu*yPQ**2 + v*(-xPA*xPQ*yPQ**2 - xPB*xPQ*yPQ**2)) + v**2*(-tv*xPQ**2 + u*(xPQ**2*yPQ*yQC + xPQ**2*yPQ*yQD))) + F4*u**2*v**2*xPQ**2*yPQ**2 + (xPA*xPB*yQC*yQD + xPA*xPB/(2*q) + yQC*yQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 7: p_x p_x p_y p_z
      if (o1 == "px" .and. o2 == "px" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(-tu*yQC*zQD + u*(xPA*xPB*yPQ*zQD + xPA*xPB*yQC*zPQ + yPQ*zQD/(2*p) + yQC*zPQ/(2*p)) + v*(-xPA*xPQ*yQC*zQD - xPB*xPQ*yQC*zQD)) + F2*(u**2*(xPA*xPB*yPQ*zPQ + yPQ*zPQ/(2*p)) + u*(-tu*yPQ*zQD - tu*yQC*zPQ + v*(-xPA*xPQ*yPQ*zQD - xPA*xPQ*yQC*zPQ - xPB*xPQ*yPQ*zQD - xPB*xPQ*yQC*zPQ)) + v**2*xPQ**2*yQC*zQD) + F3*(u**2*(-tu*yPQ*zPQ + v*(-xPA*xPQ*yPQ*zPQ - xPB*xPQ*yPQ*zPQ)) + u*v**2*(xPQ**2*yPQ*zQD + xPQ**2*yQC*zPQ)) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + (xPA*xPB*yQC*zQD + yQC*zQD/(2*p)) * F0 
      end if
    
      ! Combination 8: p_x p_x p_z p_y
      if (o1 == "px" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(-tu*yQD*zQC + u*(xPA*xPB*yPQ*zQC + xPA*xPB*yQD*zPQ + yPQ*zQC/(2*p) + yQD*zPQ/(2*p)) + v*(-xPA*xPQ*yQD*zQC - xPB*xPQ*yQD*zQC)) + F2*(u**2*(xPA*xPB*yPQ*zPQ + yPQ*zPQ/(2*p)) + u*(-tu*yPQ*zQC - tu*yQD*zPQ + v*(-xPA*xPQ*yPQ*zQC - xPA*xPQ*yQD*zPQ - xPB*xPQ*yPQ*zQC - xPB*xPQ*yQD*zPQ)) + v**2*xPQ**2*yQD*zQC) + F3*(u**2*(-tu*yPQ*zPQ + v*(-xPA*xPQ*yPQ*zPQ - xPB*xPQ*yPQ*zPQ)) + u*v**2*(xPQ**2*yPQ*zQC + xPQ**2*yQD*zPQ)) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + (xPA*xPB*yQD*zQC + yQD*zQC/(2*p)) * F0 
      end if
    
      ! Combination 9: p_x p_x p_z p_z
      if (o1 == "px" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tu*zQC*zQD - tv*xPA*xPB + u*(xPA*xPB*zPQ*zQC + xPA*xPB*zPQ*zQD + zPQ*zQC/(2*p) + zPQ*zQD/(2*p)) + v*(-xPA*xPQ*zQC*zQD - xPB*xPQ*zQC*zQD - xPA*xPQ/(2*q) - xPB*xPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(xPA*xPB*zPQ**2 + zPQ**2/(2*p)) + u*(-tu*zPQ*zQC - tu*zPQ*zQD + v*(-xPA*xPQ*zPQ*zQC - xPA*xPQ*zPQ*zQD - xPB*xPQ*zPQ*zQC - xPB*xPQ*zPQ*zQD)) + v**2*(xPQ**2*zQC*zQD + xPQ**2/(2*q)) + v*(tv*xPA*xPQ + tv*xPB*xPQ)) + F3*(u**2*(-tu*zPQ**2 + v*(-xPA*xPQ*zPQ**2 - xPB*xPQ*zPQ**2)) + v**2*(-tv*xPQ**2 + u*(xPQ**2*zPQ*zQC + xPQ**2*zPQ*zQD))) + F4*u**2*v**2*xPQ**2*zPQ**2 + (xPA*xPB*zQC*zQD + xPA*xPB/(2*q) + zQC*zQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 10: p_x p_y p_x p_x
      if (o1 == "px" .and. o2 == "py" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC*yPB + xQD*yPB) - tv*xPA*yPB + u*(xPA*xPQ*xQC*yPB + xPA*xPQ*xQD*yPB) + v*(-xPA*xQC*xQD*yPQ - xPQ*xQC*xQD*yPB - xPA*yPQ/(2*q) - xPQ*yPB/(2*q))) + F2*(u**2*xPA*xPQ**2*yPB + u*(2*pq*xPQ*yPB + v*(-xPA*xPQ*xQC*yPQ - xPA*xPQ*xQD*yPQ - xPQ**2*xQC*yPB - xPQ**2*xQD*yPB)) + v**2*(xPQ*xQC*xQD*yPQ + xPQ*yPQ/(2*q)) + v*(pq*(-xQC*yPQ - xQD*yPQ) + tv*xPA*yPQ + tv*xPQ*yPB)) + F3*(-tv*v**2*xPQ*yPQ + u**2*v*(-xPA*xPQ**2*yPQ - xPQ**3*yPB) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ**2*xQC*yPQ + xPQ**2*xQD*yPQ))) + F4*u**2*v**2*xPQ**3*yPQ + (xPA*xQC*xQD*yPB + xPA*yPB/(2*q)) * F0 
      end if
    
      ! Combination 11: p_x p_z p_x p_x
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC*zPB + xQD*zPB) - tv*xPA*zPB + u*(xPA*xPQ*xQC*zPB + xPA*xPQ*xQD*zPB) + v*(-xPA*xQC*xQD*zPQ - xPQ*xQC*xQD*zPB - xPA*zPQ/(2*q) - xPQ*zPB/(2*q))) + F2*(u**2*xPA*xPQ**2*zPB + u*(2*pq*xPQ*zPB + v*(-xPA*xPQ*xQC*zPQ - xPA*xPQ*xQD*zPQ - xPQ**2*xQC*zPB - xPQ**2*xQD*zPB)) + v**2*(xPQ*xQC*xQD*zPQ + xPQ*zPQ/(2*q)) + v*(pq*(-xQC*zPQ - xQD*zPQ) + tv*xPA*zPQ + tv*xPQ*zPB)) + F3*(-tv*v**2*xPQ*zPQ + u**2*v*(-xPA*xPQ**2*zPQ - xPQ**3*zPB) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ**2*xQC*zPQ + xPQ**2*xQD*zPQ))) + F4*u**2*v**2*xPQ**3*zPQ + (xPA*xQC*xQD*zPB + xPA*zPB/(2*q)) * F0 
      end if
    
      ! Combination 12: p_x p_y p_x p_y
      if (o1 == "px" .and. o2 == "py" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*(xPA*xQC + yPB*yQD) + u*(xPA*xPQ*yPB*yQD + xPA*xQC*yPB*yPQ) + v*(-xPA*xQC*yPQ*yQD - xPQ*xQC*yPB*yQD)) + F2*(pq**2 + pq*v*(-xPQ*xQC - yPQ*yQD) + u**2*xPA*xPQ*yPB*yPQ + u*(pq*(xPA*xPQ + yPB*yPQ) + v*(-xPA*xPQ*yPQ*yQD - xPA*xQC*yPQ**2 - xPQ**2*yPB*yQD - xPQ*xQC*yPB*yPQ)) + v**2*xPQ*xQC*yPQ*yQD) + F3*(u**2*v*(-xPA*xPQ*yPQ**2 - xPQ**2*yPB*yPQ) + u*(pq*v*(-xPQ**2 - yPQ**2) + v**2*(xPQ**2*yPQ*yQD + xPQ*xQC*yPQ**2))) + F4*u**2*v**2*xPQ**2*yPQ**2 + xPA*xQC*yPB*yQD * F0 
      end if
    
      ! Combination 13: p_x p_y p_x p_z
      if (o1 == "px" .and. o2 == "py" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*yPB*zQD + u*(xPA*xPQ*yPB*zQD + xPA*xQC*yPB*zPQ) + v*(-xPA*xQC*yPQ*zQD - xPQ*xQC*yPB*zQD)) + F2*(-pq*v*yPQ*zQD + u**2*xPA*xPQ*yPB*zPQ + u*(pq*yPB*zPQ + v*(-xPA*xPQ*yPQ*zQD - xPA*xQC*yPQ*zPQ - xPQ**2*yPB*zQD - xPQ*xQC*yPB*zPQ)) + v**2*xPQ*xQC*yPQ*zQD) + F3*(u**2*v*(-xPA*xPQ*yPQ*zPQ - xPQ**2*yPB*zPQ) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yPQ*zQD + xPQ*xQC*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPA*xQC*yPB*zQD * F0 
      end if
    
      ! Combination 14: p_x p_z p_x p_y
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*yQD*zPB + u*(xPA*xPQ*yQD*zPB + xPA*xQC*yPQ*zPB) + v*(-xPA*xQC*yQD*zPQ - xPQ*xQC*yQD*zPB)) + F2*(-pq*v*yQD*zPQ + u**2*xPA*xPQ*yPQ*zPB + u*(pq*yPQ*zPB + v*(-xPA*xPQ*yQD*zPQ - xPA*xQC*yPQ*zPQ - xPQ**2*yQD*zPB - xPQ*xQC*yPQ*zPB)) + v**2*xPQ*xQC*yQD*zPQ) + F3*(u**2*v*(-xPA*xPQ*yPQ*zPQ - xPQ**2*yPQ*zPB) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yQD*zPQ + xPQ*xQC*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPA*xQC*yQD*zPB * F0 
      end if
    
      ! Combination 15: p_x p_z p_x p_z
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*(xPA*xQC + zPB*zQD) + u*(xPA*xPQ*zPB*zQD + xPA*xQC*zPB*zPQ) + v*(-xPA*xQC*zPQ*zQD - xPQ*xQC*zPB*zQD)) + F2*(pq**2 + pq*v*(-xPQ*xQC - zPQ*zQD) + u**2*xPA*xPQ*zPB*zPQ + u*(pq*(xPA*xPQ + zPB*zPQ) + v*(-xPA*xPQ*zPQ*zQD - xPA*xQC*zPQ**2 - xPQ**2*zPB*zQD - xPQ*xQC*zPB*zPQ)) + v**2*xPQ*xQC*zPQ*zQD) + F3*(u**2*v*(-xPA*xPQ*zPQ**2 - xPQ**2*zPB*zPQ) + u*(pq*v*(-xPQ**2 - zPQ**2) + v**2*(xPQ**2*zPQ*zQD + xPQ*xQC*zPQ**2))) + F4*u**2*v**2*xPQ**2*zPQ**2 + xPA*xQC*zPB*zQD * F0 
      end if
    
      ! Combination 16: p_x p_y p_y p_x
      if (o1 == "px" .and. o2 == "py" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*(xPA*xQD + yPB*yQC) + u*(xPA*xPQ*yPB*yQC + xPA*xQD*yPB*yPQ) + v*(-xPA*xQD*yPQ*yQC - xPQ*xQD*yPB*yQC)) + F2*(pq**2 + pq*v*(-xPQ*xQD - yPQ*yQC) + u**2*xPA*xPQ*yPB*yPQ + u*(pq*(xPA*xPQ + yPB*yPQ) + v*(-xPA*xPQ*yPQ*yQC - xPA*xQD*yPQ**2 - xPQ**2*yPB*yQC - xPQ*xQD*yPB*yPQ)) + v**2*xPQ*xQD*yPQ*yQC) + F3*(u**2*v*(-xPA*xPQ*yPQ**2 - xPQ**2*yPB*yPQ) + u*(pq*v*(-xPQ**2 - yPQ**2) + v**2*(xPQ**2*yPQ*yQC + xPQ*xQD*yPQ**2))) + F4*u**2*v**2*xPQ**2*yPQ**2 + xPA*xQD*yPB*yQC * F0 
      end if
    
      ! Combination 17: p_x p_y p_z p_x
      if (o1 == "px" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*yPB*zQC + u*(xPA*xPQ*yPB*zQC + xPA*xQD*yPB*zPQ) + v*(-xPA*xQD*yPQ*zQC - xPQ*xQD*yPB*zQC)) + F2*(-pq*v*yPQ*zQC + u**2*xPA*xPQ*yPB*zPQ + u*(pq*yPB*zPQ + v*(-xPA*xPQ*yPQ*zQC - xPA*xQD*yPQ*zPQ - xPQ**2*yPB*zQC - xPQ*xQD*yPB*zPQ)) + v**2*xPQ*xQD*yPQ*zQC) + F3*(u**2*v*(-xPA*xPQ*yPQ*zPQ - xPQ**2*yPB*zPQ) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yPQ*zQC + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPA*xQD*yPB*zQC * F0 
      end if
    
      ! Combination 18: p_x p_z p_y p_x
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*yQC*zPB + u*(xPA*xPQ*yQC*zPB + xPA*xQD*yPQ*zPB) + v*(-xPA*xQD*yQC*zPQ - xPQ*xQD*yQC*zPB)) + F2*(-pq*v*yQC*zPQ + u**2*xPA*xPQ*yPQ*zPB + u*(pq*yPQ*zPB + v*(-xPA*xPQ*yQC*zPQ - xPA*xQD*yPQ*zPQ - xPQ**2*yQC*zPB - xPQ*xQD*yPQ*zPB)) + v**2*xPQ*xQD*yQC*zPQ) + F3*(u**2*v*(-xPA*xPQ*yPQ*zPQ - xPQ**2*yPQ*zPB) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yQC*zPQ + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPA*xQD*yQC*zPB * F0 
      end if
    
      ! Combination 19: p_x p_z p_z p_x
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*(xPA*xQD + zPB*zQC) + u*(xPA*xPQ*zPB*zQC + xPA*xQD*zPB*zPQ) + v*(-xPA*xQD*zPQ*zQC - xPQ*xQD*zPB*zQC)) + F2*(pq**2 + pq*v*(-xPQ*xQD - zPQ*zQC) + u**2*xPA*xPQ*zPB*zPQ + u*(pq*(xPA*xPQ + zPB*zPQ) + v*(-xPA*xPQ*zPQ*zQC - xPA*xQD*zPQ**2 - xPQ**2*zPB*zQC - xPQ*xQD*zPB*zPQ)) + v**2*xPQ*xQD*zPQ*zQC) + F3*(u**2*v*(-xPA*xPQ*zPQ**2 - xPQ**2*zPB*zPQ) + u*(pq*v*(-xPQ**2 - zPQ**2) + v**2*(xPQ**2*zPQ*zQC + xPQ*xQD*zPQ**2))) + F4*u**2*v**2*xPQ**2*zPQ**2 + xPA*xQD*zPB*zQC * F0 
      end if
    
      ! Combination 20: p_x p_y p_y p_y
      if (o1 == "px" .and. o2 == "py" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(xPA*yQC + xPA*yQD) - tv*xPA*yPB + u*(xPA*yPB*yPQ*yQC + xPA*yPB*yPQ*yQD) + v*(-xPA*yPQ*yQC*yQD - xPQ*yPB*yQC*yQD - xPA*yPQ/(2*q) - xPQ*yPB/(2*q))) + F2*(u**2*xPA*yPB*yPQ**2 + u*(2*pq*xPA*yPQ + v*(-xPA*yPQ**2*yQC - xPA*yPQ**2*yQD - xPQ*yPB*yPQ*yQC - xPQ*yPB*yPQ*yQD)) + v**2*(xPQ*yPQ*yQC*yQD + xPQ*yPQ/(2*q)) + v*(pq*(-xPQ*yQC - xPQ*yQD) + tv*xPA*yPQ + tv*xPQ*yPB)) + F3*(-tv*v**2*xPQ*yPQ + u**2*v*(-xPA*yPQ**3 - xPQ*yPB*yPQ**2) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ*yPQ**2*yQC + xPQ*yPQ**2*yQD))) + F4*u**2*v**2*xPQ*yPQ**3 + (xPA*yPB*yQC*yQD + xPA*yPB/(2*q)) * F0 
      end if
    
      ! Combination 21: p_x p_y p_y p_z
      if (o1 == "px" .and. o2 == "py" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*xPA*zQD + u*(xPA*yPB*yPQ*zQD + xPA*yPB*yQC*zPQ) + v*(-xPA*yPQ*yQC*zQD - xPQ*yPB*yQC*zQD)) + F2*(-pq*v*xPQ*zQD + u**2*xPA*yPB*yPQ*zPQ + u*(pq*xPA*zPQ + v*(-xPA*yPQ**2*zQD - xPA*yPQ*yQC*zPQ - xPQ*yPB*yPQ*zQD - xPQ*yPB*yQC*zPQ)) + v**2*xPQ*yPQ*yQC*zQD) + F3*(u**2*v*(-xPA*yPQ**2*zPQ - xPQ*yPB*yPQ*zPQ) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ**2*zQD + xPQ*yPQ*yQC*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xPA*yPB*yQC*zQD * F0 
      end if
    
      ! Combination 22: p_x p_y p_z p_y
      if (o1 == "px" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*xPA*zQC + u*(xPA*yPB*yPQ*zQC + xPA*yPB*yQD*zPQ) + v*(-xPA*yPQ*yQD*zQC - xPQ*yPB*yQD*zQC)) + F2*(-pq*v*xPQ*zQC + u**2*xPA*yPB*yPQ*zPQ + u*(pq*xPA*zPQ + v*(-xPA*yPQ**2*zQC - xPA*yPQ*yQD*zPQ - xPQ*yPB*yPQ*zQC - xPQ*yPB*yQD*zPQ)) + v**2*xPQ*yPQ*yQD*zQC) + F3*(u**2*v*(-xPA*yPQ**2*zPQ - xPQ*yPB*yPQ*zPQ) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ**2*zQC + xPQ*yPQ*yQD*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xPA*yPB*yQD*zQC * F0 
      end if
    
      ! Combination 23: p_x p_y p_z p_z
      if (o1 == "px" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*xPA*yPB + u*(xPA*yPB*zPQ*zQC + xPA*yPB*zPQ*zQD) + v*(-xPA*yPQ*zQC*zQD - xPQ*yPB*zQC*zQD - xPA*yPQ/(2*q) - xPQ*yPB/(2*q))) + F2*(u**2*xPA*yPB*zPQ**2 + v**2*(xPQ*yPQ*zQC*zQD + xPQ*yPQ/(2*q)) + v*(tv*xPA*yPQ + tv*xPQ*yPB + u*(-xPA*yPQ*zPQ*zQC - xPA*yPQ*zPQ*zQD - xPQ*yPB*zPQ*zQC - xPQ*yPB*zPQ*zQD))) + F3*(u**2*v*(-xPA*yPQ*zPQ**2 - xPQ*yPB*zPQ**2) + v**2*(-tv*xPQ*yPQ + u*(xPQ*yPQ*zPQ*zQC + xPQ*yPQ*zPQ*zQD))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + (xPA*yPB*zQC*zQD + xPA*yPB/(2*q)) * F0 
      end if
    
      ! Combination 24: p_x p_z p_y p_y
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*xPA*zPB + u*(xPA*yPQ*yQC*zPB + xPA*yPQ*yQD*zPB) + v*(-xPA*yQC*yQD*zPQ - xPQ*yQC*yQD*zPB - xPA*zPQ/(2*q) - xPQ*zPB/(2*q))) + F2*(u**2*xPA*yPQ**2*zPB + v**2*(xPQ*yQC*yQD*zPQ + xPQ*zPQ/(2*q)) + v*(tv*xPA*zPQ + tv*xPQ*zPB + u*(-xPA*yPQ*yQC*zPQ - xPA*yPQ*yQD*zPQ - xPQ*yPQ*yQC*zPB - xPQ*yPQ*yQD*zPB))) + F3*(u**2*v*(-xPA*yPQ**2*zPQ - xPQ*yPQ**2*zPB) + v**2*(-tv*xPQ*zPQ + u*(xPQ*yPQ*yQC*zPQ + xPQ*yPQ*yQD*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + (xPA*yQC*yQD*zPB + xPA*zPB/(2*q)) * F0 
      end if
    
      ! Combination 25: p_x p_z p_y p_z
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*xPA*yQC + u*(xPA*yPQ*zPB*zQD + xPA*yQC*zPB*zPQ) + v*(-xPA*yQC*zPQ*zQD - xPQ*yQC*zPB*zQD)) + F2*(-pq*v*xPQ*yQC + u**2*xPA*yPQ*zPB*zPQ + u*(pq*xPA*yPQ + v*(-xPA*yPQ*zPQ*zQD - xPA*yQC*zPQ**2 - xPQ*yPQ*zPB*zQD - xPQ*yQC*zPB*zPQ)) + v**2*xPQ*yQC*zPQ*zQD) + F3*(u**2*v*(-xPA*yPQ*zPQ**2 - xPQ*yPQ*zPB*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQD + xPQ*yQC*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xPA*yQC*zPB*zQD * F0 
      end if
    
      ! Combination 26: p_x p_z p_z p_y
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*xPA*yQD + u*(xPA*yPQ*zPB*zQC + xPA*yQD*zPB*zPQ) + v*(-xPA*yQD*zPQ*zQC - xPQ*yQD*zPB*zQC)) + F2*(-pq*v*xPQ*yQD + u**2*xPA*yPQ*zPB*zPQ + u*(pq*xPA*yPQ + v*(-xPA*yPQ*zPQ*zQC - xPA*yQD*zPQ**2 - xPQ*yPQ*zPB*zQC - xPQ*yQD*zPB*zPQ)) + v**2*xPQ*yQD*zPQ*zQC) + F3*(u**2*v*(-xPA*yPQ*zPQ**2 - xPQ*yPQ*zPB*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQC + xPQ*yQD*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xPA*yQD*zPB*zQC * F0 
      end if
    
      ! Combination 27: p_x p_z p_z p_z
      if (o1 == "px" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(xPA*zQC + xPA*zQD) - tv*xPA*zPB + u*(xPA*zPB*zPQ*zQC + xPA*zPB*zPQ*zQD) + v*(-xPA*zPQ*zQC*zQD - xPQ*zPB*zQC*zQD - xPA*zPQ/(2*q) - xPQ*zPB/(2*q))) + F2*(u**2*xPA*zPB*zPQ**2 + u*(2*pq*xPA*zPQ + v*(-xPA*zPQ**2*zQC - xPA*zPQ**2*zQD - xPQ*zPB*zPQ*zQC - xPQ*zPB*zPQ*zQD)) + v**2*(xPQ*zPQ*zQC*zQD + xPQ*zPQ/(2*q)) + v*(pq*(-xPQ*zQC - xPQ*zQD) + tv*xPA*zPQ + tv*xPQ*zPB)) + F3*(-tv*v**2*xPQ*zPQ + u**2*v*(-xPA*zPQ**3 - xPQ*zPB*zPQ**2) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ*zPQ**2*zQC + xPQ*zPQ**2*zQD))) + F4*u**2*v**2*xPQ*zPQ**3 + (xPA*zPB*zQC*zQD + xPA*zPB/(2*q)) * F0 
      end if
    
      ! Combination 28: p_y p_x p_x p_x
      if (o1 == "py" .and. o2 == "px" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC*yPA + xQD*yPA) - tv*xPB*yPA + u*(xPB*xPQ*xQC*yPA + xPB*xPQ*xQD*yPA) + v*(-xPB*xQC*xQD*yPQ - xPQ*xQC*xQD*yPA - xPB*yPQ/(2*q) - xPQ*yPA/(2*q))) + F2*(u**2*xPB*xPQ**2*yPA + u*(2*pq*xPQ*yPA + v*(-xPB*xPQ*xQC*yPQ - xPB*xPQ*xQD*yPQ - xPQ**2*xQC*yPA - xPQ**2*xQD*yPA)) + v**2*(xPQ*xQC*xQD*yPQ + xPQ*yPQ/(2*q)) + v*(pq*(-xQC*yPQ - xQD*yPQ) + tv*xPB*yPQ + tv*xPQ*yPA)) + F3*(-tv*v**2*xPQ*yPQ + u**2*v*(-xPB*xPQ**2*yPQ - xPQ**3*yPA) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ**2*xQC*yPQ + xPQ**2*xQD*yPQ))) + F4*u**2*v**2*xPQ**3*yPQ + (xPB*xQC*xQD*yPA + xPB*yPA/(2*q)) * F0 
      end if
    
      ! Combination 29: p_z p_x p_x p_x
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(pq*(xQC*zPA + xQD*zPA) - tv*xPB*zPA + u*(xPB*xPQ*xQC*zPA + xPB*xPQ*xQD*zPA) + v*(-xPB*xQC*xQD*zPQ - xPQ*xQC*xQD*zPA - xPB*zPQ/(2*q) - xPQ*zPA/(2*q))) + F2*(u**2*xPB*xPQ**2*zPA + u*(2*pq*xPQ*zPA + v*(-xPB*xPQ*xQC*zPQ - xPB*xPQ*xQD*zPQ - xPQ**2*xQC*zPA - xPQ**2*xQD*zPA)) + v**2*(xPQ*xQC*xQD*zPQ + xPQ*zPQ/(2*q)) + v*(pq*(-xQC*zPQ - xQD*zPQ) + tv*xPB*zPQ + tv*xPQ*zPA)) + F3*(-tv*v**2*xPQ*zPQ + u**2*v*(-xPB*xPQ**2*zPQ - xPQ**3*zPA) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ**2*xQC*zPQ + xPQ**2*xQD*zPQ))) + F4*u**2*v**2*xPQ**3*zPQ + (xPB*xQC*xQD*zPA + xPB*zPA/(2*q)) * F0 
      end if
    
      ! Combination 30: p_y p_x p_x p_y
      if (o1 == "py" .and. o2 == "px" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*(xPB*xQC + yPA*yQD) + u*(xPB*xPQ*yPA*yQD + xPB*xQC*yPA*yPQ) + v*(-xPB*xQC*yPQ*yQD - xPQ*xQC*yPA*yQD)) + F2*(pq**2 + pq*v*(-xPQ*xQC - yPQ*yQD) + u**2*xPB*xPQ*yPA*yPQ + u*(pq*(xPB*xPQ + yPA*yPQ) + v*(-xPB*xPQ*yPQ*yQD - xPB*xQC*yPQ**2 - xPQ**2*yPA*yQD - xPQ*xQC*yPA*yPQ)) + v**2*xPQ*xQC*yPQ*yQD) + F3*(u**2*v*(-xPB*xPQ*yPQ**2 - xPQ**2*yPA*yPQ) + u*(pq*v*(-xPQ**2 - yPQ**2) + v**2*(xPQ**2*yPQ*yQD + xPQ*xQC*yPQ**2))) + F4*u**2*v**2*xPQ**2*yPQ**2 + xPB*xQC*yPA*yQD * F0 
      end if
    
      ! Combination 31: p_y p_x p_x p_z
      if (o1 == "py" .and. o2 == "px" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*yPA*zQD + u*(xPB*xPQ*yPA*zQD + xPB*xQC*yPA*zPQ) + v*(-xPB*xQC*yPQ*zQD - xPQ*xQC*yPA*zQD)) + F2*(-pq*v*yPQ*zQD + u**2*xPB*xPQ*yPA*zPQ + u*(pq*yPA*zPQ + v*(-xPB*xPQ*yPQ*zQD - xPB*xQC*yPQ*zPQ - xPQ**2*yPA*zQD - xPQ*xQC*yPA*zPQ)) + v**2*xPQ*xQC*yPQ*zQD) + F3*(u**2*v*(-xPB*xPQ*yPQ*zPQ - xPQ**2*yPA*zPQ) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yPQ*zQD + xPQ*xQC*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPB*xQC*yPA*zQD * F0 
      end if 
    
      ! Combination 32: p_z p_x p_x p_y
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*yQD*zPA + u*(xPB*xPQ*yQD*zPA + xPB*xQC*yPQ*zPA) + v*(-xPB*xQC*yQD*zPQ - xPQ*xQC*yQD*zPA)) + F2*(-pq*v*yQD*zPQ + u**2*xPB*xPQ*yPQ*zPA + u*(pq*yPQ*zPA + v*(-xPB*xPQ*yQD*zPQ - xPB*xQC*yPQ*zPQ - xPQ**2*yQD*zPA - xPQ*xQC*yPQ*zPA)) + v**2*xPQ*xQC*yQD*zPQ) + F3*(u**2*v*(-xPB*xPQ*yPQ*zPQ - xPQ**2*yPQ*zPA) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yQD*zPQ + xPQ*xQC*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPB*xQC*yQD*zPA * F0 
      end if
    
      ! Combination 33: p_z p_x p_x p_z
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*(xPB*xQC + zPA*zQD) + u*(xPB*xPQ*zPA*zQD + xPB*xQC*zPA*zPQ) + v*(-xPB*xQC*zPQ*zQD - xPQ*xQC*zPA*zQD)) + F2*(pq**2 + pq*v*(-xPQ*xQC - zPQ*zQD) + u**2*xPB*xPQ*zPA*zPQ + u*(pq*(xPB*xPQ + zPA*zPQ) + v*(-xPB*xPQ*zPQ*zQD - xPB*xQC*zPQ**2 - xPQ**2*zPA*zQD - xPQ*xQC*zPA*zPQ)) + v**2*xPQ*xQC*zPQ*zQD) + F3*(u**2*v*(-xPB*xPQ*zPQ**2 - xPQ**2*zPA*zPQ) + u*(pq*v*(-xPQ**2 - zPQ**2) + v**2*(xPQ**2*zPQ*zQD + xPQ*xQC*zPQ**2))) + F4*u**2*v**2*xPQ**2*zPQ**2 + xPB*xQC*zPA*zQD * F0 
      end if
    
      ! Combination 34: p_y p_x p_y p_x
      if (o1 == "py" .and. o2 == "px" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*(xPB*xQD + yPA*yQC) + u*(xPB*xPQ*yPA*yQC + xPB*xQD*yPA*yPQ) + v*(-xPB*xQD*yPQ*yQC - xPQ*xQD*yPA*yQC)) + F2*(pq**2 + pq*v*(-xPQ*xQD - yPQ*yQC) + u**2*xPB*xPQ*yPA*yPQ + u*(pq*(xPB*xPQ + yPA*yPQ) + v*(-xPB*xPQ*yPQ*yQC - xPB*xQD*yPQ**2 - xPQ**2*yPA*yQC - xPQ*xQD*yPA*yPQ)) + v**2*xPQ*xQD*yPQ*yQC) + F3*(u**2*v*(-xPB*xPQ*yPQ**2 - xPQ**2*yPA*yPQ) + u*(pq*v*(-xPQ**2 - yPQ**2) + v**2*(xPQ**2*yPQ*yQC + xPQ*xQD*yPQ**2))) + F4*u**2*v**2*xPQ**2*yPQ**2 + xPB*xQD*yPA*yQC * F0 
      end if
    
      ! Combination 35: p_y p_x p_z p_x
      if (o1 == "py" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*yPA*zQC + u*(xPB*xPQ*yPA*zQC + xPB*xQD*yPA*zPQ) + v*(-xPB*xQD*yPQ*zQC - xPQ*xQD*yPA*zQC)) + F2*(-pq*v*yPQ*zQC + u**2*xPB*xPQ*yPA*zPQ + u*(pq*yPA*zPQ + v*(-xPB*xPQ*yPQ*zQC - xPB*xQD*yPQ*zPQ - xPQ**2*yPA*zQC - xPQ*xQD*yPA*zPQ)) + v**2*xPQ*xQD*yPQ*zQC) + F3*(u**2*v*(-xPB*xPQ*yPQ*zPQ - xPQ**2*yPA*zPQ) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yPQ*zQC + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPB*xQD*yPA*zQC * F0 
      end if
    
      ! Combination 36: p_z p_x p_y p_x
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*yQC*zPA + u*(xPB*xPQ*yQC*zPA + xPB*xQD*yPQ*zPA) + v*(-xPB*xQD*yQC*zPQ - xPQ*xQD*yQC*zPA)) + F2*(-pq*v*yQC*zPQ + u**2*xPB*xPQ*yPQ*zPA + u*(pq*yPQ*zPA + v*(-xPB*xPQ*yQC*zPQ - xPB*xQD*yPQ*zPQ - xPQ**2*yQC*zPA - xPQ*xQD*yPQ*zPA)) + v**2*xPQ*xQD*yQC*zPQ) + F3*(u**2*v*(-xPB*xPQ*yPQ*zPQ - xPQ**2*yPQ*zPA) + u*(-pq*v*yPQ*zPQ + v**2*(xPQ**2*yQC*zPQ + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + xPB*xQD*yQC*zPA * F0 
      end if
    
      ! Combination 37: p_z p_x p_z p_x
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*(xPB*xQD + zPA*zQC) + u*(xPB*xPQ*zPA*zQC + xPB*xQD*zPA*zPQ) + v*(-xPB*xQD*zPQ*zQC - xPQ*xQD*zPA*zQC)) + F2*(pq**2 + pq*v*(-xPQ*xQD - zPQ*zQC) + u**2*xPB*xPQ*zPA*zPQ + u*(pq*(xPB*xPQ + zPA*zPQ) + v*(-xPB*xPQ*zPQ*zQC - xPB*xQD*zPQ**2 - xPQ**2*zPA*zQC - xPQ*xQD*zPA*zPQ)) + v**2*xPQ*xQD*zPQ*zQC) + F3*(u**2*v*(-xPB*xPQ*zPQ**2 - xPQ**2*zPA*zPQ) + u*(pq*v*(-xPQ**2 - zPQ**2) + v**2*(xPQ**2*zPQ*zQC + xPQ*xQD*zPQ**2))) + F4*u**2*v**2*xPQ**2*zPQ**2 + xPB*xQD*zPA*zQC * F0 
      end if
    
      ! Combination 38: p_y p_x p_y p_y
      if (o1 == "py" .and. o2 == "px" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(xPB*yQC + xPB*yQD) - tv*xPB*yPA + u*(xPB*yPA*yPQ*yQC + xPB*yPA*yPQ*yQD) + v*(-xPB*yPQ*yQC*yQD - xPQ*yPA*yQC*yQD - xPB*yPQ/(2*q) - xPQ*yPA/(2*q))) + F2*(u**2*xPB*yPA*yPQ**2 + u*(2*pq*xPB*yPQ + v*(-xPB*yPQ**2*yQC - xPB*yPQ**2*yQD - xPQ*yPA*yPQ*yQC - xPQ*yPA*yPQ*yQD)) + v**2*(xPQ*yPQ*yQC*yQD + xPQ*yPQ/(2*q)) + v*(pq*(-xPQ*yQC - xPQ*yQD) + tv*xPB*yPQ + tv*xPQ*yPA)) + F3*(-tv*v**2*xPQ*yPQ + u**2*v*(-xPB*yPQ**3 - xPQ*yPA*yPQ**2) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ*yPQ**2*yQC + xPQ*yPQ**2*yQD))) + F4*u**2*v**2*xPQ*yPQ**3 + (xPB*yPA*yQC*yQD + xPB*yPA/(2*q)) * F0 
      end if
    
      ! Combination 39: p_y p_x p_y p_z
      if (o1 == "py" .and. o2 == "px" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*xPB*zQD + u*(xPB*yPA*yPQ*zQD + xPB*yPA*yQC*zPQ) + v*(-xPB*yPQ*yQC*zQD - xPQ*yPA*yQC*zQD)) + F2*(-pq*v*xPQ*zQD + u**2*xPB*yPA*yPQ*zPQ + u*(pq*xPB*zPQ + v*(-xPB*yPQ**2*zQD - xPB*yPQ*yQC*zPQ - xPQ*yPA*yPQ*zQD - xPQ*yPA*yQC*zPQ)) + v**2*xPQ*yPQ*yQC*zQD) + F3*(u**2*v*(-xPB*yPQ**2*zPQ - xPQ*yPA*yPQ*zPQ) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ**2*zQD + xPQ*yPQ*yQC*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xPB*yPA*yQC*zQD * F0 
      end if
    
      ! Combination 40: p_y p_x p_z p_y
      if (o1 == "py" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*xPB*zQC + u*(xPB*yPA*yPQ*zQC + xPB*yPA*yQD*zPQ) + v*(-xPB*yPQ*yQD*zQC - xPQ*yPA*yQD*zQC)) + F2*(-pq*v*xPQ*zQC + u**2*xPB*yPA*yPQ*zPQ + u*(pq*xPB*zPQ + v*(-xPB*yPQ**2*zQC - xPB*yPQ*yQD*zPQ - xPQ*yPA*yPQ*zQC - xPQ*yPA*yQD*zPQ)) + v**2*xPQ*yPQ*yQD*zQC) + F3*(u**2*v*(-xPB*yPQ**2*zPQ - xPQ*yPA*yPQ*zPQ) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ**2*zQC + xPQ*yPQ*yQD*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xPB*yPA*yQD*zQC * F0 
      end if
    
      ! Combination 41: p_y p_x p_z p_z
      if (o1 == "py" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tv*xPB*yPA + u*(xPB*yPA*zPQ*zQC + xPB*yPA*zPQ*zQD) + v*(-xPB*yPQ*zQC*zQD - xPQ*yPA*zQC*zQD - xPB*yPQ/(2*q) - xPQ*yPA/(2*q))) + F2*(u**2*xPB*yPA*zPQ**2 + v**2*(xPQ*yPQ*zQC*zQD + xPQ*yPQ/(2*q)) + v*(tv*xPB*yPQ + tv*xPQ*yPA + u*(-xPB*yPQ*zPQ*zQC - xPB*yPQ*zPQ*zQD - xPQ*yPA*zPQ*zQC - xPQ*yPA*zPQ*zQD))) + F3*(u**2*v*(-xPB*yPQ*zPQ**2 - xPQ*yPA*zPQ**2) + v**2*(-tv*xPQ*yPQ + u*(xPQ*yPQ*zPQ*zQC + xPQ*yPQ*zPQ*zQD))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + (xPB*yPA*zQC*zQD + xPB*yPA/(2*q)) * F0 
      end if
    
      ! Combination 42: p_z p_x p_y p_y
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tv*xPB*zPA + u*(xPB*yPQ*yQC*zPA + xPB*yPQ*yQD*zPA) + v*(-xPB*yQC*yQD*zPQ - xPQ*yQC*yQD*zPA - xPB*zPQ/(2*q) - xPQ*zPA/(2*q))) + F2*(u**2*xPB*yPQ**2*zPA + v**2*(xPQ*yQC*yQD*zPQ + xPQ*zPQ/(2*q)) + v*(tv*xPB*zPQ + tv*xPQ*zPA + u*(-xPB*yPQ*yQC*zPQ - xPB*yPQ*yQD*zPQ - xPQ*yPQ*yQC*zPA - xPQ*yPQ*yQD*zPA))) + F3*(u**2*v*(-xPB*yPQ**2*zPQ - xPQ*yPQ**2*zPA) + v**2*(-tv*xPQ*zPQ + u*(xPQ*yPQ*yQC*zPQ + xPQ*yPQ*yQD*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + (xPB*yQC*yQD*zPA + xPB*zPA/(2*q)) * F0 
      end if
    
      ! Combination 43: p_z p_x p_y p_z
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*xPB*yQC + u*(xPB*yPQ*zPA*zQD + xPB*yQC*zPA*zPQ) + v*(-xPB*yQC*zPQ*zQD - xPQ*yQC*zPA*zQD)) + F2*(-pq*v*xPQ*yQC + u**2*xPB*yPQ*zPA*zPQ + u*(pq*xPB*yPQ + v*(-xPB*yPQ*zPQ*zQD - xPB*yQC*zPQ**2 - xPQ*yPQ*zPA*zQD - xPQ*yQC*zPA*zPQ)) + v**2*xPQ*yQC*zPQ*zQD) + F3*(u**2*v*(-xPB*yPQ*zPQ**2 - xPQ*yPQ*zPA*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQD + xPQ*yQC*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xPB*yQC*zPA*zQD * F0 
      end if
    
      ! Combination 44: p_z p_x p_z p_y
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*xPB*yQD + u*(xPB*yPQ*zPA*zQC + xPB*yQD*zPA*zPQ) + v*(-xPB*yQD*zPQ*zQC - xPQ*yQD*zPA*zQC)) + F2*(-pq*v*xPQ*yQD + u**2*xPB*yPQ*zPA*zPQ + u*(pq*xPB*yPQ + v*(-xPB*yPQ*zPQ*zQC - xPB*yQD*zPQ**2 - xPQ*yPQ*zPA*zQC - xPQ*yQD*zPA*zPQ)) + v**2*xPQ*yQD*zPQ*zQC) + F3*(u**2*v*(-xPB*yPQ*zPQ**2 - xPQ*yPQ*zPA*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQC + xPQ*yQD*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xPB*yQD*zPA*zQC * F0 
      end if
    
      ! Combination 45: p_z p_x p_z p_z
      if (o1 == "pz" .and. o2 == "px" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(xPB*zQC + xPB*zQD) - tv*xPB*zPA + u*(xPB*zPA*zPQ*zQC + xPB*zPA*zPQ*zQD) + v*(-xPB*zPQ*zQC*zQD - xPQ*zPA*zQC*zQD - xPB*zPQ/(2*q) - xPQ*zPA/(2*q))) + F2*(u**2*xPB*zPA*zPQ**2 + u*(2*pq*xPB*zPQ + v*(-xPB*zPQ**2*zQC - xPB*zPQ**2*zQD - xPQ*zPA*zPQ*zQC - xPQ*zPA*zPQ*zQD)) + v**2*(xPQ*zPQ*zQC*zQD + xPQ*zPQ/(2*q)) + v*(pq*(-xPQ*zQC - xPQ*zQD) + tv*xPB*zPQ + tv*xPQ*zPA)) + F3*(-tv*v**2*xPQ*zPQ + u**2*v*(-xPB*zPQ**3 - xPQ*zPA*zPQ**2) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ*zPQ**2*zQC + xPQ*zPQ**2*zQD))) + F4*u**2*v**2*xPQ*zPQ**3 + (xPB*zPA*zQC*zQD + xPB*zPA/(2*q)) * F0 
      end if
    
      ! Combination 46: p_y p_y p_x p_x
      if (o1 == "py" .and. o2 == "py" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tu*xQC*xQD - tv*yPA*yPB + u*(xPQ*xQC*yPA*yPB + xPQ*xQD*yPA*yPB + xPQ*xQC/(2*p) + xPQ*xQD/(2*p)) + v*(-xQC*xQD*yPA*yPQ - xQC*xQD*yPB*yPQ - yPA*yPQ/(2*q) - yPB*yPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(xPQ**2*yPA*yPB + xPQ**2/(2*p)) + u*(-tu*xPQ*xQC - tu*xPQ*xQD + v*(-xPQ*xQC*yPA*yPQ - xPQ*xQC*yPB*yPQ - xPQ*xQD*yPA*yPQ - xPQ*xQD*yPB*yPQ)) + v**2*(xQC*xQD*yPQ**2 + yPQ**2/(2*q)) + v*(tv*yPA*yPQ + tv*yPB*yPQ)) + F3*(u**2*(-tu*xPQ**2 + v*(-xPQ**2*yPA*yPQ - xPQ**2*yPB*yPQ)) + v**2*(-tv*yPQ**2 + u*(xPQ*xQC*yPQ**2 + xPQ*xQD*yPQ**2))) + F4*u**2*v**2*xPQ**2*yPQ**2 + (xQC*xQD*yPA*yPB + yPA*yPB/(2*q) + xQC*xQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 47: p_y p_z p_x p_x
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*yPA*zPB + u*(xPQ*xQC*yPA*zPB + xPQ*xQD*yPA*zPB) + v*(-xQC*xQD*yPA*zPQ - xQC*xQD*yPQ*zPB - yPA*zPQ/(2*q) - yPQ*zPB/(2*q))) + F2*(u**2*xPQ**2*yPA*zPB + v**2*(xQC*xQD*yPQ*zPQ + yPQ*zPQ/(2*q)) + v*(tv*yPA*zPQ + tv*yPQ*zPB + u*(-xPQ*xQC*yPA*zPQ - xPQ*xQC*yPQ*zPB - xPQ*xQD*yPA*zPQ - xPQ*xQD*yPQ*zPB))) + F3*(u**2*v*(-xPQ**2*yPA*zPQ - xPQ**2*yPQ*zPB) + v**2*(-tv*yPQ*zPQ + u*(xPQ*xQC*yPQ*zPQ + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + (xQC*xQD*yPA*zPB + yPA*zPB/(2*q)) * F0 
      end if
    
      ! Combination 48: p_z p_y p_x p_x
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tv*yPB*zPA + u*(xPQ*xQC*yPB*zPA + xPQ*xQD*yPB*zPA) + v*(-xQC*xQD*yPB*zPQ - xQC*xQD*yPQ*zPA - yPB*zPQ/(2*q) - yPQ*zPA/(2*q))) + F2*(u**2*xPQ**2*yPB*zPA + v**2*(xQC*xQD*yPQ*zPQ + yPQ*zPQ/(2*q)) + v*(tv*yPB*zPQ + tv*yPQ*zPA + u*(-xPQ*xQC*yPB*zPQ - xPQ*xQC*yPQ*zPA - xPQ*xQD*yPB*zPQ - xPQ*xQD*yPQ*zPA))) + F3*(u**2*v*(-xPQ**2*yPB*zPQ - xPQ**2*yPQ*zPA) + v**2*(-tv*yPQ*zPQ + u*(xPQ*xQC*yPQ*zPQ + xPQ*xQD*yPQ*zPQ))) + F4*u**2*v**2*xPQ**2*yPQ*zPQ + (xQC*xQD*yPB*zPA + yPB*zPA/(2*q)) * F0 
      end if
    
      ! Combination 49: p_z p_z p_x p_x
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "px") then
        value = F1*(-tu*xQC*xQD - tv*zPA*zPB + u*(xPQ*xQC*zPA*zPB + xPQ*xQD*zPA*zPB + xPQ*xQC/(2*p) + xPQ*xQD/(2*p)) + v*(-xQC*xQD*zPA*zPQ - xQC*xQD*zPB*zPQ - zPA*zPQ/(2*q) - zPB*zPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(xPQ**2*zPA*zPB + xPQ**2/(2*p)) + u*(-tu*xPQ*xQC - tu*xPQ*xQD + v*(-xPQ*xQC*zPA*zPQ - xPQ*xQC*zPB*zPQ - xPQ*xQD*zPA*zPQ - xPQ*xQD*zPB*zPQ)) + v**2*(xQC*xQD*zPQ**2 + zPQ**2/(2*q)) + v*(tv*zPA*zPQ + tv*zPB*zPQ)) + F3*(u**2*(-tu*xPQ**2 + v*(-xPQ**2*zPA*zPQ - xPQ**2*zPB*zPQ)) + v**2*(-tv*zPQ**2 + u*(xPQ*xQC*zPQ**2 + xPQ*xQD*zPQ**2))) + F4*u**2*v**2*xPQ**2*zPQ**2 + (xQC*xQD*zPA*zPB + zPA*zPB/(2*q) + xQC*xQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 50: p_y p_y p_x p_y
      if (o1 == "py" .and. o2 == "py" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*(xQC*yPA + xQC*yPB) - tu*xQC*yQD + u*(xPQ*yPA*yPB*yQD + xQC*yPA*yPB*yPQ + xPQ*yQD/(2*p) + xQC*yPQ/(2*p)) + v*(-xQC*yPA*yPQ*yQD - xQC*yPB*yPQ*yQD)) + F2*(-2*pq*v*xQC*yPQ + u**2*(xPQ*yPA*yPB*yPQ + xPQ*yPQ/(2*p)) + u*(pq*(xPQ*yPA + xPQ*yPB) - tu*xPQ*yQD - tu*xQC*yPQ + v*(-xPQ*yPA*yPQ*yQD - xPQ*yPB*yPQ*yQD - xQC*yPA*yPQ**2 - xQC*yPB*yPQ**2)) + v**2*xQC*yPQ**2*yQD) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPQ*yPA*yPQ**2 - xPQ*yPB*yPQ**2)) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ*yPQ**2*yQD + xQC*yPQ**3))) + F4*u**2*v**2*xPQ*yPQ**3 + (xQC*yPA*yPB*yQD + xQC*yQD/(2*p)) * F0 
      end if
    
      ! Combination 51: p_y p_y p_x p_z
      if (o1 == "py" .and. o2 == "py" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(-tu*xQC*zQD + u*(xPQ*yPA*yPB*zQD + xQC*yPA*yPB*zPQ + xPQ*zQD/(2*p) + xQC*zPQ/(2*p)) + v*(-xQC*yPA*yPQ*zQD - xQC*yPB*yPQ*zQD)) + F2*(u**2*(xPQ*yPA*yPB*zPQ + xPQ*zPQ/(2*p)) + u*(-tu*xPQ*zQD - tu*xQC*zPQ + v*(-xPQ*yPA*yPQ*zQD - xPQ*yPB*yPQ*zQD - xQC*yPA*yPQ*zPQ - xQC*yPB*yPQ*zPQ)) + v**2*xQC*yPQ**2*zQD) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPQ*yPA*yPQ*zPQ - xPQ*yPB*yPQ*zPQ)) + u*v**2*(xPQ*yPQ**2*zQD + xQC*yPQ**2*zPQ)) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + (xQC*yPA*yPB*zQD + xQC*zQD/(2*p)) * F0 
      end if
    
      ! Combination 52: p_y p_z p_x p_y
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*xQC*zPB + u*(xPQ*yPA*yQD*zPB + xQC*yPA*yPQ*zPB) + v*(-xQC*yPA*yQD*zPQ - xQC*yPQ*yQD*zPB)) + F2*(-pq*v*xQC*zPQ + u**2*xPQ*yPA*yPQ*zPB + u*(pq*xPQ*zPB + v*(-xPQ*yPA*yQD*zPQ - xPQ*yPQ*yQD*zPB - xQC*yPA*yPQ*zPQ - xQC*yPQ**2*zPB)) + v**2*xQC*yPQ*yQD*zPQ) + F3*(u**2*v*(-xPQ*yPA*yPQ*zPQ - xPQ*yPQ**2*zPB) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ*yQD*zPQ + xQC*yPQ**2*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xQC*yPA*yQD*zPB * F0 
      end if
    
      ! Combination 53: p_y p_z p_x p_z
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*xQC*yPA + u*(xPQ*yPA*zPB*zQD + xQC*yPA*zPB*zPQ) + v*(-xQC*yPA*zPQ*zQD - xQC*yPQ*zPB*zQD)) + F2*(-pq*v*xQC*yPQ + u**2*xPQ*yPA*zPB*zPQ + u*(pq*xPQ*yPA + v*(-xPQ*yPA*zPQ*zQD - xPQ*yPQ*zPB*zQD - xQC*yPA*zPQ**2 - xQC*yPQ*zPB*zPQ)) + v**2*xQC*yPQ*zPQ*zQD) + F3*(u**2*v*(-xPQ*yPA*zPQ**2 - xPQ*yPQ*zPB*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQD + xQC*yPQ*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xQC*yPA*zPB*zQD * F0 
      end if
    
      ! Combination 54: p_z p_y p_x p_y
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(pq*xQC*zPA + u*(xPQ*yPB*yQD*zPA + xQC*yPB*yPQ*zPA) + v*(-xQC*yPB*yQD*zPQ - xQC*yPQ*yQD*zPA)) + F2*(-pq*v*xQC*zPQ + u**2*xPQ*yPB*yPQ*zPA + u*(pq*xPQ*zPA + v*(-xPQ*yPB*yQD*zPQ - xPQ*yPQ*yQD*zPA - xQC*yPB*yPQ*zPQ - xQC*yPQ**2*zPA)) + v**2*xQC*yPQ*yQD*zPQ) + F3*(u**2*v*(-xPQ*yPB*yPQ*zPQ - xPQ*yPQ**2*zPA) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ*yQD*zPQ + xQC*yPQ**2*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xQC*yPB*yQD*zPA * F0 
      end if
    
      ! Combination 55: p_z p_y p_x p_z
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*xQC*yPB + u*(xPQ*yPB*zPA*zQD + xQC*yPB*zPA*zPQ) + v*(-xQC*yPB*zPQ*zQD - xQC*yPQ*zPA*zQD)) + F2*(-pq*v*xQC*yPQ + u**2*xPQ*yPB*zPA*zPQ + u*(pq*xPQ*yPB + v*(-xPQ*yPB*zPQ*zQD - xPQ*yPQ*zPA*zQD - xQC*yPB*zPQ**2 - xQC*yPQ*zPA*zPQ)) + v**2*xQC*yPQ*zPQ*zQD) + F3*(u**2*v*(-xPQ*yPB*zPQ**2 - xPQ*yPQ*zPA*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQD + xQC*yPQ*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xQC*yPB*zPA*zQD * F0 
      end if
    
      ! Combination 56: p_z p_z p_x p_y
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "py") then
        value = F1*(-tu*xQC*yQD + u*(xPQ*yQD*zPA*zPB + xQC*yPQ*zPA*zPB + xPQ*yQD/(2*p) + xQC*yPQ/(2*p)) + v*(-xQC*yQD*zPA*zPQ - xQC*yQD*zPB*zPQ)) + F2*(u**2*(xPQ*yPQ*zPA*zPB + xPQ*yPQ/(2*p)) + u*(-tu*xPQ*yQD - tu*xQC*yPQ + v*(-xPQ*yQD*zPA*zPQ - xPQ*yQD*zPB*zPQ - xQC*yPQ*zPA*zPQ - xQC*yPQ*zPB*zPQ)) + v**2*xQC*yQD*zPQ**2) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPQ*yPQ*zPA*zPQ - xPQ*yPQ*zPB*zPQ)) + u*v**2*(xPQ*yQD*zPQ**2 + xQC*yPQ*zPQ**2)) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + (xQC*yQD*zPA*zPB + xQC*yQD/(2*p)) * F0 
      end if
    
      ! Combination 57: p_z p_z p_x p_z
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "px" .and. o4 == "pz") then
        value = F1*(pq*(xQC*zPA + xQC*zPB) - tu*xQC*zQD + u*(xPQ*zPA*zPB*zQD + xQC*zPA*zPB*zPQ + xPQ*zQD/(2*p) + xQC*zPQ/(2*p)) + v*(-xQC*zPA*zPQ*zQD - xQC*zPB*zPQ*zQD)) + F2*(-2*pq*v*xQC*zPQ + u**2*(xPQ*zPA*zPB*zPQ + xPQ*zPQ/(2*p)) + u*(pq*(xPQ*zPA + xPQ*zPB) - tu*xPQ*zQD - tu*xQC*zPQ + v*(-xPQ*zPA*zPQ*zQD - xPQ*zPB*zPQ*zQD - xQC*zPA*zPQ**2 - xQC*zPB*zPQ**2)) + v**2*xQC*zPQ**2*zQD) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPQ*zPA*zPQ**2 - xPQ*zPB*zPQ**2)) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ*zPQ**2*zQD + xQC*zPQ**3))) + F4*u**2*v**2*xPQ*zPQ**3 + (xQC*zPA*zPB*zQD + xQC*zQD/(2*p)) * F0 
      end if
    
      ! Combination 58: p_y p_y p_y p_x
      if (o1 == "py" .and. o2 == "py" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*(xQD*yPA + xQD*yPB) - tu*xQD*yQC + u*(xPQ*yPA*yPB*yQC + xQD*yPA*yPB*yPQ + xPQ*yQC/(2*p) + xQD*yPQ/(2*p)) + v*(-xQD*yPA*yPQ*yQC - xQD*yPB*yPQ*yQC)) + F2*(-2*pq*v*xQD*yPQ + u**2*(xPQ*yPA*yPB*yPQ + xPQ*yPQ/(2*p)) + u*(pq*(xPQ*yPA + xPQ*yPB) - tu*xPQ*yQC - tu*xQD*yPQ + v*(-xPQ*yPA*yPQ*yQC - xPQ*yPB*yPQ*yQC - xQD*yPA*yPQ**2 - xQD*yPB*yPQ**2)) + v**2*xQD*yPQ**2*yQC) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPQ*yPA*yPQ**2 - xPQ*yPB*yPQ**2)) + u*(-2*pq*v*xPQ*yPQ + v**2*(xPQ*yPQ**2*yQC + xQD*yPQ**3))) + F4*u**2*v**2*xPQ*yPQ**3 + (xQD*yPA*yPB*yQC + xQD*yQC/(2*p)) * F0 
      end if
    
      ! Combination 59: p_y p_y p_z p_x
      if (o1 == "py" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(-tu*xQD*zQC + u*(xPQ*yPA*yPB*zQC + xQD*yPA*yPB*zPQ + xPQ*zQC/(2*p) + xQD*zPQ/(2*p)) + v*(-xQD*yPA*yPQ*zQC - xQD*yPB*yPQ*zQC)) + F2*(u**2*(xPQ*yPA*yPB*zPQ + xPQ*zPQ/(2*p)) + u*(-tu*xPQ*zQC - tu*xQD*zPQ + v*(-xPQ*yPA*yPQ*zQC - xPQ*yPB*yPQ*zQC - xQD*yPA*yPQ*zPQ - xQD*yPB*yPQ*zPQ)) + v**2*xQD*yPQ**2*zQC) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPQ*yPA*yPQ*zPQ - xPQ*yPB*yPQ*zPQ)) + u*v**2*(xPQ*yPQ**2*zQC + xQD*yPQ**2*zPQ)) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + (xQD*yPA*yPB*zQC + xQD*zQC/(2*p)) * F0 
      end if
    
      ! Combination 60: p_y p_z p_y p_x
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*xQD*zPB + u*(xPQ*yPA*yQC*zPB + xQD*yPA*yPQ*zPB) + v*(-xQD*yPA*yQC*zPQ - xQD*yPQ*yQC*zPB)) + F2*(-pq*v*xQD*zPQ + u**2*xPQ*yPA*yPQ*zPB + u*(pq*xPQ*zPB + v*(-xPQ*yPA*yQC*zPQ - xPQ*yPQ*yQC*zPB - xQD*yPA*yPQ*zPQ - xQD*yPQ**2*zPB)) + v**2*xQD*yPQ*yQC*zPQ) + F3*(u**2*v*(-xPQ*yPA*yPQ*zPQ - xPQ*yPQ**2*zPB) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ*yQC*zPQ + xQD*yPQ**2*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xQD*yPA*yQC*zPB * F0 
      end if
    
      ! Combination 61: p_y p_z p_z p_x
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*xQD*yPA + u*(xPQ*yPA*zPB*zQC + xQD*yPA*zPB*zPQ) + v*(-xQD*yPA*zPQ*zQC - xQD*yPQ*zPB*zQC)) + F2*(-pq*v*xQD*yPQ + u**2*xPQ*yPA*zPB*zPQ + u*(pq*xPQ*yPA + v*(-xPQ*yPA*zPQ*zQC - xPQ*yPQ*zPB*zQC - xQD*yPA*zPQ**2 - xQD*yPQ*zPB*zPQ)) + v**2*xQD*yPQ*zPQ*zQC) + F3*(u**2*v*(-xPQ*yPA*zPQ**2 - xPQ*yPQ*zPB*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQC + xQD*yPQ*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xQD*yPA*zPB*zQC * F0 
      end if
    
      ! Combination 62: p_z p_y p_y p_x
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(pq*xQD*zPA + u*(xPQ*yPB*yQC*zPA + xQD*yPB*yPQ*zPA) + v*(-xQD*yPB*yQC*zPQ - xQD*yPQ*yQC*zPA)) + F2*(-pq*v*xQD*zPQ + u**2*xPQ*yPB*yPQ*zPA + u*(pq*xPQ*zPA + v*(-xPQ*yPB*yQC*zPQ - xPQ*yPQ*yQC*zPA - xQD*yPB*yPQ*zPQ - xQD*yPQ**2*zPA)) + v**2*xQD*yPQ*yQC*zPQ) + F3*(u**2*v*(-xPQ*yPB*yPQ*zPQ - xPQ*yPQ**2*zPA) + u*(-pq*v*xPQ*zPQ + v**2*(xPQ*yPQ*yQC*zPQ + xQD*yPQ**2*zPQ))) + F4*u**2*v**2*xPQ*yPQ**2*zPQ + xQD*yPB*yQC*zPA * F0 
      end if
    
      ! Combination 63: p_z p_y p_z p_x
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*xQD*yPB + u*(xPQ*yPB*zPA*zQC + xQD*yPB*zPA*zPQ) + v*(-xQD*yPB*zPQ*zQC - xQD*yPQ*zPA*zQC)) + F2*(-pq*v*xQD*yPQ + u**2*xPQ*yPB*zPA*zPQ + u*(pq*xPQ*yPB + v*(-xPQ*yPB*zPQ*zQC - xPQ*yPQ*zPA*zQC - xQD*yPB*zPQ**2 - xQD*yPQ*zPA*zPQ)) + v**2*xQD*yPQ*zPQ*zQC) + F3*(u**2*v*(-xPQ*yPB*zPQ**2 - xPQ*yPQ*zPA*zPQ) + u*(-pq*v*xPQ*yPQ + v**2*(xPQ*yPQ*zPQ*zQC + xQD*yPQ*zPQ**2))) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + xQD*yPB*zPA*zQC * F0 
      end if
    
      ! Combination 64: p_z p_z p_y p_x
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "px") then
        value = F1*(-tu*xQD*yQC + u*(xPQ*yQC*zPA*zPB + xQD*yPQ*zPA*zPB + xPQ*yQC/(2*p) + xQD*yPQ/(2*p)) + v*(-xQD*yQC*zPA*zPQ - xQD*yQC*zPB*zPQ)) + F2*(u**2*(xPQ*yPQ*zPA*zPB + xPQ*yPQ/(2*p)) + u*(-tu*xPQ*yQC - tu*xQD*yPQ + v*(-xPQ*yQC*zPA*zPQ - xPQ*yQC*zPB*zPQ - xQD*yPQ*zPA*zPQ - xQD*yPQ*zPB*zPQ)) + v**2*xQD*yQC*zPQ**2) + F3*(u**2*(-tu*xPQ*yPQ + v*(-xPQ*yPQ*zPA*zPQ - xPQ*yPQ*zPB*zPQ)) + u*v**2*(xPQ*yQC*zPQ**2 + xQD*yPQ*zPQ**2)) + F4*u**2*v**2*xPQ*yPQ*zPQ**2 + (xQD*yQC*zPA*zPB + xQD*yQC/(2*p)) * F0 
      end if
    
      ! Combination 65: p_z p_z p_z p_x
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "px") then
        value = F1*(pq*(xQD*zPA + xQD*zPB) - tu*xQD*zQC + u*(xPQ*zPA*zPB*zQC + xQD*zPA*zPB*zPQ + xPQ*zQC/(2*p) + xQD*zPQ/(2*p)) + v*(-xQD*zPA*zPQ*zQC - xQD*zPB*zPQ*zQC)) + F2*(-2*pq*v*xQD*zPQ + u**2*(xPQ*zPA*zPB*zPQ + xPQ*zPQ/(2*p)) + u*(pq*(xPQ*zPA + xPQ*zPB) - tu*xPQ*zQC - tu*xQD*zPQ + v*(-xPQ*zPA*zPQ*zQC - xPQ*zPB*zPQ*zQC - xQD*zPA*zPQ**2 - xQD*zPB*zPQ**2)) + v**2*xQD*zPQ**2*zQC) + F3*(u**2*(-tu*xPQ*zPQ + v*(-xPQ*zPA*zPQ**2 - xPQ*zPB*zPQ**2)) + u*(-2*pq*v*xPQ*zPQ + v**2*(xPQ*zPQ**2*zQC + xQD*zPQ**3))) + F4*u**2*v**2*xPQ*zPQ**3 + (xQD*zPA*zPB*zQC + xQD*zQC/(2*p)) * F0 
      end if
    
      ! Combination 66: p_y p_y p_y p_y
      if (o1 == "py" .and. o2 == "py" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(yPA*yQC + yPA*yQD + yPB*yQC + yPB*yQD) - tu*yQC*yQD - tv*yPA*yPB + u*(yPA*yPB*yPQ*yQC + yPA*yPB*yPQ*yQD + yPQ*yQC/(2*p) + yPQ*yQD/(2*p)) + v*(-yPA*yPQ*yQC*yQD - yPB*yPQ*yQC*yQD - yPA*yPQ/(2*q) - yPB*yPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(2*pq**2 + tu*tv + u**2*(yPA*yPB*yPQ**2 + yPQ**2/(2*p)) + u*(pq*(2*yPA*yPQ + 2*yPB*yPQ) - tu*yPQ*yQC - tu*yPQ*yQD + v*(-yPA*yPQ**2*yQC - yPA*yPQ**2*yQD - yPB*yPQ**2*yQC - yPB*yPQ**2*yQD)) + v**2*(yPQ**2*yQC*yQD + yPQ**2/(2*q)) + v*(pq*(-2*yPQ*yQC - 2*yPQ*yQD) + tv*yPA*yPQ + tv*yPB*yPQ)) + F3*(-tv*v**2*yPQ**2 + u**2*(-tu*yPQ**2 + v*(-yPA*yPQ**3 - yPB*yPQ**3)) + u*(-4*pq*v*yPQ**2 + v**2*(yPQ**3*yQC + yPQ**3*yQD))) + F4*u**2*v**2*yPQ**4 + (yPA*yPB*yQC*yQD + yPA*yPB/(2*q) + yQC*yQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 67: p_y p_y p_y p_z
      if (o1 == "py" .and. o2 == "py" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*(yPA*zQD + yPB*zQD) - tu*yQC*zQD + u*(yPA*yPB*yPQ*zQD + yPA*yPB*yQC*zPQ + yPQ*zQD/(2*p) + yQC*zPQ/(2*p)) + v*(-yPA*yPQ*yQC*zQD - yPB*yPQ*yQC*zQD)) + F2*(-2*pq*v*yPQ*zQD + u**2*(yPA*yPB*yPQ*zPQ + yPQ*zPQ/(2*p)) + u*(pq*(yPA*zPQ + yPB*zPQ) - tu*yPQ*zQD - tu*yQC*zPQ + v*(-yPA*yPQ**2*zQD - yPA*yPQ*yQC*zPQ - yPB*yPQ**2*zQD - yPB*yPQ*yQC*zPQ)) + v**2*yPQ**2*yQC*zQD) + F3*(u**2*(-tu*yPQ*zPQ + v*(-yPA*yPQ**2*zPQ - yPB*yPQ**2*zPQ)) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ**3*zQD + yPQ**2*yQC*zPQ))) + F4*u**2*v**2*yPQ**3*zPQ + (yPA*yPB*yQC*zQD + yQC*zQD/(2*p)) * F0 
      end if
    
      ! Combination 68: p_y p_y p_z p_y
      if (o1 == "py" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*(yPA*zQC + yPB*zQC) - tu*yQD*zQC + u*(yPA*yPB*yPQ*zQC + yPA*yPB*yQD*zPQ + yPQ*zQC/(2*p) + yQD*zPQ/(2*p)) + v*(-yPA*yPQ*yQD*zQC - yPB*yPQ*yQD*zQC)) + F2*(-2*pq*v*yPQ*zQC + u**2*(yPA*yPB*yPQ*zPQ + yPQ*zPQ/(2*p)) + u*(pq*(yPA*zPQ + yPB*zPQ) - tu*yPQ*zQC - tu*yQD*zPQ + v*(-yPA*yPQ**2*zQC - yPA*yPQ*yQD*zPQ - yPB*yPQ**2*zQC - yPB*yPQ*yQD*zPQ)) + v**2*yPQ**2*yQD*zQC) + F3*(u**2*(-tu*yPQ*zPQ + v*(-yPA*yPQ**2*zPQ - yPB*yPQ**2*zPQ)) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ**3*zQC + yPQ**2*yQD*zPQ))) + F4*u**2*v**2*yPQ**3*zPQ + (yPA*yPB*yQD*zQC + yQD*zQC/(2*p)) * F0 
      end if
    
      ! Combination 69: p_y p_y p_z p_z
      if (o1 == "py" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(-tu*zQC*zQD - tv*yPA*yPB + u*(yPA*yPB*zPQ*zQC + yPA*yPB*zPQ*zQD + zPQ*zQC/(2*p) + zPQ*zQD/(2*p)) + v*(-yPA*yPQ*zQC*zQD - yPB*yPQ*zQC*zQD - yPA*yPQ/(2*q) - yPB*yPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(yPA*yPB*zPQ**2 + zPQ**2/(2*p)) + u*(-tu*zPQ*zQC - tu*zPQ*zQD + v*(-yPA*yPQ*zPQ*zQC - yPA*yPQ*zPQ*zQD - yPB*yPQ*zPQ*zQC - yPB*yPQ*zPQ*zQD)) + v**2*(yPQ**2*zQC*zQD + yPQ**2/(2*q)) + v*(tv*yPA*yPQ + tv*yPB*yPQ)) + F3*(u**2*(-tu*zPQ**2 + v*(-yPA*yPQ*zPQ**2 - yPB*yPQ*zPQ**2)) + v**2*(-tv*yPQ**2 + u*(yPQ**2*zPQ*zQC + yPQ**2*zPQ*zQD))) + F4*u**2*v**2*yPQ**2*zPQ**2 + (yPA*yPB*zQC*zQD + yPA*yPB/(2*q) + zQC*zQD/(2*p) + 1/(4*p*q)) * F0 
      end if
    
      ! Combination 70: p_y p_z p_y p_y
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(yQC*zPB + yQD*zPB) - tv*yPA*zPB + u*(yPA*yPQ*yQC*zPB + yPA*yPQ*yQD*zPB) + v*(-yPA*yQC*yQD*zPQ - yPQ*yQC*yQD*zPB - yPA*zPQ/(2*q) - yPQ*zPB/(2*q))) + F2*(u**2*yPA*yPQ**2*zPB + u*(2*pq*yPQ*zPB + v*(-yPA*yPQ*yQC*zPQ - yPA*yPQ*yQD*zPQ - yPQ**2*yQC*zPB - yPQ**2*yQD*zPB)) + v**2*(yPQ*yQC*yQD*zPQ + yPQ*zPQ/(2*q)) + v*(pq*(-yQC*zPQ - yQD*zPQ) + tv*yPA*zPQ + tv*yPQ*zPB)) + F3*(-tv*v**2*yPQ*zPQ + u**2*v*(-yPA*yPQ**2*zPQ - yPQ**3*zPB) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ**2*yQC*zPQ + yPQ**2*yQD*zPQ))) + F4*u**2*v**2*yPQ**3*zPQ + (yPA*yQC*yQD*zPB + yPA*zPB/(2*q)) * F0 
      end if

      ! Combination 71: p_y p_z p_y p_z
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*(yPA*yQC + zPB*zQD) + u*(yPA*yPQ*zPB*zQD + yPA*yQC*zPB*zPQ) + v*(-yPA*yQC*zPQ*zQD - yPQ*yQC*zPB*zQD)) + F2*(pq**2 + pq*v*(-yPQ*yQC - zPQ*zQD) + u**2*yPA*yPQ*zPB*zPQ + u*(pq*(yPA*yPQ + zPB*zPQ) + v*(-yPA*yPQ*zPQ*zQD - yPA*yQC*zPQ**2 - yPQ**2*zPB*zQD - yPQ*yQC*zPB*zPQ)) + v**2*yPQ*yQC*zPQ*zQD) + F3*(u**2*v*(-yPA*yPQ*zPQ**2 - yPQ**2*zPB*zPQ) + u*(pq*v*(-yPQ**2 - zPQ**2) + v**2*(yPQ**2*zPQ*zQD + yPQ*yQC*zPQ**2))) + F4*u**2*v**2*yPQ**2*zPQ**2 + yPA*yQC*zPB*zQD * F0 
      end if

      ! Combination 72: p_y p_z p_z p_y
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*(yPA*yQD + zPB*zQC) + u*(yPA*yPQ*zPB*zQC + yPA*yQD*zPB*zPQ) + v*(-yPA*yQD*zPQ*zQC - yPQ*yQD*zPB*zQC)) + F2*(pq**2 + pq*v*(-yPQ*yQD - zPQ*zQC) + u**2*yPA*yPQ*zPB*zPQ + u*(pq*(yPA*yPQ + zPB*zPQ) + v*(-yPA*yPQ*zPQ*zQC - yPA*yQD*zPQ**2 - yPQ**2*zPB*zQC - yPQ*yQD*zPB*zPQ)) + v**2*yPQ*yQD*zPQ*zQC) + F3*(u**2*v*(-yPA*yPQ*zPQ**2 - yPQ**2*zPB*zPQ) + u*(pq*v*(-yPQ**2 - zPQ**2) + v**2*(yPQ**2*zPQ*zQC + yPQ*yQD*zPQ**2))) + F4*u**2*v**2*yPQ**2*zPQ**2 + yPA*yQD*zPB*zQC * F0 
      end if

      ! Combination 73: p_y p_z p_z p_z
      if (o1 == "py" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(yPA*zQC + yPA*zQD) - tv*yPA*zPB + u*(yPA*zPB*zPQ*zQC + yPA*zPB*zPQ*zQD) + v*(-yPA*zPQ*zQC*zQD - yPQ*zPB*zQC*zQD - yPA*zPQ/(2*q) - yPQ*zPB/(2*q))) + F2*(u**2*yPA*zPB*zPQ**2 + u*(2*pq*yPA*zPQ + v*(-yPA*zPQ**2*zQC - yPA*zPQ**2*zQD - yPQ*zPB*zPQ*zQC - yPQ*zPB*zPQ*zQD)) + v**2*(yPQ*zPQ*zQC*zQD + yPQ*zPQ/(2*q)) + v*(pq*(-yPQ*zQC - yPQ*zQD) + tv*yPA*zPQ + tv*yPQ*zPB)) + F3*(-tv*v**2*yPQ*zPQ + u**2*v*(-yPA*zPQ**3 - yPQ*zPB*zPQ**2) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ*zPQ**2*zQC + yPQ*zPQ**2*zQD))) + F4*u**2*v**2*yPQ*zPQ**3 + (yPA*zPB*zQC*zQD + yPA*zPB/(2*q)) * F0 
      end if

      ! Combination 74: p_z p_y p_y p_y
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(pq*(yQC*zPA + yQD*zPA) - tv*yPB*zPA + u*(yPB*yPQ*yQC*zPA + yPB*yPQ*yQD*zPA) + v*(-yPB*yQC*yQD*zPQ - yPQ*yQC*yQD*zPA - yPB*zPQ/(2*q) - yPQ*zPA/(2*q))) + F2*(u**2*yPB*yPQ**2*zPA + u*(2*pq*yPQ*zPA + v*(-yPB*yPQ*yQC*zPQ - yPB*yPQ*yQD*zPQ - yPQ**2*yQC*zPA - yPQ**2*yQD*zPA)) + v**2*(yPQ*yQC*yQD*zPQ + yPQ*zPQ/(2*q)) + v*(pq*(-yQC*zPQ - yQD*zPQ) + tv*yPB*zPQ + tv*yPQ*zPA)) + F3*(-tv*v**2*yPQ*zPQ + u**2*v*(-yPB*yPQ**2*zPQ - yPQ**3*zPA) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ**2*yQC*zPQ + yPQ**2*yQD*zPQ))) + F4*u**2*v**2*yPQ**3*zPQ + (yPB*yQC*yQD*zPA + yPB*zPA/(2*q)) * F0 
      end if

      ! Combination 75: p_z p_y p_y p_z
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*(yPB*yQC + zPA*zQD) + u*(yPB*yPQ*zPA*zQD + yPB*yQC*zPA*zPQ) + v*(-yPB*yQC*zPQ*zQD - yPQ*yQC*zPA*zQD)) + F2*(pq**2 + pq*v*(-yPQ*yQC - zPQ*zQD) + u**2*yPB*yPQ*zPA*zPQ + u*(pq*(yPB*yPQ + zPA*zPQ) + v*(-yPB*yPQ*zPQ*zQD - yPB*yQC*zPQ**2 - yPQ**2*zPA*zQD - yPQ*yQC*zPA*zPQ)) + v**2*yPQ*yQC*zPQ*zQD) + F3*(u**2*v*(-yPB*yPQ*zPQ**2 - yPQ**2*zPA*zPQ) + u*(pq*v*(-yPQ**2 - zPQ**2) + v**2*(yPQ**2*zPQ*zQD + yPQ*yQC*zPQ**2))) + F4*u**2*v**2*yPQ**2*zPQ**2 + yPB*yQC*zPA*zQD * F0 
      end if

      ! Combination 76: p_z p_y p_z p_y
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*(yPB*yQD + zPA*zQC) + u*(yPB*yPQ*zPA*zQC + yPB*yQD*zPA*zPQ) + v*(-yPB*yQD*zPQ*zQC - yPQ*yQD*zPA*zQC)) + F2*(pq**2 + pq*v*(-yPQ*yQD - zPQ*zQC) + u**2*yPB*yPQ*zPA*zPQ + u*(pq*(yPB*yPQ + zPA*zPQ) + v*(-yPB*yPQ*zPQ*zQC - yPB*yQD*zPQ**2 - yPQ**2*zPA*zQC - yPQ*yQD*zPA*zPQ)) + v**2*yPQ*yQD*zPQ*zQC) + F3*(u**2*v*(-yPB*yPQ*zPQ**2 - yPQ**2*zPA*zPQ) + u*(pq*v*(-yPQ**2 - zPQ**2) + v**2*(yPQ**2*zPQ*zQC + yPQ*yQD*zPQ**2))) + F4*u**2*v**2*yPQ**2*zPQ**2 + yPB*yQD*zPA*zQC * F0 
      end if

      ! Combination 77: p_z p_y p_z p_z
      if (o1 == "pz" .and. o2 == "py" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(yPB*zQC + yPB*zQD) - tv*yPB*zPA + u*(yPB*zPA*zPQ*zQC + yPB*zPA*zPQ*zQD) + v*(-yPB*zPQ*zQC*zQD - yPQ*zPA*zQC*zQD - yPB*zPQ/(2*q) - yPQ*zPA/(2*q))) + F2*(u**2*yPB*zPA*zPQ**2 + u*(2*pq*yPB*zPQ + v*(-yPB*zPQ**2*zQC - yPB*zPQ**2*zQD - yPQ*zPA*zPQ*zQC - yPQ*zPA*zPQ*zQD)) + v**2*(yPQ*zPQ*zQC*zQD + yPQ*zPQ/(2*q)) + v*(pq*(-yPQ*zQC - yPQ*zQD) + tv*yPB*zPQ + tv*yPQ*zPA)) + F3*(-tv*v**2*yPQ*zPQ + u**2*v*(-yPB*zPQ**3 - yPQ*zPA*zPQ**2) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ*zPQ**2*zQC + yPQ*zPQ**2*zQD))) + F4*u**2*v**2*yPQ*zPQ**3 + (yPB*zPA*zQC*zQD + yPB*zPA/(2*q)) * F0 
      end if

      ! Combination 78: p_z p_z p_y p_y
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "py") then
        value = F1*(-tu*yQC*yQD - tv*zPA*zPB + u*(yPQ*yQC*zPA*zPB + yPQ*yQD*zPA*zPB + yPQ*yQC/(2*p) + yPQ*yQD/(2*p)) + v*(-yQC*yQD*zPA*zPQ - yQC*yQD*zPB*zPQ - zPA*zPQ/(2*q) - zPB*zPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(tu*tv + u**2*(yPQ**2*zPA*zPB + yPQ**2/(2*p)) + u*(-tu*yPQ*yQC - tu*yPQ*yQD + v*(-yPQ*yQC*zPA*zPQ - yPQ*yQC*zPB*zPQ - yPQ*yQD*zPA*zPQ - yPQ*yQD*zPB*zPQ)) + v**2*(yQC*yQD*zPQ**2 + zPQ**2/(2*q)) + v*(tv*zPA*zPQ + tv*zPB*zPQ)) + F3*(u**2*(-tu*yPQ**2 + v*(-yPQ**2*zPA*zPQ - yPQ**2*zPB*zPQ)) + v**2*(-tv*zPQ**2 + u*(yPQ*yQC*zPQ**2 + yPQ*yQD*zPQ**2))) + F4*u**2*v**2*yPQ**2*zPQ**2 + (yQC*yQD*zPA*zPB + zPA*zPB/(2*q) + yQC*yQD/(2*p) + 1/(4*p*q)) * F0 
      end if

      ! Combination 79: p_z p_z p_y p_z
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "py" .and. o4 == "pz") then
        value = F1*(pq*(yQC*zPA + yQC*zPB) - tu*yQC*zQD + u*(yPQ*zPA*zPB*zQD + yQC*zPA*zPB*zPQ + yPQ*zQD/(2*p) + yQC*zPQ/(2*p)) + v*(-yQC*zPA*zPQ*zQD - yQC*zPB*zPQ*zQD)) + F2*(-2*pq*v*yQC*zPQ + u**2*(yPQ*zPA*zPB*zPQ + yPQ*zPQ/(2*p)) + u*(pq*(yPQ*zPA + yPQ*zPB) - tu*yPQ*zQD - tu*yQC*zPQ + v*(-yPQ*zPA*zPQ*zQD - yPQ*zPB*zPQ*zQD - yQC*zPA*zPQ**2 - yQC*zPB*zPQ**2)) + v**2*yQC*zPQ**2*zQD) + F3*(u**2*(-tu*yPQ*zPQ + v*(-yPQ*zPA*zPQ**2 - yPQ*zPB*zPQ**2)) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ*zPQ**2*zQD + yQC*zPQ**3))) + F4*u**2*v**2*yPQ*zPQ**3 + (yQC*zPA*zPB*zQD + yQC*zQD/(2*p)) * F0 
      end if

      ! Combination 80: p_z p_z p_z p_y
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "py") then
        value = F1*(pq*(yQD*zPA + yQD*zPB) - tu*yQD*zQC + u*(yPQ*zPA*zPB*zQC + yQD*zPA*zPB*zPQ + yPQ*zQC/(2*p) + yQD*zPQ/(2*p)) + v*(-yQD*zPA*zPQ*zQC - yQD*zPB*zPQ*zQC)) + F2*(-2*pq*v*yQD*zPQ + u**2*(yPQ*zPA*zPB*zPQ + yPQ*zPQ/(2*p)) + u*(pq*(yPQ*zPA + yPQ*zPB) - tu*yPQ*zQC - tu*yQD*zPQ + v*(-yPQ*zPA*zPQ*zQC - yPQ*zPB*zPQ*zQC - yQD*zPA*zPQ**2 - yQD*zPB*zPQ**2)) + v**2*yQD*zPQ**2*zQC) + F3*(u**2*(-tu*yPQ*zPQ + v*(-yPQ*zPA*zPQ**2 - yPQ*zPB*zPQ**2)) + u*(-2*pq*v*yPQ*zPQ + v**2*(yPQ*zPQ**2*zQC + yQD*zPQ**3))) + F4*u**2*v**2*yPQ*zPQ**3 + (yQD*zPA*zPB*zQC + yQD*zQC/(2*p)) * F0 
      end if

      ! Combination 81: p_z p_z p_z p_z
      if (o1 == "pz" .and. o2 == "pz" .and. o3 == "pz" .and. o4 == "pz") then
        value = F1*(pq*(zPA*zQC + zPA*zQD + zPB*zQC + zPB*zQD) - tu*zQC*zQD - tv*zPA*zPB + u*(zPA*zPB*zPQ*zQC + zPA*zPB*zPQ*zQD + zPQ*zQC/(2*p) + zPQ*zQD/(2*p)) + v*(-zPA*zPQ*zQC*zQD - zPB*zPQ*zQC*zQD - zPA*zPQ/(2*q) - zPB*zPQ/(2*q)) - tu/(2*q) - tv/(2*p)) + F2*(2*pq**2 + tu*tv + u**2*(zPA*zPB*zPQ**2 + zPQ**2/(2*p)) + u*(pq*(2*zPA*zPQ + 2*zPB*zPQ) - tu*zPQ*zQC - tu*zPQ*zQD + v*(-zPA*zPQ**2*zQC - zPA*zPQ**2*zQD - zPB*zPQ**2*zQC - zPB*zPQ**2*zQD)) + v**2*(zPQ**2*zQC*zQD + zPQ**2/(2*q)) + v*(pq*(-2*zPQ*zQC - 2*zPQ*zQD) + tv*zPA*zPQ + tv*zPB*zPQ)) + F3*(-tv*v**2*zPQ**2 + u**2*(-tu*zPQ**2 + v*(-zPA*zPQ**3 - zPB*zPQ**3)) + u*(-4*pq*v*zPQ**2 + v**2*(zPQ**3*zQC + zPQ**3*zQD))) + F4*u**2*v**2*zPQ**4 + (zPA*zPB*zQC*zQD + zPA*zPB/(2*q) + zQC*zQD/(2*p) + 1/(4*p*q)) * F0 
      end if




end subroutine