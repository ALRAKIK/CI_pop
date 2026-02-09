subroutine overlap_integral_ss_toroidal_3D(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI
      use bessel_functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_ss_normal

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha , beta
      double precision                 :: c1    , c2 
      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision                 :: X , Y , Z
      double precision                 :: const 
      double precision                 ::   overlap_x ,   overlap_y ,   overlap_z 
      double precision                 ::     gamma_x ,     gamma_y ,     gamma_z 
      double precision                 :: I_0_gamma_x , I_0_gamma_y , I_0_gamma_z
      double precision                 :: ax2 , ay2 , az2


      x1 = r1(1) ; x2 = r2(1)
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      ax2 = ax * ax 
      ay2 = ay * ay 
      az2 = az * az 

      !-----------------------------------------------------------------!
 
      S_ss_normal = 0.d0
      
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)

              !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
              !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              call bary_exponent(alpha,beta,X,gamma_x)
              call bary_exponent(alpha,beta,Y,gamma_y)
              call bary_exponent(alpha,beta,Z,gamma_z)

              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay2))
              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az2))


              overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
              overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
              overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z

              S_ss_normal =  S_ss_normal + overlap_x * overlap_y * overlap_z


        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_ss_toroidal_3D


subroutine overlap_integral_sp_toroidal_3D(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use classification_ERI
      use bessel_functions

      implicit none 

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_sp_normal

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha , beta
      double precision                 :: c1    , c2 
      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision                 :: X , Y , Z
      double precision                 :: const 
      double precision                 :: overlap_x  ,overlap_y ,overlap_z 
      double precision                 :: gamma_x    ,gamma_y   , gamma_z 
      double precision                 :: xp         , yp       , zp   
      double precision                 :: I_0_gamma_x , I_1_gamma_x
      double precision                 :: I_0_gamma_y , I_1_gamma_y
      double precision                 :: I_0_gamma_z , I_1_gamma_z


      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      !-----------------------------------------------------------------!

      S_sp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
        
              !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ax*(X)))
              !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ay*(Y)))
              !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(az*(Z)))

              call bary_exponent(alpha,beta,X,gamma_x)
              call bary_exponent(alpha,beta,Y,gamma_y)
              call bary_exponent(alpha,beta,Z,gamma_z)

              !xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5d0 * Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              !yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5d0 * Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))  
              !zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5d0 * Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))  

              call bary_center_toroidal(alpha,beta,x1,x2,xp)
              call bary_center_toroidal(alpha,beta,y1,y2,yp)
              call bary_center_toroidal(alpha,beta,z1,z2,zp)

              
              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
              I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))

              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay**2))
              I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay**2))
        
              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az**2))
              I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az**2))
        
            
              if (AO2%orbital=="px")  then 
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z  
                
              end if 

              
              if (AO2%orbital=="py")   then 
              
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              end if


              if (AO2%orbital=="pz")   then

                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

              end if 

              S_sp_normal =  S_sp_normal + overlap_x * overlap_y * overlap_z

        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine overlap_integral_sp_toroidal_3D

subroutine overlap_integral_pp_toroidal_3D(r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use atom_basis
      use classification_ERI
      use bessel_functions

      implicit none 
 
      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_pp_normal

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha , beta
      double precision                 :: c1    , c2 
      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision                 :: X , Y , Z
      double precision                 :: const , term 
      double precision                 :: overlap_x  ,overlap_y  ,overlap_z 
      double precision                 :: gamma_x    ,gamma_y    ,gamma_z
      double precision                 :: xp         ,yp         ,zp
      double precision                 :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision                 :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
      double precision                 :: I_2_gamma_x,I_2_gamma_y,I_2_gamma_z

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      !-----------------------------------------------------------------!

      S_pp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
        
              !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
              !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              call bary_exponent(alpha,beta,X,gamma_x)
              call bary_exponent(alpha,beta,Y,gamma_y)
              call bary_exponent(alpha,beta,Z,gamma_z)

              !xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              !yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              !zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

              call bary_center_toroidal(alpha,beta,x1,x2,xp)
              call bary_center_toroidal(alpha,beta,y1,y2,yp)
              call bary_center_toroidal(alpha,beta,z1,z2,zp)
        
              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
              I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))
              I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax**2))

              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay**2))
              I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay**2))
              I_2_gamma_y = iv_scaled(2, 2.d0*gamma_y/(ay**2))

              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az**2))
              I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az**2))
              I_2_gamma_z = iv_scaled(2, 2.d0*gamma_z/(az**2))

              
              if (AO1%orbital == "px" .and. AO2%orbital == "px") then 
                
                term        = ( (dsin(ax*(xp-x2))) * (dsin(ax*(xp-x1))) * I_0_gamma_x +  0.5d0 * dcos(ax*(2*xp-x1-x2))* (I_0_gamma_x - I_2_gamma_x) ) /ax**2

                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * term
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
                
                S_pp_normal =  S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 

              if (AO1%orbital == "px" .and. AO2%orbital == "py") then 
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
              
                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 


              if (AO1%orbital == "px" .and. AO2%orbital == "pz")  then
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z
                
                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "px")  then

                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) *  I_0_gamma_z
                
                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "py") then

                term        = ( (dsin(ay*(yp-y2))) * (dsin(ay*(yp-y1))) * I_0_gamma_y +  0.5d0 * dcos(ay*(2*yp-y1-y2))* (I_0_gamma_y - I_2_gamma_y) ) /ay**2

                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * term
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
                
                S_pp_normal =  S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "pz")  then
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if  


              if (AO1%orbital == "pz" .and. AO2%orbital == "px") then 
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * term
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if 


              if (AO1%orbital == "pz" .and. AO2%orbital == "py") then 
                
                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z
                
                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z 

              end if 

              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 

                term        = ( (dsin(az*(zp-z2))) * (dsin(az*(zp-z1))) * I_0_gamma_z +  0.5d0 * dcos(az*(2*zp-z1-z2))* (I_0_gamma_z - I_2_gamma_z) ) /az**2

                overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * term

                S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z

              end if
                
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_pp_toroidal_3D