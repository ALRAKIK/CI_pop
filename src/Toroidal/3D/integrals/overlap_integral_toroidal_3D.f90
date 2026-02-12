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

              call bary_exponent_x(alpha,beta,X,gamma_x)
              call bary_exponent_y(alpha,beta,Y,gamma_y)
              call bary_exponent_z(alpha,beta,Z,gamma_z)

              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay2))
              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az2))


              overlap_x   = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
              overlap_y   = Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
              overlap_z   = Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z

              S_ss_normal =  S_ss_normal + const * overlap_x * overlap_y * overlap_z


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
      double precision                 :: ax2 , ay2 , az2 


      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      ax2      = ax * ax 
      ay2      = ay * ay
      az2      = az * az

      !-----------------------------------------------------------------!

      S_sp_normal = 0.d0
      overlap_x   = 0.d0 
      overlap_y   = 0.d0 
      overlap_z   = 0.d0


      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              const       = c1*c2
      
              call bary_exponent_x(alpha,beta,X,gamma_x)
              call bary_exponent_y(alpha,beta,Y,gamma_y)
              call bary_exponent_z(alpha,beta,Z,gamma_z)

              call bary_center_toroidal_x(alpha,beta,x1,x2,xp)
              call bary_center_toroidal_y(alpha,beta,y1,y2,yp)
              call bary_center_toroidal_z(alpha,beta,z1,z2,zp)

              
              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
              I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax2))

              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay2))
              I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay2))
        
              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az2))
              I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az2))
        
            
              if (AO2%orbital=="px")  then 
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z  
                
              end if 

              
              if (AO2%orbital=="py")   then 
              
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              end if


              if (AO2%orbital=="pz")   then

                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

              end if 

              S_sp_normal =  S_sp_normal + const * overlap_x * overlap_y * overlap_z

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
      double precision                 :: ax2 , ay2 , az2

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      ax2  = ax * ax
      ay2  = ay * ay
      az2  = az * az


      !-----------------------------------------------------------------!

      S_pp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              const       = c1 * c2
        
              call bary_exponent_x(alpha,beta,X,gamma_x)
              call bary_exponent_y(alpha,beta,Y,gamma_y)
              call bary_exponent_z(alpha,beta,Z,gamma_z)

              call bary_center_toroidal_x(alpha,beta,x1,x2,xp)
              call bary_center_toroidal_y(alpha,beta,y1,y2,yp)
              call bary_center_toroidal_z(alpha,beta,z1,z2,zp)
        
              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
              I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax2))
              I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax2))

              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay2))
              I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay2))
              I_2_gamma_y = iv_scaled(2, 2.d0*gamma_y/(ay2))

              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az2))
              I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az2))
              I_2_gamma_z = iv_scaled(2, 2.d0*gamma_z/(az2))

              
              if (AO1%orbital == "px" .and. AO2%orbital == "px") then 
                
                term        = ( (dsin(ax*(xp-x2))) * (dsin(ax*(xp-x1))) * I_0_gamma_x +  0.5d0 * dcos(ax*(2.d0*xp-x1-x2))* (I_0_gamma_x - I_2_gamma_x) ) /ax2

                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * term
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z
                
              end if 

              if (AO1%orbital == "px" .and. AO2%orbital == "py") then 
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * (dsin(ax*(xp-x1))/ax) * I_1_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z
              
              end if 


              if (AO1%orbital == "px" .and. AO2%orbital == "pz")  then
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * (dsin(ax*(xp-x1))/ax) * I_1_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z
                
              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "px")  then

                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * (dsin(ay*(yp-y1))/ay) * I_1_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) *  I_0_gamma_z
                
              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "py") then

                term        = ( (dsin(ay*(yp-y2))) * (dsin(ay*(yp-y1))) * I_0_gamma_y +  0.5d0 * dcos(ay*(2.d0*yp-y1-y2)) * (I_0_gamma_y - I_2_gamma_y) ) /ay2

                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * term
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z
                
              end if 


              if (AO1%orbital == "py" .and. AO2%orbital == "pz")  then
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * (dsin(ay*(yp-y1))/ay) * I_1_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

              end if  


              if (AO1%orbital == "pz" .and. AO2%orbital == "px") then 
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * (dsin(az*(zp-z1))/az) * I_1_gamma_z

              end if 


              if (AO1%orbital == "pz" .and. AO2%orbital == "py") then 
                
                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * (dsin(az*(zp-z1))/az) * I_1_gamma_z
                
              end if 

              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 

                term        = ( (dsin(az*(zp-z2))) * (dsin(az*(zp-z1))) * I_0_gamma_z +  0.5d0 * dcos(az*(2.d0*zp-z1-z2))* (I_0_gamma_z - I_2_gamma_z) ) /az2

                overlap_x   =  Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
                overlap_y   =  Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
                overlap_z   =  Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * term


              end if

              S_pp_normal = S_pp_normal + const * overlap_x * overlap_y * overlap_z

                
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_pp_toroidal_3D