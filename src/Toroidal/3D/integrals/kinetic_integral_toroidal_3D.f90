subroutine kinetic_integral_ss_toroidal_3D(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule
      use bessel_functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: S_ss_normal

      integer                        :: i , j 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: gamma_x , gamma_y , gamma_z 
      double precision               :: xp      , yp      , zp 
      double precision               :: I_0_gamma_x , I_0_gamma_y , I_0_gamma_z
      double precision               :: I_1_gamma_x , I_1_gamma_y , I_1_gamma_z
      double precision               :: I_2_gamma_x , I_2_gamma_y , I_2_gamma_z
      double precision               :: X_k , Y_k , z_k
      double precision               :: ax2 , ay2 , az2
      
      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      ax2          = ax * ax 
      ay2          = ay * ay 
      az2          = az * az 

      !-----------------------------------------------------------------!

      S_ss_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

          !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ax*(X)))
          !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ay*(Y)))
          !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(az*(Z)))

          call bary_exponent(alpha,beta,X,gamma_x)
          call bary_exponent(alpha,beta,Y,gamma_y)
          call bary_exponent(alpha,beta,Z,gamma_z)

          I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
          I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax2))
          I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax2))

          I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay2))
          I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay2))
          I_2_gamma_y = iv_scaled(2, 2.d0*gamma_y/(ay2))

          I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az2))
          I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az2))
          I_2_gamma_z = iv_scaled(2, 2.d0*gamma_z/(az2))

          !xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5 * Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2))
          !yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5 * Ly * Heaviside(-alpha*dcos(ay*y1)-beta*dcos(ay*y2))
          !zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5 * Lz * Heaviside(-alpha*dcos(az*z1)-beta*dcos(az*z2))

          call bary_center_toroidal(alpha,beta,x1,x2,xp)
          call bary_center_toroidal(alpha,beta,y1,y2,yp)
          call bary_center_toroidal(alpha,beta,z1,z2,zp)

          X_k =  beta * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * Lx * (dcos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax2) * ( dcos(2*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) ) * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * Ly * Lz * I_0_gamma_y * I_0_gamma_z  
          Y_k =  beta * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * Ly * (dcos(ay*(yp-y2)) * I_1_gamma_y + (beta/ay2) * ( dcos(2*ay*(yp-y2)) * I_2_gamma_y - I_0_gamma_y ) ) * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lz * Lx * I_0_gamma_z * I_0_gamma_x 
          Z_k =  beta * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * Lz * (dcos(az*(zp-z2)) * I_1_gamma_z + (beta/az2) * ( dcos(2*az*(zp-z2)) * I_2_gamma_z - I_0_gamma_z ) ) * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * Lx * Ly * I_0_gamma_x * I_0_gamma_y 

          S_ss_normal =  S_ss_normal + c1 * c2 * (X_k+Y_k+Z_k)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal_3D


subroutine kinetic_integral_sp_toroidal_3D(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule
      use bessel_functions

      implicit none 

      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: S_sp_normal

      integer                        :: i , j 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1 , c2 
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X  , Y  , Z
      double precision               :: const

      double precision               :: gamma_x , gamma_y , gamma_z
      double precision               :: xp , yp , zp 
      double precision               :: I_0_gamma_x , I_0_gamma_y , I_0_gamma_z
      double precision               :: I_1_gamma_x , I_1_gamma_y , I_1_gamma_z
      double precision               :: I_2_gamma_x , I_2_gamma_y , I_2_gamma_z
      double precision               :: I_3_gamma_x , I_3_gamma_y , I_3_gamma_z 
      double precision               :: X_k , Y_k , Z_k 

      double precision               :: D01x , D01y , D01z
      double precision               :: S01x , S01y , S01z
      double precision               :: D00x , D00y , D00z
      double precision               :: S00x , S00y , S00z

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!

      S_sp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

              const       =  c1*c2

              !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ax*(X)))
              !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ay*(Y)))
              !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(az*(Z)))

              call bary_exponent(alpha,beta,X,gamma_x)
              call bary_exponent(alpha,beta,Y,gamma_y)
              call bary_exponent(alpha,beta,Z,gamma_z)

              !xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5d0 * Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2))
              !yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5d0 * Ly * Heaviside(-alpha*dcos(ay*y1)-beta*dcos(ay*y2))
              !zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5d0 * Lz * Heaviside(-alpha*dcos(az*z1)-beta*dcos(az*z2))

              call bary_center_toroidal(alpha,beta,x1,x2,xp)
              call bary_center_toroidal(alpha,beta,y1,y2,yp)
              call bary_center_toroidal(alpha,beta,z1,z2,zp)
            
              I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
              I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))
              I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax**2))
              I_3_gamma_x = iv_scaled(3, 2.d0*gamma_x/(ax**2))

              I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay**2))
              I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay**2))
              I_2_gamma_y = iv_scaled(2, 2.d0*gamma_y/(ay**2))
              I_3_gamma_y = iv_scaled(3, 2.d0*gamma_y/(ay**2))

              I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az**2))
              I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az**2))
              I_2_gamma_z = iv_scaled(2, 2.d0*gamma_z/(az**2))
              I_3_gamma_z = iv_scaled(3, 2.d0*gamma_z/(az**2))

              D01x = - ax*ax*dsin(ax * (xp-x2)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp-x2)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp-x2)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp-x2)) * I_3_gamma_x)
              D01x = D01x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 

              D01y = - ay*ay*dsin(ay * (yp-y2)) * I_1_gamma_y - 3.0d0 * beta * dsin(2.0d0 * ay * (yp-y2)) * I_2_gamma_y + (2.d0*beta/ay)**2 * (0.25d0) * ( 3.0d0 * dsin(ay * (yp-y2)) * I_1_gamma_y - dsin(3.0d0 * ay * (yp-y2)) * I_3_gamma_y)
              D01y = D01y * Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) / ay

              D01z = - az*az*dsin(az * (zp-z2)) * I_1_gamma_z - 3.0d0 * beta * dsin(2.0d0 * az * (zp-z2)) * I_2_gamma_z + (2.d0*beta/az)**2 * (0.25d0) * ( 3.0d0 * dsin(az * (zp-z2)) * I_1_gamma_z - dsin(3.0d0 * az * (zp-z2)) * I_3_gamma_z)
              D01z = D01z * Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) / az

              S01x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
              S01y = Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
              S01z = Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

              S00x = Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
              S00y = Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              S00z = Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              D00x = -2.d0 * beta * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * ( dcos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( dcos(2.d0*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) )
              D00y = -2.d0 * beta * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * ( dcos(ay*(yp-y2)) * I_1_gamma_y + (beta/ay**2) * ( dcos(2.d0*ay*(yp-y2)) * I_2_gamma_y - I_0_gamma_y ) )
              D00z = -2.d0 * beta * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * ( dcos(az*(zp-z2)) * I_1_gamma_z + (beta/az**2) * ( dcos(2.d0*az*(zp-z2)) * I_2_gamma_z - I_0_gamma_z ) )


              if (AO2%orbital=="px") then 
                              
                X_k  = D01x * S00y * S00z
                Y_k  = S01x * D00y * S00z
                Z_k  = S01x * S00y * D00z
            
              end if 

              if (AO2%orbital=="py") then 
                
                X_k  = D00x * S01y * S00z
                Y_k  = S00x * D01y * S00z
                Z_k  = S00x * S01y * D00z
              
              end if 


              if (AO2%orbital=="pz") then 
                
                X_k  = D00x * S00y * S01z
                Y_k  = S00x * D00y * S01z
                Z_k  = S00x * S00y * D01z

              end if 

              S_sp_normal =  S_sp_normal + (- 0.5d0) * const * (X_k+Y_k+Z_k)

        end do 
      end do


!-----------------------------------------------------------------!

end subroutine kinetic_integral_sp_toroidal_3D


subroutine kinetic_integral_pp_toroidal_3D(r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule
      use bessel_functions

      implicit none 

      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: S_pp_normal

      integer                        :: i , j
      double precision,parameter     :: pi = dacos(-1.d0)
      double precision               :: alpha , beta
      double precision               :: c1 , c2 
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: const

      double precision               :: gamma_x, gamma_y , gamma_z    
      double precision               :: xp , yp , zp
      double precision               :: I_0_gamma_x, I_0_gamma_y, I_0_gamma_z
      double precision               :: I_1_gamma_x, I_1_gamma_y, I_1_gamma_z
      double precision               :: I_2_gamma_x, I_2_gamma_y, I_2_gamma_z
      double precision               :: I_3_gamma_x, I_3_gamma_y, I_3_gamma_z
      double precision               :: I_4_gamma_x, I_4_gamma_y, I_4_gamma_z

      double precision               :: X_k , Y_k , Z_k

      double precision               :: D11x , D11y , D11z
      double precision               :: S11x , S11y , S11z
      double precision               :: D01x , D01y , D01z
      double precision               :: S01x , S01y , S01z
      double precision               :: D10x , D10y , D10z
      double precision               :: S10x , S10y , S10z
      double precision               :: D00x , D00y , D00z
      double precision               :: S00x , S00y , S00z


      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!

      S_pp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

            const       =  c1*c2
            
            !gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ax*(X)))
            !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ay*(Y)))
            !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(az*(Z)))

            call bary_exponent(alpha,beta,X,gamma_x)
            call bary_exponent(alpha,beta,Y,gamma_y)
            call bary_exponent(alpha,beta,Z,gamma_z)

            ! xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5d0 * Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2))
            ! yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5d0 * Ly * Heaviside(-alpha*dcos(ay*y1)-beta*dcos(ay*y2))
            ! zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5d0 * Lz * Heaviside(-alpha*dcos(az*z1)-beta*dcos(az*z2))

            call bary_center_toroidal(alpha,beta,x1,x2,xp)
            call bary_center_toroidal(alpha,beta,y1,y2,yp)
            call bary_center_toroidal(alpha,beta,z1,z2,zp)

            
            I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
            I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))
            I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax**2))
            I_3_gamma_x = iv_scaled(3, 2.d0*gamma_x/(ax**2))
            I_4_gamma_x = iv_scaled(4, 2.d0*gamma_x/(ax**2))

            I_0_gamma_y = iv_scaled(0, 2.d0*gamma_y/(ay**2))
            I_1_gamma_y = iv_scaled(1, 2.d0*gamma_y/(ay**2))
            I_2_gamma_y = iv_scaled(2, 2.d0*gamma_y/(ay**2))
            I_3_gamma_y = iv_scaled(3, 2.d0*gamma_y/(ay**2))
            I_4_gamma_y = iv_scaled(4, 2.d0*gamma_y/(ay**2))

            I_0_gamma_z = iv_scaled(0, 2.d0*gamma_z/(az**2))
            I_1_gamma_z = iv_scaled(1, 2.d0*gamma_z/(az**2))
            I_2_gamma_z = iv_scaled(2, 2.d0*gamma_z/(az**2))
            I_3_gamma_z = iv_scaled(3, 2.d0*gamma_z/(az**2))
            I_4_gamma_z = iv_scaled(4, 2.d0*gamma_z/(az**2))


            D11x = dcos(ax*(x2-x1))*I_0_gamma_x - dcos(ax*(2.d0*xp-x1-x2))*I_2_gamma_x + (3.d0*beta/ax**2) * ( dcos(ax*(xp-2.d0*x2+x1)) * I_1_gamma_x - dcos(ax*(3.d0*xp-2.d0*x2-x1)) * I_3_gamma_x ) - 0.25d0 * (1.d0 /(ax*ax)) * (2.d0*beta/ax)**2 * ( 3.d0 * dcos(ax*(x2-x1))*I_0_gamma_x - 3.d0*dcos(ax*(2.d0*xp-x1-x2))*I_2_gamma_x - dcos(ax*(2.d0*xp-3.d0*x2+x1)) * I_2_gamma_x + dcos(ax*(4.d0*xp-3.d0*x2-x1)) * I_4_gamma_x  )
            D11x = D11x * (- 0.5d0) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2)

            D11y = dcos(ay*(y2-y1))*I_0_gamma_y - dcos(ay*(2.d0*yp-y1-y2))*I_2_gamma_y + (3.d0*beta/ay**2) * ( dcos(ay*(yp-2.d0*y2+y1)) * I_1_gamma_y - dcos(ay*(3.d0*yp-2.d0*y2-y1)) * I_3_gamma_y ) - 0.25d0 * (1.d0 /(ay*ay)) * (2.d0*beta/ay)**2 * ( 3.d0 * dcos(ay*(y2-y1))*I_0_gamma_y - 3.d0*dcos(ay*(2.d0*yp-y1-y2))*I_2_gamma_y - dcos(ay*(2.d0*yp-3.d0*y2+y1)) * I_2_gamma_y + dcos(ay*(4.d0*yp-3.d0*y2-y1)) * I_4_gamma_y  )
            D11y = D11y * (- 0.5d0) * Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2)

            D11z = dcos(az*(z2-z1))*I_0_gamma_z - dcos(az*(2.d0*zp-z1-z2))*I_2_gamma_z + (3.d0*beta/az**2) * ( dcos(az*(zp-2.d0*z2+z1)) * I_1_gamma_z - dcos(az*(3.d0*zp-2.d0*z2-z1)) * I_3_gamma_z ) - 0.25d0 * (1.d0 /(az*az)) * (2.d0*beta/az)**2 * ( 3.d0 * dcos(az*(z2-z1))*I_0_gamma_z - 3.d0*dcos(az*(2.d0*zp-z1-z2))*I_2_gamma_z - dcos(az*(2.d0*zp-3.d0*z2+z1)) * I_2_gamma_z + dcos(az*(4.d0*zp-3.d0*z2-z1)) * I_4_gamma_z  )
            D11z = D11z * (- 0.5d0) * Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2)

            S11x = ( (dsin(ax*(xp-x2))) * (dsin(ax*(xp-x1))) * I_0_gamma_x +  0.5d0 * dcos(ax*(2.d0*xp-x1-x2))* (I_0_gamma_x - I_2_gamma_x) ) /ax**2
            S11x = S11x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2)

            S11y = ( (dsin(ay*(yp-y2))) * (dsin(ay*(yp-y1))) * I_0_gamma_y +  0.5d0 * dcos(ay*(2.d0*yp-y1-y2))* (I_0_gamma_y - I_2_gamma_y) ) /ay**2
            S11y = S11y * Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2)

            S11z = ( (dsin(az*(zp-z2))) * (dsin(az*(zp-z1))) * I_0_gamma_z +  0.5d0 * dcos(az*(2.d0*zp-z1-z2))* (I_0_gamma_z - I_2_gamma_z) ) /az**2
            S11z = S11z * Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2)

            D01x = - ax*ax*dsin(ax * (xp-x2)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp-x2)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp-x2)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp-x2)) * I_3_gamma_x)
            D01x = D01x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 

            D01y = - ay*ay*dsin(ay * (yp-y2)) * I_1_gamma_y - 3.0d0 * beta * dsin(2.0d0 * ay * (yp-y2)) * I_2_gamma_y + (2.d0*beta/ay)**2 * (0.25d0) * ( 3.0d0 * dsin(ay * (yp-y2)) * I_1_gamma_y - dsin(3.0d0 * ay * (yp-y2)) * I_3_gamma_y)
            D01y = D01y * Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) / ay

            D01z = - az*az*dsin(az * (zp-z2)) * I_1_gamma_z - 3.0d0 * beta * dsin(2.0d0 * az * (zp-z2)) * I_2_gamma_z + (2.d0*beta/az)**2 * (0.25d0) * ( 3.0d0 * dsin(az * (zp-z2)) * I_1_gamma_z - dsin(3.0d0 * az * (zp-z2)) * I_3_gamma_z)
            D01z = D01z * Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) / az

            D10x = - ax*ax*dsin(ax * (xp-x1)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp-x1)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp-x1)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp-x1)) * I_3_gamma_x)
            D10x = D10x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 

            D10y = - ay*ay*dsin(ay * (yp-y1)) * I_1_gamma_y - 3.0d0 * beta * dsin(2.0d0 * ay * (yp-y1)) * I_2_gamma_y + (2.d0*beta/ay)**2 * (0.25d0) * ( 3.0d0 * dsin(ay * (yp-y1)) * I_1_gamma_y - dsin(3.0d0 * ay * (yp-y1)) * I_3_gamma_y)
            D10y = D10y * Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) / ay

            D10z = - az*az*dsin(az * (zp-z1)) * I_1_gamma_z - 3.0d0 * beta * dsin(2.0d0 * az * (zp-z1)) * I_2_gamma_z + (2.d0*beta/az)**2 * (0.25d0) * ( 3.0d0 * dsin(az * (zp-z1)) * I_1_gamma_z - dsin(3.0d0 * az * (zp-z1)) * I_3_gamma_z)
            D10z = D10z * Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) / az

            S01x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_1_gamma_x
            S01y = Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y2))/ay) * I_1_gamma_y
            S01z = Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z2))/az) * I_1_gamma_z

            S10x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x1))/ax) * I_1_gamma_x
            S10y = Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dsin(ay*(yp-y1))/ay) * I_1_gamma_y
            S10z = Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dsin(az*(zp-z1))/az) * I_1_gamma_z

            S00x = Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
            S00y = Ly * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
            S00z = Lz * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

            D00x = -2.d0 * beta * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dcos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( dcos(2.d0*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) )
            D00y = -2.d0 * beta * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay**2) * (dcos(ay*(yp-y2)) * I_1_gamma_y + (beta/ay**2) * ( dcos(2.d0*ay*(yp-y2)) * I_2_gamma_y - I_0_gamma_y ) )
            D00z = -2.d0 * beta * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az**2) * (dcos(az*(zp-z2)) * I_1_gamma_z + (beta/az**2) * ( dcos(2.d0*az*(zp-z2)) * I_2_gamma_z - I_0_gamma_z ) )


            if (AO1%orbital == "px" .and. AO2%orbital == "px") then 

              X_k  = D11x * S00y * S00z
              Y_k  = S11x * D00y * S00z
              Z_k  = S11x * S00y * D00z

            end if

            if (AO1%orbital == "px" .and. AO2%orbital == "py") then 

              X_k  = D10x * S01y * S00z
              Y_k  = S10x * D01y * S00z
              Z_k  = S10x * S01y * D00z
              
            end if

            if (AO1%orbital == "px" .and. AO2%orbital == "pz") then

              X_k  = D10x * S00y * S01z
              Y_k  = S10x * D00y * S01z
              Z_k  = S10x * S00y * D01z
              
            end if 

            if (AO1%orbital == "py" .and. AO2%orbital == "px") then 

              X_k  = D01x * S10y * S00z
              Y_k  = S01x * D10y * S00z
              Z_k  = S01x * S10y * D00z

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "py") then 

              X_k  = D00x * S11y * S00z
              Y_k  = S00x * D11y * S00z
              Z_k  = S00x * S11y * D00z
              
            end if 


            if (AO1%orbital == "py" .and. AO2%orbital == "pz") then 

              X_k  = D00x * S10y * S01z
              Y_k  = S00x * D10y * S01z
              Z_k  = S00x * S10y * D01z
              
            end if 

            if (AO1%orbital == "pz" .and. AO2%orbital == "px") then 

              X_k  = D01x * S00y * S10z
              Y_k  = S01x * D00y * S10z
              Z_k  = S01x * S00y * D10z
              
            end if 

            if (AO1%orbital == "pz" .and. AO2%orbital == "py") then 

              X_k  = D00x * S01y * S10z
              Y_k  = S00x * D01y * S10z
              Z_k  = S00x * S01y * D10z
              
            end if 

            if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 

              X_k  = D00x * S00y * S11z
              Y_k  = S00x * D00y * S11z
              Z_k  = S00x * S00y * D11z

            end if 

            S_pp_normal =  S_pp_normal + (- 0.5d0) * const * (X_k+Y_k+Z_k)


        end do 
      end do 

!-----------------------------------------------------------------!

end subroutine kinetic_integral_pp_toroidal_3D