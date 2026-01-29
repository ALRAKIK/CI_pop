subroutine kinetic_integral_ss_toroidal(r1,r2,AO1,AO2,S_ss_normal)

      use gsl_bessel_mod
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
      double precision               :: gamma_x
      double precision               :: xp_C 
      double precision               :: yp_R , zp_R
      double precision               :: albe , inv_albe , mu
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x
      double precision               :: I_2_gamma_x
      double precision               :: X_k , Y_k , z_k
      double precision               :: ax2
      double precision               :: const
      double precision               :: D00x , D00y , D00z
      double precision               :: S00x , S00y , S00z
      double precision                 :: theta , theta2 , cos_theta , sum_ab

      
      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      ax2         = ax * ax 

      !-----------------------------------------------------------------!

      S_ss_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

          const       = Lx*c1*c2
          
          ! Clifford Gaussian ! 

          !gamma_x     = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          call bary_exponent(alpha,beta,X,gamma_x)

          I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax2))
          I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax2))
          I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax2))

          !xp_C        = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2))
          call bary_center_toroidal(alpha,beta,x1,x2,xp_C)


          !   Real Gaussian   !

          albe        = alpha + beta 
          inv_albe    = 1.d0 / albe
          mu          = alpha * beta * inv_albe

          yp_R        = ( alpha * y1 + beta * y2 ) * inv_albe
          zp_R        = ( alpha * z1 + beta * z2 ) * inv_albe

          ! ----------------------------------------------------------- !

          D00x = -2.d0 * beta * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * ( dcos(ax*(xp_C-x2))       * I_1_gamma_x + (beta/ax2) * ( dcos(2.d0*ax*(xp_C-x2)) * I_2_gamma_x - I_0_gamma_x ) )
          D00y =                     dexp(- mu * (Y * Y))                 * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe) 
          D00z =                     dexp(- mu * (Z * Z))                 * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe) 

          S00x =      dexp(-2.d0*(alpha+beta-gamma_x)/ax2)   * I_0_gamma_x
          S00y =      dexp(- mu * (Y * Y))                   * dsqrt(pi*inv_albe)
          S00z =      dexp(- mu * (Z * Z))                   * dsqrt(pi*inv_albe)

          X_k = D00x * S00y * S00z 
          Y_k = S00x * D00y * S00z 
          Z_k = S00x * S00y * D00z 
          
          S_ss_normal =  S_ss_normal + (-0.5d0) * const  * (X_k+Y_k+Z_k)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal

subroutine kinetic_integral_sp_toroidal(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use gsl_bessel_mod
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

      double precision               :: gamma_x    
      double precision               :: xp_C
      double precision               :: I_0_gamma_x , I_1_gamma_x
      double precision               :: I_2_gamma_x , I_3_gamma_x
      double precision               :: X_k , Y_k , Z_k 
      double precision               :: albe , inv_albe
      double precision               :: mu
      double precision               :: yp_R , zp_R

      double precision               :: D01x , S01x , D01y , S01y , D01z , S01z
      double precision               :: S00x , D00x
      double precision               :: S00y , D00y
      double precision               :: S00z , D00z

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

          ! Clifford Gaussian !

          !gamma_x     = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          call bary_exponent(alpha,beta,X,gamma_x)

          !xp_C          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
          call bary_center_toroidal(alpha,beta,x1,x2,xp_C)


          I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
          I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))
          I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax**2))
          I_3_gamma_x = iv_scaled(3, 2.d0*gamma_x/(ax**2))

          !   Real Gaussian  !

          albe        = alpha + beta
          inv_albe    = 1.d0 / albe
          mu          = alpha * beta * inv_albe

          yp_R        = ( alpha * y1 + beta * y2 ) * inv_albe
          zp_R        = ( alpha * z1 + beta * z2 ) * inv_albe

          ! ----------------------------------------------------------- !

            

          if (AO2%orbital=="px") then 

          D01x = - ax*ax*dsin(ax * (xp_C-x2)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp_C-x2)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp_C-x2)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp_C-x2)) * I_3_gamma_x)
          D01x = D01x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 


          S01x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2)   * I_1_gamma_x            * (dsin(ax*(xp_C-x2))/ax) 

          D00y = dexp(- mu * (Y * Y)) * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe) 
          S00y = dexp(- mu * (Y * Y)) * dsqrt(pi*inv_albe)

          D00z = dexp(- mu * (Z * Z)) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe)
          S00z = dexp(- mu * (Z * Z)) * dsqrt(pi*inv_albe)

          X_K = D01x * S00y * S00z
          Y_K = S01x * D00y * S00z
          Z_K = S01x * S00y * D00z

          end if 

          if (AO2%orbital=="py")  then 
            
            D00x = -2.d0 * beta * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * ( dcos(ax*(xp_C-x2)) * I_1_gamma_x + (beta/ax**2) * ( dcos(2.d0*ax*(xp_C-x2)) * I_2_gamma_x - I_0_gamma_x ) )
            S00x = Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x

            D01y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y2) * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta )
            S01y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y2)

            D00z = dexp(- mu * (Z * Z)) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe)
            S00z = dexp(- mu * (Z * Z)) * dsqrt(pi*inv_albe)

            X_k  = D00x * S01y * S00z
            Y_k  = S00x * D01y * S00z
            Z_k  = S00x * S01y * D00z
            
          end if

          if (AO2%orbital=="pz") then 

            D01z = dexp(- mu * (Z * Z))   * dsqrt(pi*inv_albe) * (zp_R-z2) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta )
            S01z = dexp(- mu * (Z * Z))   * dsqrt(pi*inv_albe) * (zp_R-z2)

            X_k  = D00x * S00y * S01z
            Y_k  = S00x * D00y * S01z
            Z_k  = S00x * S00y * D01z

          end if 


          S_sp_normal =  S_sp_normal + (-0.5d0) * const * (X_k+Y_k+Z_k)

        end do 
      end do


!-----------------------------------------------------------------!

end subroutine kinetic_integral_sp_toroidal


subroutine kinetic_integral_pp_toroidal(r1,r2,AO1,AO2,S_pp_normal)

      use gsl_bessel_mod
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

      double precision               :: gamma_x    
      double precision               :: xp_C         
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x
      double precision               :: I_2_gamma_x
      double precision               :: I_3_gamma_x
      double precision               :: I_4_gamma_x
      double precision               :: X_k , Y_K , Z_K 

      double precision               :: albe , inv_albe
      double precision               :: mu
      double precision               :: yp_R , zp_R

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

          ! Clifford Gaussian !
            
          !gamma_x     = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          call bary_exponent(alpha,beta,X,gamma_x)

          !xp_C          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2)) 
          call bary_center_toroidal(alpha,beta,x1,x2,xp_C)

          I_0_gamma_x = iv_scaled(0, 2.d0*gamma_x/(ax**2))
          I_1_gamma_x = iv_scaled(1, 2.d0*gamma_x/(ax**2))
          I_2_gamma_x = iv_scaled(2, 2.d0*gamma_x/(ax**2))
          I_3_gamma_X = iv_scaled(3, 2.d0*gamma_x/(ax**2))
          I_4_gamma_X = iv_scaled(4, 2.d0*gamma_x/(ax**2))

          !   Real Gaussian  !

          albe     = alpha + beta 
          inv_albe = 1.d0 / albe
          mu       = alpha * beta * inv_albe
          yp_R     = ( alpha * y1 + beta * y2 ) * inv_albe
          zp_R     = ( alpha * z1 + beta * z2 ) * inv_albe

          ! --------------------------------------------------------- !

          S11x = ( (dsin(ax*(xp_C-x2))) * (dsin(ax*(xp_C-x1))) * I_0_gamma_x +  0.5d0 * dcos(ax*(2.d0*xp_C-x1-x2))* (I_0_gamma_x - I_2_gamma_x) ) /ax**2
          S11x = S11x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2)

          S11y = dexp(- mu * (Y*Y) ) * dsqrt(pi*inv_albe) * ( (yp_R-y1)*(yp_R-y2) + 0.5d0 * inv_albe )
          S11z = dexp(- mu * (Z*Z) ) * dsqrt(pi*inv_albe) * ( (zp_R-z1)*(zp_R-z2) + 0.5d0 * inv_albe )

          S01x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp_C-x2))/ax) * I_1_gamma_x
          S01y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y2)
          S01z = dexp(- mu * (Z*Z))   * dsqrt(pi*inv_albe) * (zp_R-z2)

          S10x = Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp_C-x1))/ax) * I_1_gamma_x
          S10y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y1)
          S10z = dexp(- mu * (Z*Z))   * dsqrt(pi*inv_albe) * (zp_R-z1)

          S00x = Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
          S00y = dexp(- mu * (Y * Y)) * dsqrt(pi*inv_albe)
          S00z = dexp(- mu * (Z * Z)) * dsqrt(pi*inv_albe)

          D11x = dcos(ax*(x2-x1))*I_0_gamma_x - dcos(ax*(2.d0*xp_C-x1-x2))*I_2_gamma_x + (3.d0*beta/ax**2) * ( dcos(ax*(xp_C-2.d0*x2+x1)) * I_1_gamma_x - dcos(ax*(3.d0*xp_C-2.d0*x2-x1)) * I_3_gamma_x ) - 0.25d0 * (1.d0 /(ax*ax)) * (2.d0*beta/ax)**2 * ( 3.d0 * dcos(ax*(x2-x1))*I_0_gamma_x - 3.d0*dcos(ax*(2.d0*xp_C-x1-x2))*I_2_gamma_x - dcos(ax*(2.d0*xp_C-3.d0*x2+x1)) * I_2_gamma_x + dcos(ax*(4.d0*xp_C-3.d0*x2-x1)) * I_4_gamma_x  )
          D11x = D11x * (- 0.5d0) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2)

          D11y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * ( (yp_R-y1)*(yp_R-y2) + 0.5d0 * inv_albe ) * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) - 4.d0 * alpha * beta * inv_albe * (yp_R-y2 * s01y  + (s11y) ) 
          D11z = dexp(- mu * (Z*Z))   * dsqrt(pi*inv_albe) * ( (zp_R-z1)*(zp_R-z2) + 0.5d0 * inv_albe ) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) - 4.d0 * alpha * beta * inv_albe * (zp_R-z2 * s01z  + (s11z) )

          D01x = - ax*ax*dsin(ax * (xp_C-x2)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp_C-x2)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp_C-x2)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp_C-x2)) * I_3_gamma_x)
          D01x = D01x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 

          D01y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y2) * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta )
          D01z = dexp(- mu * (Z*Z))   * dsqrt(pi*inv_albe) * (zp_R-z2) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta )

          D10x = - ax*ax*dsin(ax * (xp_C-x1)) * I_1_gamma_x - 3.0d0 * beta * dsin(2.0d0 * ax * (xp_C-x1)) * I_2_gamma_x + (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * dsin(ax * (xp_C-x1)) * I_1_gamma_x - dsin(3.0d0 * ax * (xp_C-x1)) * I_3_gamma_x)
          D10x = D10x * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) / ax 

          D10y = dexp(- mu * (Y*Y))   * dsqrt(pi*inv_albe) * (yp_R-y1) * ( 4.d0 * alpha * alpha * ( (yp_R-y1) * (yp_R-y1) + 0.5d0 * inv_albe )  - 2.d0 * alpha )
          D10z = dexp(- mu * (Z*Z))   * dsqrt(pi*inv_albe) * (zp_R-z1) * ( 4.d0 * alpha * alpha * ( (zp_R-z1) * (zp_R-z1) + 0.5d0 * inv_albe )  - 2.d0 * alpha )

          D00x = -2.d0 * beta * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dcos(ax*(xp_C-x2)) * I_1_gamma_x + (beta/ax**2) * ( dcos(2.d0*ax*(xp_C-x2)) * I_2_gamma_x - I_0_gamma_x ) )
          D00y = dexp(- mu * (Y * Y)) * ( 4.d0 * beta * beta * ( (yp_R-y2) * (yp_R-y2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe)
          D00z = dexp(- mu * (Z * Z)) * ( 4.d0 * beta * beta * ( (zp_R-z2) * (zp_R-z2) + 0.5d0 * inv_albe )  - 2.d0 * beta ) * dsqrt(pi*inv_albe)


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

end subroutine kinetic_integral_pp_toroidal