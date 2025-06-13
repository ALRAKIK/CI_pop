subroutine kinetic_integral_ss_toroidal(r1,r2,AO1,AO2,S_ss_normal)

      use gsl_bessel_mod
      use torus_init
      use classification_ERI
      use HeavisideModule

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
      double precision               :: eta = 1e-9
      double precision               :: gamma_x
      double precision               :: xp
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x
      double precision               :: I_2_gamma_x
      double precision               :: X_k , Y_k , z_k
      
      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!

      S_ss_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

          gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta

          I_0_gamma_x = bessi_scaled(0, 2.d0*gamma_x/(ax**2))
          I_1_gamma_x = bessi_scaled(1, 2.d0*gamma_x/(ax**2))
          I_2_gamma_x = bessi_scaled(2, 2.d0*gamma_x/(ax**2))

          xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))

          X_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (cos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( cos(2*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) )
          Y_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (1-beta/(alpha+beta)) * I_0_gamma_x
          Z_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (1-beta/(alpha+beta)) * I_0_gamma_x

          S_ss_normal =  S_ss_normal + c1 * c2 * (X_k+Y_k+Z_k)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal

subroutine kinetic_integral_sp_toroidal(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use gsl_bessel_mod
      use classification_ERI
      use HeavisideModule

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
      double precision               :: eta = 1e-9 , const

      double precision               :: gamma_x    
      double precision               :: xp         
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x
      double precision               :: I_2_gamma_x
      double precision               :: I_3_gamma_x 
      double precision               :: X_k , Y_k , Z_k 

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

            
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
            
              I_0_gamma_x = bessi_scaled(0, 2.d0*gamma_x/(ax**2))
              I_1_gamma_x = bessi_scaled(1, 2.d0*gamma_x/(ax**2))
              I_2_gamma_x = bessi_scaled(2, 2.d0*gamma_x/(ax**2))
              I_3_gamma_X = bessi_scaled(3, 2.d0*gamma_x/(ax**2))

              const       =  c1*c2

              X_k = ax*ax*sin(ax * (xp-x2)) * I_1_gamma_x + 3.0d0 * beta * sin(2.0d0 * ax * (xp-x2)) * I_2_gamma_x
              X_k = X_k - (2.d0*beta/ax)**2 * (0.25d0) * ( 3.0d0 * sin(ax * (xp-x2)) * I_1_gamma_x - sin(3.0d0 * ax * (xp-x2)) * I_3_gamma_x)
              X_k = X_k  * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2)
              X_k = X_k  * 0.5d0 * (pi/(alpha+beta)) / ax 

              Y_K = (alpha*beta/(alpha+beta)) * (pi/(alpha+beta)) * (dsin(ax*(xp-x2))/ax) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_1_gamma_x
              Z_K = (alpha*beta/(alpha+beta)) * (pi/(alpha+beta)) * (dsin(ax*(xp-x2))/ax) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_1_gamma_x
              

              if (AO2%orbital=="px") S_sp_normal =  S_sp_normal + const*(X_k+Y_k+Z_K)
              if (AO2%orbital=="py") S_sp_normal =  0.d0
              if (AO2%orbital=="pz") S_sp_normal =  0.d0

        end do 
      end do


!-----------------------------------------------------------------!

end subroutine kinetic_integral_sp_toroidal


subroutine kinetic_integral_pp_toroidal(r1,r2,AO1,AO2,S_pp_normal)

      use gsl_bessel_mod
      use torus_init
      use classification_ERI
      use HeavisideModule

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
      double precision               :: eta = 1e-9
      double precision               :: const

      double precision               :: gamma_x    
      double precision               :: xp         
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x
      double precision               :: I_2_gamma_x
      double precision               :: I_3_gamma_x
      double precision               :: I_4_gamma_x
      double precision               :: X_k , Y_K , Z_K 


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
            
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2)) 
            
              I_0_gamma_x = bessi_scaled(0, 2.d0*gamma_x/(ax**2))
              I_1_gamma_x = bessi_scaled(1, 2.d0*gamma_x/(ax**2))
              I_2_gamma_x = bessi_scaled(2, 2.d0*gamma_x/(ax**2))
              I_3_gamma_X = bessi_scaled(3, 2.d0*gamma_x/(ax**2))
              I_4_gamma_X = bessi_scaled(4, 2.d0*gamma_x/(ax**2))


              X_K = 3.d0*cos(ax*(x2-x1))*I_0_gamma_x - 3.d0*cos(ax*(2.d0*xp-x1-x2))*I_2_gamma_x
              X_K = X_K - cos(ax*(2.d0*xp-3.d0*x2-x1)) * I_2_gamma_x + cos(ax*(4.d0*xp-3.d0*x2-x1)) * I_4_gamma_x
              X_K = X_K * 0.25d0 * (2.d0*beta/ax)**2 * (-1.d0) /(ax*ax)   
              X_K = X_K + (3.d0*beta/ax**2) * ( cos(ax*(xp-2.d0*x2+x1)) * I_1_gamma_x - cos(ax*(3.d0*xp-2.d0*x2-x1)) * I_3_gamma_x )
              X_K = X_K + cos(ax*(x2-x1)) * I_0_gamma_x - cos(ax*(2*xp-x1-x2)) * I_2_gamma_x
              X_k = X_k * 0.25d0
              X_K = X_K * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (pi/(alpha+beta))
               

              Y_K = (dsin(ax*(xp-x2)))*(dsin(ax*(xp-x1)))*I_0_gamma_x + dcos(ax*(2*xp-x1-x2))*(ax**2/(2*gamma_x))*I_1_gamma_x
              Y_K = Y_K * ((alpha*beta)/(alpha+beta)) /ax/ax 
              Y_K = Y_K * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (pi/(alpha+beta))
              
              Z_K = (dsin(ax*(xp-x2)))*(dsin(ax*(xp-x1)))*I_0_gamma_x + dcos(ax*(2*xp-x1-x2))*(ax**2/(2*gamma_x))*I_1_gamma_x
              Z_K = Z_K * ((alpha*beta)/(alpha+beta)) /ax/ax 
              Z_K = Z_K * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (pi/(alpha+beta))

              if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal + const*(X_k+Y_k+Z_K) 
              if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal =  0.d0
              if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal =  0.d0

              X_K = 0.d0 
              Y_K = 0.d0 
              Z_K = 0.d0 

              X_K = 0.5d0 * beta                      * (pi/(alpha+beta)) * (1.d0/(alpha+beta)) * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (cos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( cos(2*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) )
              Y_K = 1.5d0 * (alpha*beta/(alpha+beta)) * (pi/(alpha+beta)) * (1.d0/(alpha+beta)) * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
              Z_K = 0.5d0 * (alpha*beta/(alpha+beta)) * (pi/(alpha+beta)) * (1.d0/(alpha+beta)) * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 


              if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal =  0.d0
              if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal + const*(X_k+Y_k+Z_K)
              if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal =  0.d0

              if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal =  0.d0
              if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal =  0.d0
              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal + const*(X_k+Y_k+Z_K)
            
        end do 
      end do 

!-----------------------------------------------------------------!

end subroutine kinetic_integral_pp_toroidal