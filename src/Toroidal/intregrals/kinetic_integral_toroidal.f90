subroutine kinetic_integral_ss_toroidal(r1,r2,AO1,AO2,S_ss_normal)

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
      double precision               :: const
      double precision               :: gamma_x    ,gamma_y    ,gamma_z
      double precision               :: xp  , yp  , zp
      double precision               :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision               :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
      double precision               :: I_2_gamma_x,I_2_gamma_y,I_2_gamma_z
      double precision               :: overlap_x , overlap_y , overlap_z
      double precision               :: X_k , Y_k , z_k
      double precision               :: overlap_ss

      INTERFACE
      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
      END FUNCTION gsl_sf_bessel_I0_scaled

      FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
      END FUNCTION gsl_sf_bessel_I1_scaled

      END INTERFACE

      

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

          const       = c1*c2
          const       = sign(dabs(const)**(1.0D0/3.0D0),const)

          gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
          gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))+eta
          gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))+eta

          I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
          I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
          I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))

          I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
          I_1_gamma_y = gsl_sf_bessel_I1_scaled(2.d0*gamma_y/(ay**2))
          I_1_gamma_z = gsl_sf_bessel_I1_scaled(2.d0*gamma_z/(az**2))

          I_2_gamma_x =  I_0_gamma_x - (2.d0/(gamma_x * (Lx**2/(2.d0*pi**2))))*I_1_gamma_x
          I_2_gamma_y =  I_0_gamma_y - (2.d0/(gamma_y * (Ly**2/(2.d0*pi**2))))*I_1_gamma_y
          I_2_gamma_z =  I_0_gamma_z - (2.d0/(gamma_z * (Lz**2/(2.d0*pi**2))))*I_1_gamma_z

          xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
          yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
          zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

          overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
          overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
          overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

          overlap_ss  = overlap_x * overlap_y * overlap_z 

          X_k = (beta * (cos(ax*(xp-x2)) - beta/gamma_x ) * I_1_gamma_x - 2.d0*beta**2*(1/ax*sin(ax*(xp-x2)))**2*I_2_gamma_x)  / I_0_gamma_x * overlap_ss
          Y_k = (beta * (cos(ay*(yp-y2)) - beta/gamma_y ) * I_1_gamma_y - 2.d0*beta**2*(1/ay*sin(ay*(yp-y2)))**2*I_2_gamma_y)  / I_0_gamma_y * overlap_ss
          Z_k = (beta * (cos(az*(zp-z2)) - beta/gamma_z ) * I_1_gamma_z - 2.d0*beta**2*(1/az*sin(az*(zp-z2)))**2*I_2_gamma_z)  / I_0_gamma_z * overlap_ss

          S_ss_normal =  S_ss_normal + (X_k+Y_k+Z_K)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal

subroutine kinetic_integral_ss_toroidal_11(r1,r2,AO1,AO2,S_ss_normal)

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
      double precision               :: const
      double precision               :: gamma_x    ,gamma_y    ,gamma_z
      double precision               :: xp  , yp  , zp
      double precision               :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision               :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
      double precision               :: I_2_gamma_x,I_2_gamma_y,I_2_gamma_z
      double precision               :: overlap_x , overlap_y , overlap_z
      double precision               :: X_k , Y_k , z_k
      double precision               :: overlap_ss

      INTERFACE
      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
      END FUNCTION gsl_sf_bessel_I0_scaled

      FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
      END FUNCTION gsl_sf_bessel_I1_scaled

      END INTERFACE

      

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

          !const       = c1*c2
          !const       = sign(dabs(const)**(1.0D0/3.0D0),const)

          gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
          !gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))+eta
          !gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))+eta

          I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
          !I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
          !I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))

          I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
          !I_1_gamma_y = gsl_sf_bessel_I1_scaled(2.d0*gamma_y/(ay**2))
          !I_1_gamma_z = gsl_sf_bessel_I1_scaled(2.d0*gamma_z/(az**2))

          I_2_gamma_x =  I_0_gamma_x - (2.d0/(gamma_x * (Lx**2/(2.d0*pi**2))))*I_1_gamma_x
          !I_2_gamma_y =  I_0_gamma_y - (2.d0/(gamma_y * (Ly**2/(2.d0*pi**2))))*I_1_gamma_y
          !I_2_gamma_z =  I_0_gamma_z - (2.d0/(gamma_z * (Lz**2/(2.d0*pi**2))))*I_1_gamma_z

          xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
          !yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
          !zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

          !overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
          !overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
          !overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

          !overlap_ss  = overlap_x * overlap_y * overlap_z 

          !X_k = (beta * (cos(ax*(xp-x2)) - beta/gamma_x ) * I_1_gamma_x - 2.d0*beta**2*(1/ax*sin(ax*(xp-x2)))**2*I_2_gamma_x)  / I_0_gamma_x * overlap_ss
          !Y_k = (beta * (cos(ay*(yp-y2)) - beta/gamma_y ) * I_1_gamma_y - 2.d0*beta**2*(1/ay*sin(ay*(yp-y2)))**2*I_2_gamma_y)  / I_0_gamma_y * overlap_ss
          !Z_k = (beta * (cos(az*(zp-z2)) - beta/gamma_z ) * I_1_gamma_z - 2.d0*beta**2*(1/az*sin(az*(zp-z2)))**2*I_2_gamma_z)  / I_0_gamma_z * overlap_ss

          X_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (cos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( cos(2*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) )
          Y_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (1-beta/(alpha+beta)) * I_0_gamma_x
          Z_k =  beta * (pi/(alpha+beta)) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (1-beta/(alpha+beta)) * I_0_gamma_x

          S_ss_normal =  S_ss_normal + c1 * c2 * (X_k+Y_k+Z_K)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal_11


subroutine kinetic_integral_sp_toroidal(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
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

      double precision               :: gamma_x    ,gamma_y    ,gamma_z
      double precision               :: xp         ,yp         ,zp
      double precision               :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision               :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
      double precision               :: I_2_gamma_x,I_2_gamma_y,I_2_gamma_z
      double precision               :: I_3_gamma_x 
      double precision               :: X_k , Y_k , Z_k 


      INTERFACE

        FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0
        END FUNCTION gsl_sf_bessel_I0
      
        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
        END FUNCTION gsl_sf_bessel_I0_scaled
      
        FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1
        END FUNCTION gsl_sf_bessel_I1

        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
        END FUNCTION gsl_sf_bessel_I1_scaled
      
      END INTERFACE

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
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))+eta
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))+eta

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
            
              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))

              I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
              I_1_gamma_y = gsl_sf_bessel_I1_scaled(2.d0*gamma_y/(ay**2))
              I_1_gamma_z = gsl_sf_bessel_I1_scaled(2.d0*gamma_z/(az**2))

              I_2_gamma_x = I_0_gamma_x - (2.d0/(gamma_x * (Lx**2/(2.d0*pi**2))))*I_1_gamma_x
              if (I_2_gamma_x < 1e-9)  I_2_gamma_x = 0.0d0 
              I_2_gamma_y = I_0_gamma_y - (2.d0/(gamma_y * (Ly**2/(2.d0*pi**2))))*I_1_gamma_y
              I_2_gamma_z = I_0_gamma_z - (2.d0/(gamma_z * (Lz**2/(2.d0*pi**2))))*I_1_gamma_z

              I_3_gamma_X = I_1_gamma_x - (4.0d0/(gamma_x *(Lx**2/(2.d0*pi**2))))*I_2_gamma_x 

              const       =  c1*c2
              const       =  sign(dabs(const)**(1.0d0/3.0d0),const)

              X_k = sin(ax * (xp-x2)) * I_1_gamma_x + 3.0d0 * (beta/ax**2) * sin(2.0d0 * ax * (xp-x2)) * I_2_gamma_x

              X_k = X_k - (beta/ax**2)**2 * ( 3.0d0 * sin(ax * (xp-x2)) * I_1_gamma_x - sin(3.0d0 * ax * (xp-x2)) * I_3_gamma_x)

              X_k = X_k * const * Lx   * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * ax 
                           
              X_k = X_k * const * Ly   * dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y

              X_k = X_k * const * Lz   * dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z 
              
              X_k = X_k / 2.d0 

              Y_k =       const * (beta) * Ly * dexp(-2.d0 * (alpha+beta-gamma_y) / ay**2 ) * ( cos(ay*(yp-y2)) * I_1_gamma_y + (beta/ay**2) * (dcos(2.0d0*ay*(yp-y2)) * I_2_gamma_y - I_0_gamma_y ) )
              Z_k =       const * (beta) * Lz * dexp(-2.d0 * (alpha+beta-gamma_z) / az**2 ) * ( cos(az*(zp-z2)) * I_1_gamma_z + (beta/az**2) * (dcos(2.0d0*az*(zp-z2)) * I_2_gamma_z - I_0_gamma_z ) )

              Y_k = Y_k * const * (dsin(ax*(xp-x2))/ax) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_1_gamma_x
              Z_k = Z_k * const * (dsin(ax*(xp-x2))/ax) * Lx * dexp(-2.0d0*(alpha+beta-gamma_x)/ax**2) * I_1_gamma_x

              Y_k = Y_k * const * Lz  *  dexp(-2.0d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
              Z_k = Z_k * const * Ly  *  dexp(-2.0d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y

              if (AO2%orbital=="px") S_sp_normal =  S_sp_normal + (X_k+Y_k+Z_K)
              if (AO2%orbital=="py") S_sp_normal =  0.d0 !S_sp_normal + (X_k+Y_k+Z_K)
              if (AO2%orbital=="pz") S_sp_normal =  0.d0 !S_sp_normal + (X_k+Y_k +Z_K)

        end do 
      end do


!-----------------------------------------------------------------!

end subroutine kinetic_integral_sp_toroidal


subroutine kinetic_integral_pp_toroidal(r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule

      implicit none 

      double precision  ,intent(in)  :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2
      double precision  ,intent(out) :: S_pp_normal

      integer                        :: i , j , k , l 
      double precision,parameter     :: pi = dacos(-1.d0)
      double precision               :: alpha , beta
      double precision               :: c1 , c2 
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: eta = 1e-9
      double precision               :: const , term 

      double precision               :: gamma_x    ,gamma_y    ,gamma_z
      double precision               :: xp         ,yp         ,zp
      double precision               :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision               :: I_1_gamma_x

      INTERFACE

        FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0
        END FUNCTION gsl_sf_bessel_I0
      
        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
        END FUNCTION gsl_sf_bessel_I0_scaled
      
        FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1
        END FUNCTION gsl_sf_bessel_I1

        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
        END FUNCTION gsl_sf_bessel_I1_scaled
      
      END INTERFACE

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
            
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax
              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay
              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az
            
              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
              I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
            
              const       =  c1*c2*Lx*Ly*Lz*exp(-2.d0*(alpha+beta-gamma_x)/ax**2)*exp(-2.d0*(alpha+beta-gamma_y)/ay**2)*exp(-2.d0*(alpha+beta-gamma_z)/az**2)
              term        = (dsin(ax*(xp-x2))/ax)*(dsin(ax*(xp-x1))/ax)*I_0_gamma_x + dcos(ax*(2*xp-x1-x2))*(ax**2/(2*gamma_x))*I_1_gamma_x/ax**2

              if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term

              if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term

              if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal =  0.d0 !S_pp_normal + const*I_0_gamma_y*I_0_gamma_z*term
            
        end do 
      end do 

!-----------------------------------------------------------------!

end subroutine kinetic_integral_pp_toroidal