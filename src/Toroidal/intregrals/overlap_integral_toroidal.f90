!subroutine overlap_integral_ss_toroidal(r1,r2,AO1,AO2,S_ss_normal)
!
!      use torus_init
!      use classification_ERI
!
!      implicit none 
!
!      !-----------------------------------------------------------------!
!
!      double precision  ,intent(in)    :: r1(3) , r2(3)
!      type(ERI_function),intent(in)    :: AO1 , AO2
!      double precision  ,intent(out)   :: S_ss_normal
!
!      integer                      :: i , j 
!      double precision,parameter   :: pi = 3.14159265358979323846D00
!      double precision             :: alpha , beta
!      double precision             :: c1    , c2 
!      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
!      double precision             :: X , Y , Z
!      double precision             :: eta = 1e-9
!      double precision             :: const 
!      double precision             :: overlap_x , overlap_y , overlap_z 
!      double precision             :: gamma_x    ,gamma_y    ,gamma_z
!      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
!
!
!      INTERFACE
!
!      FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0
!      END FUNCTION gsl_sf_bessel_I0
!
!      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
!      END FUNCTION gsl_sf_bessel_I0_scaled
!
!      END INTERFACE
!
!      x1 = r1(1) ; x2 = r2(1) 
!      y1 = r1(2) ; y2 = r2(2)
!      z1 = r1(3) ; z2 = r2(3)
!
!      X        = (x1 - x2)
!      Y        = (y1 - y2)
!      Z        = (z1 - z2)
!
!      !-----------------------------------------------------------------!
! 
!      S_ss_normal = 0.d0
!      
!      do i = 1 , size(AO1%exponent)
!        alpha = AO1%exponent(i)
!        c1    = AO1%coefficient(i)
!        do j = 1 , size(AO2%exponent)
!          beta = AO2%exponent(j)
!          c2   = AO2%coefficient(j)
!
!              const       = c1*c2
!              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
!
!              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
!              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
!              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))
!
!              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
!              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
!              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
!
!              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
!              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
!              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z  
!
!              S_ss_normal =  S_ss_normal + overlap_x * overlap_y * overlap_z
!
!        end do 
!      end do
!
!      !-----------------------------------------------------------------!
!
!end subroutine overlap_integral_ss_toroidal


subroutine overlap_integral_ss_toroidal(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_ss_normal

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: overlap_x , overlap_y , overlap_z 
      double precision             :: gamma_x    
      double precision             :: I_0_gamma_x


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

      END INTERFACE

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

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

              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta


              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
              overlap_y   = const * dsqrt(pi/(alpha+beta))
              overlap_z   = const * dsqrt(pi/(alpha+beta))

              S_ss_normal =  S_ss_normal + overlap_x * overlap_y * overlap_z

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_ss_toroidal

!subroutine overlap_integral_sp_toroidal(r1,r2,AO1,AO2,S_sp_normal)
!
!      use torus_init
!      use classification_ERI
!      use HeavisideModule
!
!      implicit none 
!
!      double precision  ,intent(in)    :: r1(3) , r2(3)
!      type(ERI_function),intent(in)    :: AO1 , AO2
!      double precision  ,intent(out)   :: S_sp_normal
!
!      integer                      :: i , j 
!      double precision,parameter   :: pi = 3.14159265358979323846D00
!      double precision             :: alpha , beta
!      double precision             :: c1    , c2 
!      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
!      double precision             :: X , Y , Z
!      double precision             :: eta = 1e-9
!      double precision             :: const 
!      double precision             :: overlap_x  ,overlap_y  ,overlap_z 
!      double precision             :: gamma_x    ,gamma_y    ,gamma_z
!      double precision             :: xp         ,yp         ,zp
!      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
!
!
!      INTERFACE
!
!      FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0
!      END FUNCTION gsl_sf_bessel_I0
!    
!      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
!      END FUNCTION gsl_sf_bessel_I0_scaled
!
!      FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I1
!      END FUNCTION gsl_sf_bessel_I1
!
!      FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
!      END FUNCTION gsl_sf_bessel_I1_scaled
!    
!      END INTERFACE
!
!      x1 = r1(1) ; x2 = r2(1) 
!      y1 = r1(2) ; y2 = r2(2)
!      z1 = r1(3) ; z2 = r2(3)
!
!      X        = (x1 - x2)
!      Y        = (y1 - y2)
!      Z        = (z1 - z2)
!
!      !-----------------------------------------------------------------!
!
!      S_sp_normal = 0.d0
!
!      do i = 1 , size(AO1%exponent)
!        alpha = AO1%exponent(i)
!        c1    = AO1%coefficient(i)
!        do j = 1 , size(AO2%exponent)
!          beta = AO2%exponent(j)
!          c2   = AO2%coefficient(j)
!
!              const       = c1*c2
!              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
!        
!              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
!              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
!              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))
!
!              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
!              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
!              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
!        
!              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
!              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
!              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
!        
!              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_0_gamma_x
!              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2)                         * I_0_gamma_y
!              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2)                         * I_0_gamma_z
!            
!              if (AO2%orbital=="px")   S_sp_normal =  S_sp_normal + overlap_x * overlap_y * overlap_z  
!              if (AO2%orbital=="py")   S_sp_normal = 0.d0 
!              if (AO2%orbital=="pz")   S_sp_normal = 0.d0 
!
!                
!        end do 
!      end do
!
!  !-----------------------------------------------------------------!
!
!end subroutine overlap_integral_sp_toroidal

subroutine overlap_integral_sp_toroidal(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule

      implicit none 

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_sp_normal

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: overlap_x  ,overlap_y  ,overlap_z 
      double precision             :: gamma_x    
      double precision             :: xp        
      double precision             :: I_0_gamma_x


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
        
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
        
              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
        
              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_0_gamma_x
              overlap_y   = const * dsqrt(pi/(alpha+beta)) 
              overlap_z   = const * dsqrt(pi/(alpha+beta)) 
            
              if (AO2%orbital=="px")   S_sp_normal =  S_sp_normal + overlap_x * overlap_y * overlap_z  
              if (AO2%orbital=="py")   S_sp_normal =  0.d0 
              if (AO2%orbital=="pz")   S_sp_normal =  0.d0 

                
        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine overlap_integral_sp_toroidal


!subroutine overlap_integral_pp_toroidal(r1,r2,AO1,AO2,S_pp_normal)
!
!      use torus_init
!      use atom_basis
!      use classification_ERI
!      use HeavisideModule
!
!      implicit none 
! 
!      double precision  ,intent(in)    :: r1(3) , r2(3)
!      type(ERI_function),intent(in)    :: AO1 , AO2
!      double precision  ,intent(out)   :: S_pp_normal
!
!      integer                          :: i , j 
!      double precision,parameter       :: pi = 3.14159265358979323846D00
!      double precision                 :: alpha , beta
!      double precision                 :: c1    , c2 
!      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 
!      double precision                 :: X , Y , Z
!      double precision                 :: eta = 1e-9
!      double precision                 :: const , term 
!      double precision                 :: overlap_x  ,overlap_y  ,overlap_z 
!      double precision                 :: gamma_x    ,gamma_y    ,gamma_z
!      double precision                 :: xp         ,yp         ,zp
!      double precision                 :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
!      double precision                 :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
!
!      INTERFACE
!
!      FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0
!      END FUNCTION gsl_sf_bessel_I0
!    
!      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
!      END FUNCTION gsl_sf_bessel_I0_scaled
!
!      FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I1
!      END FUNCTION gsl_sf_bessel_I1
!
!      FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
!        USE iso_c_binding
!        REAL(C_DOUBLE), VALUE :: x_val
!        REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
!      END FUNCTION gsl_sf_bessel_I1_scaled
!    
!      END INTERFACE
!
!      x1 = r1(1) ; x2 = r2(1) 
!      y1 = r1(2) ; y2 = r2(2)
!      z1 = r1(3) ; z2 = r2(3)
!
!      X        = (x1 - x2)
!      Y        = (y1 - y2)
!      Z        = (z1 - z2)
!
!      !-----------------------------------------------------------------!
!
!      S_pp_normal = 0.d0
!
!      do i = 1 , size(AO1%exponent)
!        alpha = AO1%exponent(i)
!        c1    = AO1%coefficient(i)
!        do j = 1 , size(AO2%exponent)
!          beta = AO2%exponent(j)
!          c2   = AO2%coefficient(j)
!
!              const       = c1*c2
!              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
!        
!              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
!              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
!              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))
!
!              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
!              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
!              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
!        
!              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
!              I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
!              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
!              I_1_gamma_y = gsl_sf_bessel_I1_scaled(2.d0*gamma_y/(ay**2))
!              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
!              I_1_gamma_z = gsl_sf_bessel_I1_scaled(2.d0*gamma_z/(az**2))
!        
!              term        = (dsin(ax*(xp-x2))/ax)*(dsin(ax*(xp-x1))/ax)*I_0_gamma_x + dcos(ax*(2*xp-x1-x2))*(ax**2/(2*gamma_x))*I_1_gamma_x/ax**2
!
!              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * term 
!              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
!              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
!
!              if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal + overlap_x * overlap_y * overlap_z 
!              if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal = 0.d0  
!              if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal = 0.d0 
!
!              term        = (dsin(ay*(yp-y2))/ay)*(dsin(ay*(yp-y1))/ay)*I_0_gamma_y + dcos(ay*(2*yp-y1-y2))*(ay**2/(2*gamma_y))*I_1_gamma_y/ay**2
!
!              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
!              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * term 
!              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z
!
!              if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal = 0.d0
!              if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z
!              if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal = 0.d0 
!
!              term        = (dsin(az*(zp-z2))/az)*(dsin(az*(zp-z1))/az)*I_0_gamma_z + dcos(az*(2*zp-z1-z2))*(az**2/(2*gamma_z))*I_1_gamma_z/az**2
!
!              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
!              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
!              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * term 
!
!
!              if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal = 0.d0
!              if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal = 0.d0  
!              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z
!                
!        end do 
!      end do
!
!  !-----------------------------------------------------------------!
!
!end subroutine overlap_integral_pp_toroidal




subroutine overlap_integral_pp_toroidal(r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use atom_basis
      use classification_ERI
      use HeavisideModule

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
      double precision                 :: eta = 1e-9
      double precision                 :: const , term 
      double precision                 :: overlap_x  ,overlap_y  ,overlap_z 
      double precision                 :: gamma_x    ,gamma_y    ,gamma_z
      double precision                 :: xp         ,yp         ,zp
      double precision                 :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision                 :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z

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
        
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
        
              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
              I_1_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_1_gamma_y = gsl_sf_bessel_I1_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
              I_1_gamma_z = gsl_sf_bessel_I1_scaled(2.d0*gamma_z/(az**2))
        
              term        = (dsin(ax*(xp-x2))/ax)*(dsin(ax*(xp-x1))/ax)*I_0_gamma_x + dcos(ax*(2*xp-x1-x2))*(ax**2/(2*gamma_x))*I_1_gamma_x/ax**2

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * term 
              overlap_y   = const * dsqrt(pi/(alpha+beta)) !Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              overlap_z   = const * dsqrt(pi/(alpha+beta)) !Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal + overlap_x * overlap_y * overlap_z 
              if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal = 0.d0  
              if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal = 0.d0 

              term        = (dsin(ay*(yp-y2))/ay)*(dsin(ay*(yp-y1))/ay)*I_0_gamma_y + dcos(ay*(2*yp-y1-y2))*(ay**2/(2*gamma_y))*I_1_gamma_y/ay**2

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
              overlap_y   = const * dsqrt(pi/(alpha+beta)) / (2.d0*(alpha+beta)) !Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * term 
              overlap_z   = const * dsqrt(pi/(alpha+beta))                       !Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal = 0.d0
              if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z
              if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal = 0.d0 

              term        = (dsin(az*(zp-z2))/az)*(dsin(az*(zp-z1))/az)*I_0_gamma_z + dcos(az*(2*zp-z1-z2))*(az**2/(2*gamma_z))*I_1_gamma_z/az**2

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
              overlap_y   = const * dsqrt(pi/(alpha+beta))                       !Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              overlap_z   = const * dsqrt(pi/(alpha+beta)) / (2.d0*(alpha+beta)) !Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * term 


              if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal = 0.d0
              if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal = 0.d0  
              if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal = S_pp_normal + overlap_x * overlap_y * overlap_z
                
        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine overlap_integral_pp_toroidal


      ! --------------------------------------------------------------- !
subroutine overlap_integral_ss_toroidal_num(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI
      use HeavisideModule

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_ss_normal

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: K 
      double precision             :: xp , yp , zp 
      double precision             :: overlap
      double precision             :: gamma_x    ,gamma_y    ,gamma_z


      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X        = (x1 - x2)
      Y        = (y1 - y2)
      Z        = (z1 - z2)

      !-----------------------------------------------------------------!
 
      S_ss_normal = 0.d0
      
      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

              k = dexp(-(alpha+beta)*1.d0/(2.d0*pi**2)*(Lx**2+Ly**2+Lz**2))

              const    = c1*c2

              gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))+eta
              gamma_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))+eta

              xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              yp       = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              zp       = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))


              call mc_integrate_overlap(gamma_x,gamma_y,gamma_z,xp,yp,zp,overlap)

              S_ss_normal =  S_ss_normal + k * const * overlap

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_ss_toroidal_num


subroutine mc_integrate_overlap(gamma_x, gamma_y, gamma_z, xp, yp, zp, result)

      use omp_lib
      use torus_init
      implicit none

      ! Input parameters

      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
      double precision, intent(in) :: xp, yp, zp
      integer, parameter           :: n_samples = 1000000                ! Number of Monte Carlo samples
    
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables

      integer                       :: i
      double precision              :: sum_f
      double precision, allocatable :: x_rand(:), y_rand(:), z_rand(:)
      double precision              :: factor_x, factor_y, factor_z
    
      ! Pre-compute constant factors

      factor_x = gamma_x * (2.0d0/ax**2)
      factor_y = gamma_y * (2.0d0/ay**2)
      factor_z = gamma_z * (2.0d0/az**2)
    
      ! Allocate arrays for random numbers

      allocate(x_rand(n_samples), y_rand(n_samples), z_rand(n_samples))
    
      ! Generate all random numbers at once

      call random_seed()
      call random_number(x_rand)
      call random_number(y_rand)
      call random_number(z_rand)
    
      ! Scale random numbers to the appropriate ranges

      x_rand = x_rand * Lx
      y_rand = y_rand * Ly
      z_rand = z_rand * Lz
    
      ! Initialize sum
      
      sum_f = 0.0d0
    
      !$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i) SCHEDULE(static)
      do i = 1, n_samples

          sum_f = sum_f + exp(factor_x * cos(ax * (x_rand(i) - xp)) + &
                              factor_y * cos(ay * (y_rand(i) - yp)) + &
                              factor_z * cos(az * (z_rand(i) - zp)))
      end do
      !$OMP END PARALLEL DO
    
      ! Calculate the final result
      result = Lx * Ly * Lz * sum_f / n_samples
      
      ! Clean up
      deallocate(x_rand, y_rand, z_rand)
    
end subroutine mc_integrate_overlap

