subroutine overlap_integral_ss_tor(r1,r2,atom1,atom2,index1,index2,S_ss_normal)

      use torus_init
      use atom_basis

      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: overlap_x , overlap_y , overlap_z 
      double precision             :: gamma_x    ,gamma_y    ,gamma_z
      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision,intent(out) :: S_ss_normal


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
      
      do i = 1 , atom1%num_exponent_s
        alpha = atom1%exponent_s(i)
        c1    = atom1%coefficient_s(i,index1)
        do j = 1 , atom2%num_exponent_s
          beta = atom2%exponent_s(j)
          c2   = atom2%coefficient_s(j,index2)
            if (c1*c2 == 0.d0) cycle  
              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)

              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z  

              S_ss_normal =  S_ss_normal + overlap_x * overlap_y * overlap_z

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine

subroutine overlap_integral_sp_tor(r1,r2,atom1,atom2,index1,index2,S_sp_normal)

      use torus_init
      use atom_basis
      USE HeavisideModule
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: overlap_x  ,overlap_y  ,overlap_z 
      double precision             :: gamma_x    ,gamma_y    ,gamma_z
      double precision             :: xp         ,yp         ,zp
      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision,intent(out) :: S_sp_normal(3)


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

      S_sp_normal(:) = 0.d0

      do i = 1 , atom1%num_exponent_s
        alpha = atom1%exponent_s(i)
        c1    = atom1%coefficient_s(i,index1)
        do j = 1 , atom2%num_exponent_p
          beta = atom2%exponent_p(j)
          c2   = atom2%coefficient_p(j,index2)
            if (c1*c2 == 0.d0) cycle  
              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
        
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
        
              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
        
              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x2))/ax) * I_0_gamma_x
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2)                         * I_0_gamma_y
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2)                         * I_0_gamma_z
            
              S_sp_normal(1) =  S_sp_normal(1) + overlap_x * overlap_y * overlap_z  
              S_sp_normal(2) = 0.d0 
              S_sp_normal(3) = 0.d0 

                
        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine

subroutine overlap_integral_ps_tor(r1,r2,atom1,atom2,index1,index2,S_ps_normal)

      use torus_init
      use atom_basis
      USE HeavisideModule
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const 
      double precision             :: overlap_x  ,overlap_y  ,overlap_z 
      double precision             :: gamma_x    ,gamma_y    ,gamma_z
      double precision             :: xp         ,yp         ,zp
      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision,intent(out) :: S_ps_normal(3)


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

      S_ps_normal(:) = 0.d0

      do i = 1 , atom1%num_exponent_p
        alpha = atom1%exponent_p(i)
        c1    = atom1%coefficient_p(i,index1)
        do j = 1 , atom2%num_exponent_s
          beta = atom2%exponent_s(j)
          c2   = atom2%coefficient_s(j,index2)
            if (c1*c2 == 0.d0) cycle  
              const       = c1*c2
              const       = sign(dabs(const)**(1.0D0/3.0D0),const)
        
              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))  
              yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
              zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))
        
              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma_x/(ax**2))
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma_y/(ay**2))
              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma_z/(az**2))
        
              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * (dsin(ax*(xp-x1))/ax) * I_0_gamma_x
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2)                         * I_0_gamma_y
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2)                         * I_0_gamma_z
            
              S_ps_normal(1) =  S_ps_normal(1) + overlap_x * overlap_y * overlap_z  
              S_ps_normal(2) = 0.d0 
              S_ps_normal(3) = 0.d0 

                
        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine


subroutine overlap_integral_pp_tor(r1,r2,atom1,atom2,index1,index2,S_pp_normal)

      use torus_init
      use atom_basis
      USE HeavisideModule
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: eta = 1e-9
      double precision             :: const , term 
      double precision             :: overlap_x  ,overlap_y  ,overlap_z 
      double precision             :: gamma_x    ,gamma_y    ,gamma_z
      double precision             :: xp         ,yp         ,zp
      double precision             :: I_0_gamma_x,I_0_gamma_y,I_0_gamma_z
      double precision             :: I_1_gamma_x,I_1_gamma_y,I_1_gamma_z
      double precision,intent(out) :: S_pp_normal(3,3)


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

      S_pp_normal(:,:) = 0.d0

      do i = 1 , atom1%num_exponent_p
        alpha = atom1%exponent_p(i)
        c1    = atom1%coefficient_p(i,index1)
        do j = 1 , atom2%num_exponent_p
          beta = atom2%exponent_p(j)
          c2   = atom2%coefficient_p(j,index2)
            if (c1*c2 == 0.d0) cycle  
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
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              S_pp_normal(1,1) =  S_pp_normal(1,1) + overlap_x * overlap_y * overlap_z 
              S_pp_normal(1,2) = 0.d0  
              S_pp_normal(1,3) = 0.d0 

              term        = (dsin(ay*(yp-y2))/ay)*(dsin(ay*(yp-y1))/ay)*I_0_gamma_y + dcos(ay*(2*yp-y1-y2))*(ay**2/(2*gamma_y))*I_1_gamma_y/ay**2

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * term 
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * I_0_gamma_z

              S_pp_normal(2,1) = 0.d0
              S_pp_normal(2,2) = S_pp_normal(2,2) + overlap_x * overlap_y * overlap_z
              S_pp_normal(2,3) = 0.d0 

              term        = (dsin(az*(zp-z2))/az)*(dsin(az*(zp-z1))/az)*I_0_gamma_z + dcos(az*(2*zp-z1-z2))*(az**2/(2*gamma_z))*I_1_gamma_z/az**2

              overlap_x   = const * Lx * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * I_0_gamma_x 
              overlap_y   = const * Ly * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * I_0_gamma_y
              overlap_z   = const * Lz * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * term 


              S_pp_normal(3,1) = 0.d0
              S_pp_normal(3,2) = 0.d0  
              S_pp_normal(3,3) = S_pp_normal(3,3) + overlap_x * overlap_y * overlap_z
                
        end do 
      end do

  !-----------------------------------------------------------------!

end subroutine