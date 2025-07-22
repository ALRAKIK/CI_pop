subroutine kinetic_integral_ss_toroidal_3D(r1,r2,AO1,AO2,S_ss_normal)

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
      double precision               :: gamma_x , gamma_y , gamma_z 
      double precision               :: xp      , yp      , zp 
      double precision               :: I_0_gamma_x , I_0_gamma_y , I_0_gamma_z
      double precision               :: I_1_gamma_x , I_1_gamma_y , I_1_gamma_z
      double precision               :: I_2_gamma_x , I_2_gamma_y , I_2_gamma_z
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

          gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
          gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
          gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

          I_0_gamma_x = bessi_scaled(0, 2.d0*gamma_x/(ax**2))
          I_1_gamma_x = bessi_scaled(1, 2.d0*gamma_x/(ax**2))
          I_2_gamma_x = bessi_scaled(2, 2.d0*gamma_x/(ax**2))

          I_0_gamma_y = bessi_scaled(0, 2.d0*gamma_y/(ay**2))
          I_1_gamma_y = bessi_scaled(1, 2.d0*gamma_y/(ay**2))
          I_2_gamma_y = bessi_scaled(2, 2.d0*gamma_y/(ay**2))

          I_0_gamma_z = bessi_scaled(0, 2.d0*gamma_z/(az**2))
          I_1_gamma_z = bessi_scaled(1, 2.d0*gamma_z/(az**2))
          I_2_gamma_z = bessi_scaled(2, 2.d0*gamma_z/(az**2))

          xp          = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5 * Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
          yp          = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5 * Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
          zp          = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/ax + 0.5 * Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

          X_k =  beta * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lx * (cos(ax*(xp-x2)) * I_1_gamma_x + (beta/ax**2) * ( cos(2*ax*(xp-x2)) * I_2_gamma_x - I_0_gamma_x ) ) * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * Ly * Lz * I_0_gamma_y * I_0_gamma_z  
          Y_k =  beta * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * Ly * (cos(ay*(yp-y2)) * I_1_gamma_y + (beta/ay**2) * ( cos(2*ay*(yp-y2)) * I_2_gamma_y - I_0_gamma_y ) ) * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * Lz * Lx * I_0_gamma_z * I_0_gamma_x 
          Z_k =  beta * exp(-2.d0*(alpha+beta-gamma_z)/az**2) * Lz * (cos(az*(zp-z2)) * I_1_gamma_z + (beta/az**2) * ( cos(2*az*(zp-z2)) * I_2_gamma_z - I_0_gamma_z ) ) * exp(-2.d0*(alpha+beta-gamma_x)/ax**2) * exp(-2.d0*(alpha+beta-gamma_y)/ay**2) * Lx * Ly * I_0_gamma_x * I_0_gamma_y 

          S_ss_normal =  S_ss_normal + c1 * c2 * (X_k+Y_k+Z_k)

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_ss_toroidal_3D