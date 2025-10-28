subroutine overlap_integral_ss_toroidal_3D(r1,r2,AO1,AO2,S_ss_normal)

      use gsl_bessel_mod
      use torus_init
      use classification_ERI

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

              gamma_x     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
              gamma_y     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
              gamma_z     = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

              I_0_gamma_x = bessi_scaled(0, 2.d0*gamma_x/(ax2))
              I_0_gamma_y = bessi_scaled(0, 2.d0*gamma_y/(ay2))
              I_0_gamma_z = bessi_scaled(0, 2.d0*gamma_z/(az2))


              overlap_x   = const * Lx * dexp(-2.d0*(alpha+beta-gamma_x)/ax2) * I_0_gamma_x
              overlap_y   = const * Ly * dexp(-2.d0*(alpha+beta-gamma_y)/ay2) * I_0_gamma_y
              overlap_z   = const * Lz * dexp(-2.d0*(alpha+beta-gamma_z)/az2) * I_0_gamma_z

              S_ss_normal =  S_ss_normal + overlap_x * overlap_y * overlap_z


        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_ss_toroidal_3D