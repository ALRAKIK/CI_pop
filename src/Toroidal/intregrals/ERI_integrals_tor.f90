subroutine ERI_integral_4_function_toroidal(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI
      use HeavisideModule

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 

      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision,parameter      :: pi     = 3.14159265358979323846D00
      double precision,parameter      :: R2PI52 = 2.0d0*pi**(5.0d0/2.0d0)
      
      double precision                :: xa  , ya , za 
      double precision                :: xb  , yb , zb
      double precision                :: xc  , yc , zc
      double precision                :: xd  , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: zAB , zCD
      double precision                :: xp  , xq 
      double precision                :: yp  , yq
      double precision                :: zp  , zq
      double precision                :: c1  , c2 , c3 , c4 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: const  
      double precision                :: value_s 
      double precision                :: eta = 1e-9
      double precision                :: kc1 , kc2 
      double precision                :: mu_x, mu_y , mu_z
      double precision                :: nu_x, nu_y , nu_z
      double precision                :: mu  , nu 



      !-----------------------------------------------------------------!



      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital


          !kc1   = dexp(-(alpha+beta)*1.d0/(2.d0*pi**2)*(Lx**2+Ly**2+Lz**2))
          kc1   = dexp(-(alpha+beta)*(Lx**2)/(2.d0*pi**2))

          mu_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB)))+eta
          !mu_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(YAB)))+eta
          !mu_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(ZAB)))+eta

          mu_y  = alpha + beta
          mu_z  = alpha + beta

          xp    = datan((alpha*dsin(ax*xa)+beta*dsin(ax*xb))/(alpha*dcos(ax*xa)+beta*dcos(ax*xb)))/ax + 0.5*Lx * Heaviside(-alpha*dcos(ax*xa)-beta*dcos(ax*xb))  
          !yp    = datan((alpha*dsin(ay*ya)+beta*dsin(ay*yb))/(alpha*dcos(ay*ya)+beta*dcos(ay*yb)))/ay + 0.5*Ly * Heaviside(-alpha*dcos(ay*ya)-beta*dcos(ay*yb))
          !zp    = datan((alpha*dsin(az*za)+beta*dsin(az*zb))/(alpha*dcos(az*za)+beta*dcos(az*zb)))/az + 0.5*Lz * Heaviside(-alpha*dcos(az*za)-beta*dcos(az*zb))

          yp    = 0.d0 
          zp    = 0.d0 

          mu = alpha+beta 

          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital


              !kc2   = dexp(-(gamma+delta)*1.d0/(2.d0*pi**2)*(Lx**2+Ly**2+Lz**2))
              kc2   = dexp(-(gamma+delta)*(Lx**2)/(2.d0*pi**2))

              nu_x  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD)))+eta
              !nu_y  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ay*(YCD)))+eta
              !nu_z  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(az*(ZCD)))+eta

              nu_y  =  gamma+delta
              nu_z  =  gamma+delta

              xq    = datan((gamma*dsin(ax*xC)+delta*dsin(ax*xD))/(gamma*dcos(ax*xC)+delta*dcos(ax*xD)))/ax + 0.5*Lx * Heaviside(-gamma*dcos(ax*xC)-delta*dcos(ax*xD))  
              !yq    = datan((gamma*dsin(ay*yC)+delta*dsin(ay*yD))/(gamma*dcos(ay*yC)+delta*dcos(ay*yD)))/ay + 0.5*Ly * Heaviside(-gamma*dcos(ay*yC)-delta*dcos(ay*yD))
              !zq    = datan((gamma*dsin(az*zC)+delta*dsin(az*zD))/(gamma*dcos(az*zC)+delta*dcos(az*zD)))/az + 0.5*Lz * Heaviside(-gamma*dcos(az*zC)-delta*dcos(az*zD))

              yq    = 0.d0 
              zq    = 0.d0 

              !const  = (c1*c2*c3*c4) * kc1 * kc2 
              ! const  = (c1*c2*c3*c4) * kc1 * kc2 * Lx**2 * 2.d0 * pi**(3.0d0/2.0d0)

              !const  = kc1 * kc2 * Lx**2 
               

               nu    = gamma + delta

              !const  = (c1*c2*c3*c4) * kc1 * kc2
              !call grid_integrate_ERI(mu_x,mu_y,mu_z,xp,yp,zp,&
              !                     &nu_x,nu_y,nu_z,xq,yq,zq,&
              !                     &value_s)

              const   = (c1*c2*c3*c4) * 2.d0 * dsqrt(pi)*pi * Lx*Lx 
              
              !call grid_integrate_ERI_mod(mu,nu,mu_x,nu_x,value_s)
              ! call integrate_ERI_mod(mu,nu,mu_x,nu_x,value_s)

              call integrate_ERI_mod_mod(mu,nu,mu_x,nu_x,xp-xq,value_s)

              value = value + const * value_s
              
            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal

subroutine mc_integrate_ERI(gamma_x_1, gamma_y_1, gamma_z_1, xp, yp, zp,&
                          & gamma_x_2, gamma_y_2, gamma_z_2, xq, yq, zq,&
                          & result)

      use omp_lib
      use torus_init
      implicit none

      ! Input parameters

      double precision, intent(in) :: gamma_x_1, gamma_y_1, gamma_z_1
      double precision, intent(in) :: gamma_x_2, gamma_y_2, gamma_z_2
      double precision, intent(in) :: xp, yp, zp
      double precision, intent(in) :: xq, yq, zq 
      integer, parameter           :: n_samples = 1000000                ! Number of Monte Carlo samples
    
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables

      integer                       :: i
      double precision              :: sum_f
      double precision, allocatable :: x_rand_1(:), y_rand_1(:), z_rand_1(:)
      double precision, allocatable :: x_rand_2(:), y_rand_2(:), z_rand_2(:)
      double precision              :: factor_x_1, factor_y_1, factor_z_1
      double precision              :: factor_x_2, factor_y_2, factor_z_2
      double precision              :: denominator
    
      ! Pre-compute constant factors

      factor_x_1 = gamma_x_1 * (2.0d0/ax**2)
      factor_y_1 = gamma_y_1 * (2.0d0/ay**2)
      factor_z_1 = gamma_z_1 * (2.0d0/az**2)

      factor_x_2 = gamma_x_2 * (2.0d0/ax**2)
      factor_y_2 = gamma_y_2 * (2.0d0/ay**2)
      factor_z_2 = gamma_z_2 * (2.0d0/az**2)
    
      ! Allocate arrays for random numbers

      allocate(x_rand_1(n_samples), y_rand_1(n_samples), z_rand_1(n_samples))
      allocate(x_rand_2(n_samples), y_rand_2(n_samples), z_rand_2(n_samples))
    
      ! Generate all random numbers at once

      call random_seed()
      call random_number(x_rand_1)
      call random_number(y_rand_1)
      call random_number(z_rand_1)
      call random_number(x_rand_2)
      call random_number(y_rand_2)
      call random_number(z_rand_2)
    
      ! Scale random numbers to the appropriate ranges

      x_rand_1 = x_rand_1 * Lx
      y_rand_1 = y_rand_1 * Ly
      z_rand_1 = z_rand_1 * Lz

      x_rand_2 = x_rand_2 * Lx
      y_rand_2 = y_rand_2 * Ly
      z_rand_2 = z_rand_2 * Lz
    
      ! Initialize sum
      
      sum_f = 0.0d0

      !$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,denominator) COLLAPSE(1) SCHEDULE(static)
      do i = 1, n_samples

         denominator = dsqrt((1.0d0/ax**2)*(2.0d0 - 2.0d0*dcos(ax*(x_rand_1(i)-x_rand_2(i)))) + &
                            (1.0d0/ay**2)*(2.0d0 - 2.0d0*dcos(ay*(y_rand_1(i)-y_rand_2(i))))  + &
                            (1.0d0/az**2)*(2.0d0 - 2.0d0*dcos(az*(z_rand_1(i)-z_rand_2(i)))))

        if (denominator > 1.0d-9) then


          sum_f = sum_f + (dexp(factor_x_1 * dcos(ax * (x_rand_1(i) - xp))  + &
                                factor_y_1 * dcos(ay * (y_rand_1(i) - yp))  + &
                                factor_z_1 * dcos(az * (z_rand_1(i) - zp)))   &
          &             *  dexp(factor_x_2 * dcos(ax * (x_rand_2(i) - xq))  + &
                                factor_y_2 * dcos(ay * (y_rand_2(i) - yq))  + &
                                factor_z_2 * dcos(az * (z_rand_2(i) - zq))))  / denominator
        end if


      end do
      !$OMP END PARALLEL DO
    
      ! Calculate the final result
      result = Lx * Lx * Ly * Ly * Lz * Lz * sum_f / n_samples
      
      ! Clean up
      deallocate(x_rand_1, y_rand_1, z_rand_1)
      deallocate(x_rand_2, y_rand_2, z_rand_2)
    
end subroutine mc_integrate_ERI


subroutine grid_integrate_ERI(gamma_x_1, gamma_y_1, gamma_z_1, xp, yp, zp, &
                            & gamma_x_2, gamma_y_2, gamma_z_2, xq, yq, zq, &
                            & result)
      use omp_lib
      use torus_init

      implicit none

      ! Input parameters
      double precision, intent(in) :: gamma_x_1, gamma_y_1, gamma_z_1
      double precision, intent(in) :: gamma_x_2, gamma_y_2, gamma_z_2
      double precision, intent(in) :: xp, yp, zp
      double precision, intent(in) :: xq, yq, zq

      ! Grid parameters - use fewer points per dimension for 6D grid to manage memory
      integer, parameter :: nx1 = 20, ny1 = 20, nz1 = 20  ! First set grid dimensions
      integer, parameter :: nx2 = 20, ny2 = 20, nz2 = 20  ! Second set grid dimensions

      ! Output parameters
      double precision, intent(out) :: result

      ! Local variables
      integer          :: i1, j1, k1, i2, j2, k2
      double precision :: sum_f
      double precision :: x1, y1, z1, x2, y2, z2
      double precision :: dx1, dy1, dz1, dx2, dy2, dz2
      double precision :: factor_x_1, factor_y_1, factor_z_1
      double precision :: factor_x_2, factor_y_2, factor_z_2
      double precision :: x1_range, y1_range, z1_range
      double precision :: x2_range, y2_range, z2_range
      double precision :: denominator

      ! Pre-compute constant factors
      factor_x_1 = gamma_x_1 * (2.0d0/ax**2)
      factor_y_1 = gamma_y_1 
      factor_z_1 = gamma_z_1 

      factor_x_2 = gamma_x_2 * (2.0d0/ax**2)
      factor_y_2 = gamma_y_2 
      factor_z_2 = gamma_z_2 


      x1_range = Lx    
      y1_range = 5.0d0 / dsqrt(gamma_y_1) 
      z1_range = 5.0d0 / dsqrt(gamma_z_1)

      x2_range = Lx    
      y2_range = 5.0d0 / dsqrt(gamma_y_2) 
      z2_range = 5.0d0 / dsqrt(gamma_z_2)

      ! Calculate grid spacing
      dx1 = Lx / dble(nx1)
      dy1 = (2.0d0 * y1_range) / dble(ny1)
      dz1 = (2.0d0 * z1_range) / dble(nz1)

      dx2 = Lx / dble(nx2)
      dy2 = (2.0d0 * y2_range) / dble(ny2)
      dz2 = (2.0d0 * z2_range) / dble(nz2)

      ! Initialize sum
      sum_f = 0.0d0

      !$OMP PARALLEL DO REDUCTION(+:sum_f) &
      !$OMP PRIVATE(i1,j1,k1,i2,j2,k2,x1,y1,z1,x2,y2,z2,denominator) &
      !$OMP COLLAPSE(6) SCHEDULE(static)

      do k2 = 0, nz2-1
          do j2 = 0, ny2-1
              do i2 = 0, nx2-1
                  do k1 = 0, nz1-1
                      do j1 = 0, ny1-1
                          do i1 = 0, nx1-1
                              ! Calculate grid point coordinates (center of each cell)

                              x1 = i1 * dx1
                              y1 = - y1_range + j1 * dy1 
                              z1 = - z1_range + k1 * dz1 

                              x2 = i2 * dx2
                              y2 = - y2_range + j2 * dy2
                              z2 = - z2_range + k2 * dz2 

                              ! Calculate denominator
                              denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x1-x2))) + (y1-y2)**2 + (z1-z2)**2

                              if (denominator > 1.0d-10) then 
                                sum_f = sum_f + dexp((factor_x_1)*cos(ax*(x1-xp))) * dexp(-factor_y_1*y1**2) * dexp(-factor_z_1*z1**2) * dexp((factor_x_2)*cos(ax*(x2-xq))) * dexp(-factor_y_2*y2**2) * dexp(-factor_z_2*z2**2) / sqrt(denominator)
                              end if 
                              
                          end do
                      end do
                  end do
              end do
          end do
      end do
      !$OMP END PARALLEL DO

      ! Calculate the final result with the appropriate volume element
      result = dx1 * dy1 * dz1 * dx2 * dy2 * dz2 * sum_f
    
end subroutine grid_integrate_ERI


subroutine grid_integrate_ERI_mod(sigma,nu,sigma_x,nu_x,result)
    
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: sigma_x , sigma
      double precision, intent(in)  :: nu_x , nu
      

      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables
      integer                       :: i
      
      double precision              :: t
      double precision              :: dt
      double precision              :: I_0_sigma_x
      double precision              :: I_0_nu_x
      double precision              :: I_0_t_x
      double precision              :: t_range

      double precision              :: sum_f

      integer, parameter            :: nt = 10000

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
      

      ! Initialize sum

      sum_f = 0.0d0

      t_range = 10.d0 
      dt = t_range/dble(nt)


!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,t,I_0_sigma_x,I_0_nu_x,I_0_t_x) SCHEDULE(static)

      do i = 0, nt-1

        t =  i * dt

        I_0_sigma_x = gsl_sf_bessel_I0_scaled(2.d0*sigma_x/ax**2)
        !I_0_sigma_x = I_0_sigma_x * exp(2.d0*sigma_x/ax**2)

        I_0_nu_x    = gsl_sf_bessel_I0_scaled(2.d0*nu_x/ax**2)
        !I_0_nu_x    = I_0_nu_x * exp(2.d0*nu_x/ax**2)

        I_0_t_x     = gsl_sf_bessel_I0_scaled(2.d0*t**2/ax**2)
        !I_0_t_x     = I_0_t_x * exp(2.d0*t**2/ax**2)

        !sum_f = sum_f + 1.d0/((sigma*nu)+(sigma+nu)*t**2) * dexp(-2.d0*t**2/ax**2) * I_0_sigma_x * I_0_nu_x * I_0_t_x

        sum_f = sum_f + 1.d0/((sigma*nu)+(sigma+nu)*t**2) * I_0_sigma_x * I_0_nu_x * I_0_t_x  * exp(2.d0*(nu_x+sigma_x-nu-sigma)/ax**2)

      end do

!$OMP END PARALLEL DO

      ! Calculate the final result - multiply by volume element
      result =  dt *  sum_f
    
end subroutine grid_integrate_ERI_mod



subroutine integrate_ERI_mod(sigma,nu,sigma_x,nu_x,result)
      
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic

      implicit none

      
    
      ! Input parameters
      double precision, intent(in)  :: sigma_x , sigma
      double precision, intent(in)  :: nu_x , nu
      
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables

      double precision              :: epsabs, epsrel
      integer,parameter             :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer, parameter            :: limit = 100
      integer, parameter            :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      epsabs = 1.0e-8    ! Absolute error tolerance
      epsrel = 1.0e-6    ! Relative error tolerance

      
      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(x) result(fx)
        double precision, intent(in) :: x
        double precision     :: fx

        double precision :: I_0_x , I_0_sigma_x , I_0_nu_x

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

        I_0_x        = gsl_sf_bessel_I0_scaled(2.d0*x**2/ax**2)
        I_0_sigma_x  = gsl_sf_bessel_I0_scaled(2.d0*sigma_x/ax**2)
        I_0_nu_x     = gsl_sf_bessel_I0_scaled(2.d0*nu_x/ax**2)

        fx  = 1.d0/((sigma*nu)+(sigma+nu)*x**2) * I_0_x * I_0_sigma_x * I_0_nu_x * dexp(2.d0*(nu_x+sigma_x-nu-sigma)/ax**2)

      end function f_decay

end subroutine integrate_ERI_mod


subroutine integrate_ERI_mod_mod(sigma,nu,sigma_x,nu_x,xpq,result)
      
      use quadpack , only : dqagi
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod
      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)  :: sigma_x , sigma
      double precision, intent(in)  :: nu_x , nu
      double precision, intent(in)  :: xpq

      
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 100
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      integer                            :: n 
      integer,parameter                  :: Nmax = 40
      double precision,dimension(0:Nmax) :: IAB
      double precision                   :: A , B 

      INTERFACE
      FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
      END FUNCTION gsl_sf_bessel_I0_scaled
    
      END INTERFACE
      
      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)
        double precision, intent(in) :: t
        double precision             :: ft

        ft  = 1.d0/((sigma*nu)+(sigma+nu)*t**2) * S(t)

      end function f_decay

      double precision function S(t) result(sum)

      use gsl_bessel_mod
      implicit none
      double precision, intent(in) :: t
      double precision             :: A, B, C, term
      double precision             :: tol     = 1.0d-12
      integer                      :: n, Nmax = 40 

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
      
      A = 2.d0*sigma_x/(ax*ax)
      B = 2.d0*nu_x/(ax*ax)
      C = 2.d0*t*t/(ax*ax)

      sum = gsl_sf_bessel_I0_scaled(A)*gsl_sf_bessel_I0_scaled(B)*gsl_sf_bessel_I0_scaled(C) * exp(A+B-2.d0*(sigma+nu)/(ax*ax))

      do n = 1, Nmax
         term = bessi_scaled(n, A)*bessi_scaled(n, B)*bessi_scaled(n, C) * exp(A+B-2.d0*(sigma+nu)/(ax*ax))
         if (term < tol) exit
         sum = sum + term * 2.d0 * cos(dble(n)*ax*xpq) 
      end do

end function S

end subroutine integrate_ERI_mod_mod






subroutine check_openmp_enabled()
      
      use omp_lib
      
      implicit none
      
      logical :: openmp_enabled
      integer :: num_threads
  
      ! Check if OpenMP is available and enabled
      openmp_enabled = .false.
  
      !$OMP PARALLEL
      !$OMP MASTER
        openmp_enabled = .true.
      !$OMP END MASTER
      !$OMP END PARALLEL
  
      if (openmp_enabled) then
        print *, "OpenMP is enabled in this build"

        ! Check environment variable
        call omp_set_num_threads(omp_get_num_procs())
        num_threads = omp_get_max_threads()
        print *, "OMP_NUM_THREADS environment variable sets threads to:", num_threads
        print *, "Number of available processors:", omp_get_num_procs()

        ! Force a specific number to verify override works
        call omp_set_num_threads(11)

        !$OMP PARALLEL
        !$OMP MASTER
        print *, "After explicit setting, running with", omp_get_num_threads(), "threads"
        !$OMP END MASTER
        !$OMP END PARALLEL

        ! Reset to environment variable
        call omp_set_dynamic(.false.)
        num_threads = omp_get_max_threads()
        call omp_set_num_threads(num_threads)
      else
        print *, "ERROR: OpenMP is NOT enabled in this build!"
        print *, "Recompile with the appropriate OpenMP flag (-fopenmp for gfortran)"
      endif

end subroutine check_openmp_enabled