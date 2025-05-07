subroutine nuclear_attraction_integral_ss_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use atom_basis
      use classification_ERI
      use HeavisideModule
      use omp_lib

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)  :: AO1 , AO2
      type(atom)        ,intent(in)  :: atoms(number_of_atoms)
      integer           ,intent(in)  :: number_of_atoms
      double precision  ,intent(in)  :: r1(3) , r2(3)
      double precision  ,intent(in)  :: geometry(number_of_atoms,3)
      double precision  ,intent(out) :: S_ss_normal

      integer                        :: i , j , k 
      integer                        :: charge_atom
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: kc
      double precision               :: eta = 1e-9
      double precision               :: x1 , x2 , X
      double precision               :: y1 , y2 , Y 
      double precision               :: z1 , z2 , Z
      double precision               :: xc , yc , zc 
      double precision               :: xp , yp , zp 
      double precision               :: gamma_x , gamma_y , gamma_z
      double precision               :: NA


      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X  = (x1 - x2)
      Y  = (y1 - y2)
      Z  = (z1 - z2)

      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      
      S_ss_normal = 0.d0
      NA          = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)
            
            !kc       = dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2))
            !kc       = 2.d0*sqrt(pi)*dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2)) * Lx 

            gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
            gamma_y  = alpha+beta
            gamma_z  = alpha+beta

            xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            yp       = 0.d0
            zp       = 0.d0

            kc       = 2.d0*sqrt(pi) * Lx

            do k = 1 , number_of_atoms

              xc = geometry(k,1) 
              yc = geometry(k,2)
              zc = geometry(k,3)
              
              charge_atom = (-1)*atoms(k)%charge

              !call grid_integrate_NA_ss(gamma_x,gamma_y,gamma_z,xp,yp,zp,xc,yc,zc,NA)
              
              !call grid_integrate_NA_ss_mod(gamma_x,xp-xc,gamma_y, NA)

               call integrate_NA_ss_mod(gamma_x,xp-xc,gamma_y, NA)
              S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * kc * NA
            
            end do 
            
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_ss_toroidal


subroutine nuclear_attraction_integral_sp_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use atom_basis
      use classification_ERI
      use HeavisideModule
      use omp_lib

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)  :: AO1 , AO2
      type(atom)        ,intent(in)  :: atoms(number_of_atoms)
      integer           ,intent(in)  :: number_of_atoms
      double precision  ,intent(in)  :: r1(3) , r2(3)
      double precision  ,intent(in)  :: geometry(number_of_atoms,3)
      double precision  ,intent(out) :: S_sp_normal


      integer                       :: i , j , k 
      integer                       :: charge_atom
      double precision,parameter    :: pi = 3.14159265358979323846D00
      double precision              :: alpha , beta
      double precision              :: c1    , c2 
      double precision              :: kc
      double precision              :: eta = 1e-9
      double precision              :: xc , yc , zc
      double precision              :: x1 , x2 , X
      double precision              :: y1 , y2 , Y 
      double precision              :: z1 , z2 , Z 
      double precision              :: xp , yp , zp 
      double precision              :: gamma_x , gamma_y , gamma_z
      double precision              :: NA
      

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      
      S_sp_normal = 0.d0
      NA          = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)
            
            !kc       = dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2)) 

            gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
            gamma_y  = alpha+beta
            gamma_z  = alpha+beta

            xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            yp       = 0.d0
            zp       = 0.d0

            kc       = dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2)) * c1 * c2  * 2.d0 * dsqrt(pi) / ax * Lx 

            if (AO2%orbital=="px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                !call grid_integrate_NA_sp(gamma_x,gamma_y,gamma_z,xp,x1,x2,xc,yc,zc,NA)
                
                call integrate_NA_sp_mod(gamma_x,xp,xc,x2,alpha+beta,NA)
              
                !S_sp_normal =  S_sp_normal +  c1 * c2 * charge_atom * kc * NA / ax

                S_sp_normal =  S_sp_normal + charge_atom * kc * NA
              
              end do

            end if 

        end do 
      end do



      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_sp_toroidal

subroutine nuclear_attraction_integral_pp_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use atom_basis
      use classification_ERI
      use HeavisideModule
      use omp_lib

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)  :: AO1 , AO2
      type(atom)        ,intent(in)  :: atoms(number_of_atoms)
      integer           ,intent(in)  :: number_of_atoms
      double precision  ,intent(in)  :: r1(3) , r2(3)
      double precision  ,intent(in)  :: geometry(number_of_atoms,3)
      double precision  ,intent(out) :: S_pp_normal


      integer                       :: i , j , k 
      integer                       :: charge_atom
      double precision,parameter    :: pi = 3.14159265358979323846D00
      double precision              :: alpha , beta
      double precision              :: c1    , c2 
      double precision              :: kc
      double precision              :: eta = 1e-9
      double precision              :: xc , yc , zc
      double precision              :: x1 , x2 , X
      double precision              :: y1 , y2 , Y 
      double precision              :: z1 , z2 , Z 
      double precision              :: xp , yp , zp 
      double precision              :: gamma_x , gamma_y , gamma_z
      double precision              :: NA
      

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      
      S_pp_normal = 0.d0
      NA          = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)
            
            !kc       = dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2))

            gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
            gamma_y  = alpha+beta
            gamma_z  = alpha+beta

            xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            yp       = 0.d0
            zp       = 0.d0


            kc       = c1 * c2 * dsqrt(pi) * Lx * dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2)) / (ax**2) 

            if (AO1%orbital == "px" .and. AO2%orbital == "px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                !call grid_integrate_NA_pp_px(gamma_x,gamma_y,gamma_z,xp,x1,x2,xc,yc,zc,NA)
                call integrate_NA_pp_px_mod(gamma_x,xp,alpha+beta,xc,x1,x2,NA)
              
                !S_pp_normal =  S_pp_normal +  c1 * c2 * charge_atom * kc * NA / (ax**2)
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA 

              
              end do

            end if 

            kc       = c1 * c2 * dsqrt(pi) * Lx 


            if (AO1%orbital == "py" .and. AO2%orbital == "py") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                !call grid_integrate_NA_pp_py(gamma_x,gamma_y,gamma_z,xp,x1,x2,xc,yc,zc,NA)
                call integrate_NA_pp_py_mod(gamma_x,xp-xc,alpha+beta, NA)
              
                
                !S_pp_normal =  S_pp_normal +  c1 * c2 * charge_atom * kc * NA
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA
              
              end do

            end if

            kc       = c1 * c2 * dsqrt(pi) * Lx
            
            if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                 
                 charge_atom = (-1)*atoms(k)%charge
   
                !call grid_integrate_NA_pp_pz(gamma_x,gamma_y,gamma_z,xp,x1,x2,xc,yc,zc,NA)
                call integrate_NA_pp_py_mod(gamma_x,xp-xc,alpha+beta, NA)
              
                !S_pp_normal =  S_pp_normal +  c1 * c2 * charge_atom * kc * NA
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA
              
              end do
!
            end if

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_pp_toroidal


subroutine grid_integrate_NA_ss(gamma_x, gamma_y, gamma_z, xp, yp, zp, xc, yc, zc, result)
    
      use omp_lib
      use torus_init
    
      implicit none
    
      ! Input parameters
      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
      double precision, intent(in) :: xp, yp, zp
      double precision, intent(in) :: xc, yc, zc 

      ! Grid parameters
      integer, parameter :: nx = 200      ! Number of grid points in x direction
      integer, parameter :: ny = 100      ! Number of grid points in y direction
      integer, parameter :: nz = 100      ! Number of grid points in z direction

      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables
      integer          :: i, j, k
      double precision :: sum_f
      double precision :: x, y, z
      double precision :: dx, dy, dz
      double precision :: x_range, y_range, z_range
      double precision :: factor_x, factor_y, factor_z
      double precision :: denominator

      ! Pre-compute constant factors
      factor_x = gamma_x * (2.0d0/ax**2)
      factor_y = gamma_y 
      factor_z = gamma_z

      ! Calculate grid spacing

      x_range = Lx    
      y_range = 5.0d0 / dsqrt(gamma_y) 
      z_range = 5.0d0 / dsqrt(gamma_z) 
      
      dx = Lx   / dble(nx)
      dy = (2.0d0 * y_range) / dble(ny)
      dz = (2.0d0 * z_range) / dble(nz)

      ! Initialize sum

      sum_f = 0.0d0


!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,j,k,x,y,z,denominator) COLLAPSE(3) SCHEDULE(static)
      do k = 0, nz-1
          do j = 0, ny-1
              do i = 0, nx-1

                  x =  i * dx
                  y =  - y_range + j * dy 
                  z =  - z_range + k * dz

                  ! Calculate denominator

                  denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x-xc))) + (y-yc)**2 + (z-zc)**2

                  if (denominator > 1.0d-10) then
                    sum_f = sum_f + dexp((factor_x)*cos(ax*(x-xp))) * dexp(-factor_y*y**2) * dexp(-factor_z*z**2) / sqrt(denominator)
                  end if

              end do
          end do
      end do
!$OMP END PARALLEL DO

      ! Calculate the final result - multiply by volume element
      result = dx * dy * dz * sum_f
    
end subroutine grid_integrate_NA_ss



subroutine grid_integrate_NA_ss_mod(gamma_x,x,p, result)
    
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , X , p 
      

      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables
      integer                       :: i
      double precision              :: sum_f
      double precision              :: t
      double precision              :: dt
      double precision              :: dx 
      double precision              :: I_0_gamma_x
      double precision              :: t_range
      integer, parameter            :: nt = 1000


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

      t_range = 20.0d0 
      dt = t_range/dble(nt)


!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,t,dx,I_0_gamma_x) COLLAPSE(1) SCHEDULE(static)

      do i = 0, nt-1

        t =  i * dt

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(X))) / ax**2


        I_0_gamma_x = gsl_sf_bessel_I0_scaled(dx)
        
        !I_0_gamma_x = I_0_gamma_x * exp(dx)

        !if (ieee_is_finite(I_0_gamma_x) .eqv. .false.) then
        !  exit
        !end if

        sum_f = sum_f + 1.d0/(p+t**2) * dexp(-2.d0*(t**2+p)/ax**2 + dx)  * I_0_gamma_x

        !print*, t , dexp(-2.d0*(t**2+p)/ax**2 + dx)  , I_0_gamma_x

      end do

!$OMP END PARALLEL DO

      ! Calculate the final result - multiply by volume element
      result =  dt *  sum_f
    
end subroutine grid_integrate_NA_ss_mod

subroutine integrate_NA_ss_mod(gamma_x,xAB,p, result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , XAB , p 
      
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

        double precision :: I_0_x ,dx 

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

        dx = 2.d0*dsqrt( gamma_x**2 + x**4 + 2.d0 * gamma_x * x**2 * dcos(ax*(XAB))) / ax**2
        I_0_x = gsl_sf_bessel_I0_scaled(dx)

        fx  = 1.d0/(p+x**2) * dexp(-2.d0*(x**2+p)/ax**2 + dx)  * I_0_x

      end function f_decay

end subroutine integrate_NA_ss_mod



subroutine grid_integrate_NA_sp(gamma_x, gamma_y, gamma_z, xp, x1, x2, xc, yc, zc, result)
    
      use omp_lib
      use torus_init
    
      implicit none
    
      ! Input parameters
      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
      double precision, intent(in) :: xp, x1, x2
      double precision, intent(in) :: xc, yc, zc 

      ! Grid parameters
      integer, parameter :: nx = 250      ! Number of grid points in x direction
      integer, parameter :: ny = 100      ! Number of grid points in y direction
      integer, parameter :: nz = 100      ! Number of grid points in z direction

      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables
      integer          :: i, j, k
      double precision :: sum_f
      double precision :: x, y, z
      double precision :: dx, dy, dz
      double precision :: x_range, y_range, z_range
      double precision :: factor_x, factor_y, factor_z
      double precision :: denominator

      ! Pre-compute constant factors
      factor_x = gamma_x * (2.0d0/ax**2)
      factor_y = gamma_y 
      factor_z = gamma_z

      ! Calculate grid spacing

      x_range = Lx    
      y_range = 5.0d0 / dsqrt(gamma_y) 
      z_range = 5.0d0 / dsqrt(gamma_z) 
      
      dx = Lx   / dble(nx)
      dy = (2.0d0 * y_range) / dble(ny)
      dz = (2.0d0 * z_range) / dble(nz)

      ! Initialize sum

      sum_f = 0.0d0


!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,j,k,x,y,z,denominator) COLLAPSE(3) SCHEDULE(static)
      do k = 0, nz-1
          do j = 0, ny-1
              do i = 0, nx-1

                  x =  i * dx
                  y =  - y_range + j * dy 
                  z =  - z_range + k * dz

                  ! Calculate denominator

                  denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x-xc))) + (y-yc)**2 + (z-zc)**2

                  if (denominator > 1.0d-10) then
                    sum_f = sum_f + sin(ax*(x-x2)) * dexp((factor_x)*cos(ax*(x-xp))) * dexp(-factor_y*y**2) * dexp(-factor_z*z**2) / sqrt(denominator)
                  end if

              end do
          end do
      end do
!$OMP END PARALLEL DO

      ! Calculate the final result - multiply by volume element
      result = dx * dy * dz * sum_f
    
end subroutine grid_integrate_NA_sp


subroutine integrate_NA_sp_mod(gamma_x,xP,xc,xb,p, result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , xP , p , xc , xb 
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(x) result(fx)

        double precision, intent(in) :: x
        double precision             :: fx

        double precision             :: I_1_x , xD , nu_x 

        INTERFACE
      
      FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
        USE iso_c_binding
        REAL(C_DOUBLE), VALUE :: x_val
        REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
      END FUNCTION gsl_sf_bessel_I1_scaled
      
        END INTERFACE

        xD   = datan((gamma_x*dsin(ax*xp)+x**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+x**2.d0*dcos(ax*xc)))/ax + 0.5*Lx * Heaviside(-gamma_x*cos(ax*xp)-x**2.d0*cos(ax*xc))

        nu_x = dsqrt(gamma_x**2 + x**4 + 2.d0 * gamma_x * x**2 * dcos(ax*(xp-xc)))

        I_1_x = gsl_sf_bessel_I1_scaled(2.d0*nu_x/ax**2)

        fx  = 1.d0/(p+x**2) * dexp(-2.d0*(x**2-nu_x)/ax**2) * dsin(ax*(xD-xB)) * I_1_x

      end function f_decay
    
end subroutine integrate_NA_sp_mod


!subroutine grid_integrate_NA_pp_px(gamma_x, gamma_y, gamma_z, xp, x1, x2, xc, yc, zc, result)
!    
!      use omp_lib
!      use torus_init
!    
!      implicit none
!    
!      ! Input parameters
!      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
!      double precision, intent(in) :: xp, x1, x2
!      double precision, intent(in) :: xc, yc, zc 
!
!      ! Grid parameters
!      integer, parameter :: nx = 200      ! Number of grid points in x direction
!      integer, parameter :: ny = 100      ! Number of grid points in y direction
!      integer, parameter :: nz = 100      ! Number of grid points in z direction
!
!      ! Output parameters
!      double precision, intent(out) :: result
!    
!      ! Local variables
!      integer          :: i, j, k
!      double precision :: sum_f
!      double precision :: x, y, z
!      double precision :: dx, dy, dz
!      double precision :: x_range, y_range, z_range
!      double precision :: factor_x, factor_y, factor_z
!      double precision :: denominator
!
!      ! Pre-compute constant factors
!      factor_x = gamma_x * (2.0d0/ax**2)
!      factor_y = gamma_y 
!      factor_z = gamma_z
!
!      ! Calculate grid spacing
!
!      x_range = Lx    
!      y_range = 5.0d0 / dsqrt(gamma_y) 
!      z_range = 5.0d0 / dsqrt(gamma_z) 
!      
!      dx = Lx   / dble(nx)
!      dy = (2.0d0 * y_range) / dble(ny)
!      dz = (2.0d0 * z_range) / dble(nz)
!
!      ! Initialize sum
!
!      sum_f = 0.0d0
!
!
!!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,j,k,x,y,z,denominator) COLLAPSE(3) SCHEDULE(static)
!      do k = 0, nz-1
!          do j = 0, ny-1
!              do i = 0, nx-1
!
!                  x =  i * dx
!                  y =  - y_range + j * dy 
!                  z =  - z_range + k * dz
!
!                  ! Calculate denominator
!
!                  denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x-xc))) + (y-yc)**2 + (z-zc)**2
!
!                  if (denominator > 1.0d-10) then
!                    sum_f = sum_f + sin(ax*(x-x2)) * sin(ax*(x-x1)) * dexp((factor_x)*cos(ax*(x-xp))) * dexp(-factor_y*y**2) * dexp(-factor_z*z**2) / sqrt(denominator)
!                  end if
!
!              end do
!          end do
!      end do
!!$OMP END PARALLEL DO
!
!      ! Calculate the final result - multiply by volume element
!      result = dx * dy * dz * sum_f
!    
!end subroutine grid_integrate_NA_pp_px


subroutine integrate_NA_pp_px_mod(gamma_x,xP,p,xc,xa,xb,result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , xP , p , xc , xb , xa  
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(x) result(fx)

        use gsl_bessel_mod
        double precision, intent(in) :: x
        double precision             :: fx

        double precision             :: xD    , nu_x 
        double precision             :: I_0_x , I_2_x


        xD   = datan((gamma_x*dsin(ax*xp)+x**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+x**2.d0*dcos(ax*xc)))/ax + 0.5*Lx * Heaviside(-gamma_x*cos(ax*xp)-x**2.d0*cos(ax*xc))

        nu_x = dsqrt(gamma_x**2 + x**4 + 2.d0 * gamma_x * x**2 * dcos(ax*(xp-xc)))
        
        I_0_x = bessi_scaled(0, 2.d0*nu_x/ax**2)
        I_2_x = bessi_scaled(2, 2.d0*nu_x/ax**2)

        fx  = 1.d0/(p+x**2) * dexp(-2.d0*(x**2-nu_x)/ax**2) * ( cos(ax*(xb-xa)) * I_0_x - cos(ax*(2.d0*xD-xa-xb)) *  I_2_x  )

      end function f_decay
    
end subroutine integrate_NA_pp_px_mod



subroutine integrate_NA_pp_py_mod(gamma_x,xpc,p, result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , XPC , p 
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables
      double precision,parameter    :: epsabs = 1.0e-8 , epsrel = 1.0e-6
      integer,parameter             :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer, parameter            :: limit = 100
      integer, parameter            :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(x) result(fx)
        double precision, intent(in) :: x
        double precision     :: fx

        double precision :: I_0_x ,dx 

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

        dx = 2.d0*dsqrt( gamma_x**2 + x**4 + 2.d0 * gamma_x * x**2 * dcos(ax*(XPC))) / ax**2
        I_0_x = gsl_sf_bessel_I0_scaled(dx)

        fx  = (1.d0/(p+x**2))**2  * dexp(-2.d0*(x**2+p)/ax**2 + dx)  * I_0_x

      end function f_decay
    
end subroutine integrate_NA_pp_py_mod






!subroutine grid_integrate_NA_pp_py(gamma_x, gamma_y, gamma_z, xp, x1, x2, xc, yc, zc, result)
!    
!      use omp_lib
!      use torus_init
!    
!      implicit none
!    
!      ! Input parameters
!      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
!      double precision, intent(in) :: xp, x1, x2
!      double precision, intent(in) :: xc, yc, zc 
!
!      ! Grid parameters
!      integer, parameter :: nx = 200      ! Number of grid points in x direction
!      integer, parameter :: ny = 100      ! Number of grid points in y direction
!      integer, parameter :: nz = 100      ! Number of grid points in z direction
!
!      ! Output parameters
!      double precision, intent(out) :: result
!    
!      ! Local variables
!      integer          :: i, j, k
!      double precision :: sum_f
!      double precision :: x, y, z
!      double precision :: dx, dy, dz
!      double precision :: x_range, y_range, z_range
!      double precision :: factor_x, factor_y, factor_z
!      double precision :: denominator
!
!      ! Pre-compute constant factors
!      factor_x = gamma_x * (2.0d0/ax**2)
!      factor_y = gamma_y 
!      factor_z = gamma_z
!
!      ! Calculate grid spacing
!
!      x_range = Lx    
!      y_range = 5.0d0 / dsqrt(gamma_y) 
!      z_range = 5.0d0 / dsqrt(gamma_z) 
!      
!      dx = Lx   / dble(nx)
!      dy = (2.0d0 * y_range) / dble(ny)
!      dz = (2.0d0 * z_range) / dble(nz)
!
!      ! Initialize sum
!
!      sum_f = 0.0d0
!
!
!!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,j,k,x,y,z,denominator) COLLAPSE(3) SCHEDULE(static)
!      do k = 0, nz-1
!          do j = 0, ny-1
!              do i = 0, nx-1
!
!                  x =  i * dx
!                  y =  - y_range + j * dy 
!                  z =  - z_range + k * dz
!
!                  ! Calculate denominator
!
!                  denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x-xc))) + (y-yc)**2 + (z-zc)**2
!
!                  if (denominator > 1.0d-10) then
!                    sum_f = sum_f +  y**2 * dexp((factor_x)*cos(ax*(x-xp))) * dexp(-factor_y*y**2) * dexp(-factor_z*z**2) / sqrt(denominator)
!                  end if
!
!              end do
!          end do
!      end do
!!$OMP END PARALLEL DO
!
!      ! Calculate the final result - multiply by volume element
!      result = dx * dy * dz * sum_f
!    
!end subroutine grid_integrate_NA_pp_py
!
!
!subroutine grid_integrate_NA_pp_pz(gamma_x, gamma_y, gamma_z, xp, x1, x2, xc, yc, zc, result)
!    
!      use omp_lib
!      use torus_init
!    
!      implicit none
!    
!      ! Input parameters
!      double precision, intent(in) :: gamma_x, gamma_y, gamma_z
!      double precision, intent(in) :: xp, x1, x2
!      double precision, intent(in) :: xc, yc, zc 
!
!      ! Grid parameters
!      integer, parameter :: nx = 200      ! Number of grid points in x direction
!      integer, parameter :: ny = 100      ! Number of grid points in y direction
!      integer, parameter :: nz = 100      ! Number of grid points in z direction
!
!      ! Output parameters
!      double precision, intent(out) :: result
!    
!      ! Local variables
!      integer          :: i, j, k
!      double precision :: sum_f
!      double precision :: x, y, z
!      double precision :: dx, dy, dz
!      double precision :: x_range, y_range, z_range
!      double precision :: factor_x, factor_y, factor_z
!      double precision :: denominator
!
!      ! Pre-compute constant factors
!      factor_x = gamma_x * (2.0d0/ax**2)
!      factor_y = gamma_y 
!      factor_z = gamma_z
!
!      ! Calculate grid spacing
!
!      x_range = Lx    
!      y_range = 5.0d0 / dsqrt(gamma_y) 
!      z_range = 5.0d0 / dsqrt(gamma_z) 
!      
!      dx = Lx   / dble(nx)
!      dy = (2.0d0 * y_range) / dble(ny)
!      dz = (2.0d0 * z_range) / dble(nz)
!
!      ! Initialize sum
!
!      sum_f = 0.0d0
!
!
!!$OMP PARALLEL DO REDUCTION(+:sum_f) PRIVATE(i,j,k,x,y,z,denominator) COLLAPSE(3) SCHEDULE(static)
!      do k = 0, nz-1
!          do j = 0, ny-1
!              do i = 0, nx-1
!
!                  x =  i * dx
!                  y =  - y_range + j * dy 
!                  z =  - z_range + k * dz
!
!                  ! Calculate denominator
!
!                  denominator = (1.0d0/ax**2)*(2.0d0 - 2.0d0*cos(ax*(x-xc))) + (y-yc)**2 + (z-zc)**2
!
!                  if (denominator > 1.0d-10) then
!                    sum_f = sum_f +  z**2 * dexp((factor_x)*cos(ax*(x-xp))) * dexp(-factor_y*y**2) * dexp(-factor_z*z**2) / sqrt(denominator)
!                  end if
!
!              end do
!          end do
!      end do
!!$OMP END PARALLEL DO
!
!      ! Calculate the final result - multiply by volume element
!      result = dx * dy * dz * sum_f
!    
!end subroutine grid_integrate_NA_pp_pz
