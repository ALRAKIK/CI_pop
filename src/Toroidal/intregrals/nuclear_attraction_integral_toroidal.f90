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
            
            kc       = dexp(-(alpha+beta)*Lx**2/(2.d0*pi**2))

            gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))+eta
            gamma_y  = alpha+beta
            gamma_z  = alpha+beta

            xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            yp       = 0.d0
            zp       = 0.d0


            do k = 1 , number_of_atoms

              xc = geometry(k,1) 
              yc = geometry(k,2)
              zc = geometry(k,3)
              
              charge_atom = (-1)*atoms(k)%charge

              call grid_integrate_NA(gamma_x,gamma_y,gamma_z,xp,yp,zp,xc,yc,zc,NA)
            
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
      double precision              :: p,mu
      double precision              :: Two_PIP , P_R
      double precision              :: x1 , x2 , X
      double precision              :: y1 , y2 , Y 
      double precision              :: z1 , z2 , Z 
      double precision              :: xp , yp , zp 
      double precision              :: x_PB_normal , y_PB_normal , z_PB_normal
      double precision              :: D_normal 
      double precision              :: x_PC , y_PC , z_PC 
      double precision              :: Boys_func , R2PC
      

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      call PBC(x1,x2,X)
!      call EUC(x1,x2,X)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!
 
      S_sp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)
            
            p       = alpha + beta 
            P_R     = 1.d0 / p
            mu      = alpha*beta*p_R
            Two_PIP = 2.0d0*pi*p_R

            xp = (alpha*x1+beta*x2)*p_R
            
            call bary_center(alpha,x1,beta,x2,p_R,xp)

            yp = (alpha*y1+beta*y2)*p_R
            zp = (alpha*z1+beta*z2)*p_R

            X = x1 - x2

            call SSD(x1,x2,X)
!            call EUC(x1,x2,X)

            X_PB_normal  =  (alpha/p)*(X)
            Y_PB_normal  =  (alpha/p)*(Y)
            Z_PB_normal  =  (alpha/p)*(Z)

            X = x1 - x2 

!            if (abs(X)  < 1e-10) then 

!            else 

              do k = 1 , number_of_atoms
              
                x_PC = xp - geometry(k,1) 
                y_PC = yp - geometry(k,2)
                z_PC = zp - geometry(k,3)
  
                if (torus) call euc(xp,geometry(k,1),x_PC)


                R2PC = x_PC*x_PC + y_PC*y_PC + z_PC*z_PC
              
                charge_atom = (-1)*atoms(k)%charge
              
                if (AO2%orbital=="px") S_sp_normal =  S_sp_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (X_PB_normal * Boys_func(0,p*R2PC) - X_PC * Boys_func(1,p*R2PC) ) 
                if (AO2%orbital=="py") S_sp_normal =  S_sp_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Y_PB_normal * Boys_func(0,p*R2PC) - Y_PC * Boys_func(1,p*R2PC) ) 
                if (AO2%orbital=="pz") S_sp_normal =  S_sp_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Z_PB_normal * Boys_func(0,p*R2PC) - Z_PC * Boys_func(1,p*R2PC) ) 
            
              end do  
!            end if 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_sp_toroidal


subroutine grid_integrate_NA(gamma_x, gamma_y, gamma_z, xp, yp, zp, xc, yc, zc, result)
    
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
    
end subroutine grid_integrate_NA
