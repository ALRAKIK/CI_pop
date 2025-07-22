subroutine nuclear_attraction_integral_ss_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_ss_normal)

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
            
            gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
            gamma_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
            gamma_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

            xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            yp       = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
            zp       = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

            kc       = 2.d0/dsqrt(pi) * Lx * Ly * Lz 

            do k = 1 , number_of_atoms

              xc = geometry(k,1) 
              yc = geometry(k,2)
              zc = geometry(k,3)
              
              charge_atom = (-1)*atoms(k)%charge
              
               call integrate_NA_ss_Toroidal_3D(xp-xc,yp-yc,zp-zc,&
                                                gamma_x, gamma_y , gamma_z ,&
                                                alpha+beta,NA)

              S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * kc * NA
            
            end do 
            
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_ss_toroidal_3D

subroutine integrate_NA_ss_Toroidal_3D(xpC,ypC,zpC, &
                                       gamma_x, gamma_y, gamma_z, &
                                       albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use gsl_bessel_mod
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe
      
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
        write(*,'(A,//)')   'Error code from ss NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        double precision, intent(in) :: t
        double precision             :: fx
        double precision             :: I_0_x ,dx 
        double precision             :: I_0_y ,dy
        double precision             :: I_0_z ,dz

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        I_0_x = bessi_scaled(0, dx)
        I_0_y = bessi_scaled(0, dy)
        I_0_z = bessi_scaled(0, dz)
        
        fx  =  dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x * dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y * dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z 

      end function f_decay

end subroutine integrate_NA_ss_Toroidal_3D