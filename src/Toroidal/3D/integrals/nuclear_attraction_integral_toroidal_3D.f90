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
            
            !gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
            !gamma_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
            !gamma_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

            call bary_exponent(alpha,beta,X,gamma_x)
            call bary_exponent(alpha,beta,Y,gamma_y)
            call bary_exponent(alpha,beta,Z,gamma_z)

            !xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            !yp       = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5*Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
            !zp       = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5*Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

            call bary_center_toroidal(alpha,beta,x1,x2,xp)
            call bary_center_toroidal(alpha,beta,y1,y2,yp)
            call bary_center_toroidal(alpha,beta,z1,z2,zp)

            kc       = c1 * c2 * 2.d0/dsqrt(pi) * Lx * Ly * Lz

            do k = 1 , number_of_atoms

              xc = geometry(k,1) 
              yc = geometry(k,2)
              zc = geometry(k,3)
              
              charge_atom = (-1)*atoms(k)%charge
              
               call integrate_NA_ss_Toroidal_3D(xp-xc,yp-yc,zp-zc,&
                                                gamma_x, gamma_y , gamma_z ,&
                                                alpha+beta,NA)

              S_ss_normal =  S_ss_normal +  charge_atom * kc * NA
            
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
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters

      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe
      
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables
      
      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer,parameter             :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer, parameter            :: limit = 50
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
        double precision             :: I_0_x , I_0_y , I_0_z  
        double precision             :: dx , dy , dz  
        double precision             :: Nax , Nay , Naz  
        
        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)

        I_0_x = iv_scaled(0, dx)
        I_0_y = iv_scaled(0, dy)
        I_0_z = iv_scaled(0, dz)


        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x 
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z
        
        fx    =  Nax * Nay * Naz  

      end function f_decay

end subroutine integrate_NA_ss_Toroidal_3D



subroutine nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_sp_normal)

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
            
            !gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ax*(X)))
            !gamma_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(ay*(Y)))
            !gamma_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*dcos(az*(Z)))

            call bary_exponent(alpha,beta,X,gamma_x)
            call bary_exponent(alpha,beta,Y,gamma_y)
            call bary_exponent(alpha,beta,Z,gamma_z)

            !xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5d0 * Lx * Heaviside(-alpha*dcos(ax*x1)-beta*dcos(ax*x2))
            !yp       = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5d0 * Ly * Heaviside(-alpha*dcos(ay*y1)-beta*dcos(ay*y2))
            !zp       = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5d0 * Lz * Heaviside(-alpha*dcos(az*z1)-beta*dcos(az*z2))

            call bary_center_toroidal(alpha,beta,x1,x2,xp)
            call bary_center_toroidal(alpha,beta,y1,y2,yp)
            call bary_center_toroidal(alpha,beta,z1,z2,zp)

            kc       = c1 * c2  * 2.d0/dsqrt(pi) * Lx * Ly * Lz

            if (AO2%orbital=="px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
                  
                call integrate_NA_spx_Toroidal_3D(xp-xc,yp-yc,zp-zc,&
                                                  gamma_x, gamma_y , gamma_z ,&
                                                  xp , xc , x2 ,&
                                                  alpha+beta,NA)
              
                S_sp_normal =  S_sp_normal + charge_atom * kc * NA / ax

              end do

            end if

            if (AO2%orbital=="py") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
                  
                call integrate_NA_spy_Toroidal_3D(xp-xc,yp-yc,zp-zc,&
                                                  gamma_x, gamma_y , gamma_z ,&
                                                  yp , yc , y2 ,&
                                                  alpha+beta,NA)
              
                S_sp_normal =  S_sp_normal + charge_atom * kc * NA / ay

              end do

            end if

            if (AO2%orbital=="pz") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
                  
                call integrate_NA_spz_Toroidal_3D(xp-xc,yp-yc,zp-zc,&
                                                  gamma_x, gamma_y , gamma_z ,&
                                                  zp , zc , z2 ,&
                                                  alpha+beta,NA)
              
                S_sp_normal =  S_sp_normal + charge_atom * kc * NA / az

              end do

            end if


        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_sp_toroidal_3D


subroutine integrate_NA_spx_Toroidal_3D(xpC,ypC,zpC, &
                                       gamma_x, gamma_y, gamma_z, &
                                       xp , xc , xB ,&
                                       albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: XpC , YpC , ZpC
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z
      double precision, intent(in)  :: albe
      double precision, intent(in)  :: xp , xc , xB 
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: dx    , dy    , dz
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz

        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)

        I_0_x = iv_scaled(0, dx)
        I_0_y = iv_scaled(0, dy)
        I_0_z = iv_scaled(0, dz)

        I_1_x = iv_scaled(1, dx)
        I_1_y = iv_scaled(1, dy)
        I_1_z = iv_scaled(1, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xB))
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx    =  Nax * Nay * Naz 

      end function f_decay
    
end subroutine integrate_NA_spx_Toroidal_3D

subroutine integrate_NA_spy_Toroidal_3D(xpC,ypC,zpC, &
                                       gamma_x, gamma_y, gamma_z, &
                                       yp , yc , yB ,&
                                       albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: XpC , YpC , ZpC
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z
      double precision, intent(in)  :: albe
      double precision, intent(in)  :: yp , yc , yB 
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: dx    , dy    , dz
        double precision             :: yD
        double precision             :: Nax   , Nay   , NAz

        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))

        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)

        I_0_x = iv_scaled(0, dx)
        I_0_y = iv_scaled(0, dy)
        I_0_z = iv_scaled(0, dz)

        I_1_x = iv_scaled(1, dx)
        I_1_y = iv_scaled(1, dy)
        I_1_z = iv_scaled(1, dz)

        Nax = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_1_y * dsin(ay*(yD-yB))
        Naz = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx  =  Nax * Nay * Naz 

      end function f_decay
    
end subroutine integrate_NA_spy_Toroidal_3D


subroutine integrate_NA_spz_Toroidal_3D(xpC,ypC,zpC, &
                                       gamma_x, gamma_y, gamma_z, &
                                       zp , zc , zB ,&
                                       albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: XpC , YpC , ZpC
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z
      double precision, intent(in)  :: albe
      double precision, intent(in)  :: zp , zc , zB 
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: dx    , dy    , dz
        double precision             :: zD 
        double precision             :: Nax   , Nay   , NAz

        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))

        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)

        I_0_x = iv_scaled(0, dx)
        I_0_y = iv_scaled(0, dy)
        I_0_z = iv_scaled(0, dz)

        I_1_x = iv_scaled(1, dx)
        I_1_y = iv_scaled(1, dy)
        I_1_z = iv_scaled(1, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_1_z * dsin(az*(zD-zB))

        fx    =  Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_spz_Toroidal_3D

subroutine nuclear_attraction_integral_pp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use atom_basis
      use classification_ERI
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
            
            !gamma_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X)))
            !gamma_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(Y)))
            !gamma_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(Z)))

            call bary_exponent(alpha,beta,X,gamma_x)
            call bary_exponent(alpha,beta,Y,gamma_y)
            call bary_exponent(alpha,beta,Z,gamma_z)

            !xp       = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5d0 * Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))
            !yp       = datan((alpha*dsin(ay*y1)+beta*dsin(ay*y2))/(alpha*dcos(ay*y1)+beta*dcos(ay*y2)))/ay + 0.5d0 * Ly * Heaviside(-alpha*cos(ay*y1)-beta*cos(ay*y2))
            !zp       = datan((alpha*dsin(az*z1)+beta*dsin(az*z2))/(alpha*dcos(az*z1)+beta*dcos(az*z2)))/az + 0.5d0 * Lz * Heaviside(-alpha*cos(az*z1)-beta*cos(az*z2))

            call bary_center_toroidal(alpha,beta,x1,x2,xp)
            call bary_center_toroidal(alpha,beta,y1,y2,yp)
            call bary_center_toroidal(alpha,beta,z1,z2,zp)

            kc       = c1 * c2 * 2.d0/dsqrt(pi) * Lx * Ly * Lz

            if (AO1%orbital == "px" .and. AO2%orbital == "px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_px_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax / ax 
              
              end do

            end if 

            if (AO1%orbital == "px" .and. AO2%orbital == "py") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_py_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax / ay 
              
              end do

            end if

            if (AO1%orbital == "px" .and. AO2%orbital == "pz") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_pz_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax / az
              
              end do

            end if


            if (AO1%orbital == "py" .and. AO2%orbital == "px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_px_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / ay / ax 
              
              end do

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "py") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_py_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / ay / ay 
              
              end do

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "pz") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_pz_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / ay / az
              
              end do

            end if


            if (AO1%orbital == "pz" .and. AO2%orbital == "px") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                 
                charge_atom = (-1)*atoms(k)%charge
   
                call integrate_NA_pz_px_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / az / ax 
              
              end do

            end if


            if (AO1%orbital == "pz" .and. AO2%orbital == "py") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                 
                charge_atom = (-1)*atoms(k)%charge
   
                call integrate_NA_pz_py_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / az / ay 
              
              end do

            end if
            
            if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 
            
              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
                 
                charge_atom = (-1)*atoms(k)%charge
   
                call integrate_NA_pz_pz_Toroidal_3D(xp-xc,yp-yc,zp-zc, &
                                                    gamma_x, gamma_y, gamma_z, &
                                                    xp   , xc  , x1  , x2 , &
                                                    yp   , yc  , y1  , y2 , &
                                                    zp   , zc  , z1  , z2 , &
                                                    alpha+beta , NA)
              
                S_pp_normal =  S_pp_normal +  charge_atom * kc * NA / az / az 
              
              end do

            end if

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_pp_toroidal_3D

subroutine integrate_NA_px_px_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions 
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   = 0.5d0 * dexp(-2.d0*(t**2+albe)/ax**2 + dx) * ( dcos(ax*(xb-xa)) * I_0_x - dcos(ax*(2.d0*xD-xa-xb)) *  I_2_x  )
        Nay   =         dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   =         dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_px_Toroidal_3D

subroutine integrate_NA_px_py_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))
        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))

        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)
        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)


        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xA))
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_1_y * dsin(ay*(yD-yB))
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_py_Toroidal_3D

subroutine integrate_NA_px_pz_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))
        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))
        
        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)
        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xA))
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_1_z * dsin(az*(zD-zB))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_pz_Toroidal_3D

subroutine integrate_NA_py_px_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))
        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5*Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)
        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)
        
        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        
        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xB))
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_1_y * dsin(ay*(yD-yA))
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_px_Toroidal_3D

subroutine integrate_NA_py_py_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))

        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)


        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)
        
        Nax   =         dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = 0.5d0 * dexp(-2.d0*(t**2+albe)/ay**2 + dy) * ( dcos(ay*(yb-ya)) * I_0_y - dcos(ay*(2.d0*yD-ya-yb)) *  I_2_y  )
        Naz   =         dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_0_z

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_py_Toroidal_3D

subroutine integrate_NA_py_pz_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))
        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))

        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)
        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        
        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_1_y * dsin(ay*(yD-yA))
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_1_z * dsin(az*(zD-zB))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_pz_Toroidal_3D

subroutine integrate_NA_pz_px_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))
        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))

        call bary_center_toroidal(gamma_x,t*t,xp,xc,xD)
        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)


        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xB))
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_1_z * dsin(az*(zD-zA))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_px_Toroidal_3D


subroutine integrate_NA_pz_py_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 


        !yD   = datan((gamma_y*dsin(ay*yp)+t**2.d0*dsin(ay*yc))/(gamma_y*dcos(ay*yp)+t**2.d0*dcos(ay*yc)))/ay + 0.5d0 * Ly * Heaviside(-gamma_y*dcos(ay*yp)-t**2.d0*dcos(ay*yc))
        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))

        call bary_center_toroidal(gamma_y,t*t,yp,yc,yD)
        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)
        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_1_y * dsin(ay*(yD-yB))
        Naz   = dexp(-2.d0*(t**2+albe)/az**2 + dz) * I_1_z * dsin(az*(zD-zA))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_py_Toroidal_3D


subroutine integrate_NA_pz_pz_Toroidal_3D(xpC,ypC,zpC, &
                                          gamma_x, gamma_y, gamma_z, &
                                          xp   , xc  , xa  , xb , &
                                          yp   , yc  , ya  , yb , &
                                          zp   , zc  , za  , zb , &
                                          albe , result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use HeavisideModule
      use, intrinsic :: ieee_arithmetic
      
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  ::     XpC ,     YpC ,     ZpC 
      double precision, intent(in)  :: gamma_x , gamma_y , gamma_z 
      double precision, intent(in)  :: albe 
      double precision, intent(in)  :: xp , xc , xa , xb
      double precision, intent(in)  :: yp , yc , ya , yb
      double precision, intent(in)  :: zp , zc , za , zb
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-10 , epsrel = 1.0e-8
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 50
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(fx)

        use bessel_functions

        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x , I_0_y , I_0_z
        double precision             :: I_1_x , I_1_y , I_1_z
        double precision             :: I_2_x , I_2_y , I_2_z

        double precision             :: dx    , dy    , dz
        double precision             :: xD    , yD    , zD
        double precision             :: Nax   , Nay   , NAz 

        !zD   = datan((gamma_z*dsin(az*zp)+t**2.d0*dsin(az*zc))/(gamma_z*dcos(az*zp)+t**2.d0*dcos(az*zc)))/az + 0.5d0 * Lz * Heaviside(-gamma_z*dcos(az*zp)-t**2.d0*dcos(az*zc))

        call bary_center_toroidal(gamma_z,t*t,zp,zc,zD)

        !dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        !dy = 2.d0*dsqrt( gamma_y**2 + t**4 + 2.d0 * gamma_y * t**2 * dcos(ay*(ypC))) / ay**2
        !dz = 2.d0*dsqrt( gamma_z**2 + t**4 + 2.d0 * gamma_z * t**2 * dcos(az*(zpC))) / az**2

        call bary_exponent(gamma_x,t*t,XpC,dx)
        call bary_exponent(gamma_y,t*t,YpC,dy)
        call bary_exponent(gamma_z,t*t,ZpC,dz)

        dx = 2.d0 * dx / (ax*ax)
        dy = 2.d0 * dy / (ay*ay)
        dz = 2.d0 * dz / (az*az)


        
        I_0_x = iv_scaled(0, dx)
        I_1_x = iv_scaled(1, dx)
        I_2_x = iv_scaled(2, dx)

        I_0_y = iv_scaled(0, dy)
        I_1_y = iv_scaled(1, dy)
        I_2_y = iv_scaled(2, dy)

        I_0_z = iv_scaled(0, dz)
        I_1_z = iv_scaled(1, dz)
        I_2_z = iv_scaled(2, dz)

        Nax   =         dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   =         dexp(-2.d0*(t**2+albe)/ay**2 + dy) * I_0_y
        Naz   = 0.5d0 * dexp(-2.d0*(t**2+albe)/az**2 + dz) * ( dcos(az*(zb-za)) * I_0_z - dcos(az*(2.d0*zD-za-zb)) *  I_2_z  )

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_pz_Toroidal_3D