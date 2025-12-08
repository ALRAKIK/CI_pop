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
      double precision               :: x1 , x2 , X
      double precision               :: y1 , y2 , Y 
      double precision               :: z1 , z2 , Z
      double precision               :: xc , yc , zc 
      double precision               :: xp_C 
      double precision               :: yp_R , zp_R
      double precision               :: gamma_x
      double precision               :: NA

      double precision               :: albe 
      double precision               :: inv_albe
      double precision               :: mu

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X  = (x1 - x2)
      Y  = (y1 - y2)
      Z  = (z1 - z2)

      !-----------------------------------------------------------------!
      
      S_ss_normal = 0.d0
      NA          = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

          ! Clifford gaussian ! 

          gamma_x = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          xp_C    = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))

          !   Real Gaussian   !
            
          albe     = alpha + beta 
          inv_albe = 1.d0 / albe 
          mu       = alpha * beta * inv_albe
            
          yp_R       = (alpha * y1 + beta * y2) * inv_albe
          zp_R       = (alpha * z1 + beta * z2) * inv_albe

          kc       = 2.d0*sqrt(pi) * Lx * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))

          do k = 1 , number_of_atoms

            xc = geometry(k,1) 
            yc = geometry(k,2)
            zc = geometry(k,3)
              
            charge_atom = (-1)*atoms(k)%charge
              
            call integrate_NA_ss_Toroidal(gamma_x,xp_C-xc,yp_R-yc,zp_R-zc,albe, NA)

            S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * kc * NA 
            
          end do 
            
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_ss_toroidal

subroutine integrate_NA_ss_Toroidal(gamma_x,xPC,yPC,zPC,albe, result)
    
      use quadpack, only : dqagi
      use omp_lib
      use torus_init
      use gsl_bessel_mod
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters

      double precision, intent(in)  :: gamma_x , albe
      double precision, intent(in)  :: xPC , yPC , zPC
      
      ! Output parameters

      double precision, intent(out) :: result
    
      ! Local variables
      
      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
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

      function f_decay(t) result(ft)

        double precision, intent(in) :: t
        double precision             :: ft
        double precision             :: I_0_x ,dx 

        double precision             :: t2 ,t4 

        t2 = t  * t 
        t4 = t2 * t2 

        dx = 2.d0*dsqrt(dabs(gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XPC)))) / ax**2

        I_0_x = bessi_scaled(0, dx)
        
        ft  = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe)/ax**2 + dx)  * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) + (zPC*zPC) )  *  I_0_x 

      end function f_decay

end subroutine integrate_NA_ss_Toroidal


!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !


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
      double precision              :: xc , yc , zc
      double precision              :: x1 , x2 , X
      double precision              :: y1 , y2 , Y 
      double precision              :: z1 , z2 , Z 
      double precision              :: xp_C 
      double precision              :: yp_R , zp_R
      double precision              :: gamma_x , gamma_y , gamma_z
      double precision              :: NA
      double precision              :: albe 
      double precision              :: inv_albe
      double precision              :: mu
      
      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      !-----------------------------------------------------------------!
      
      S_sp_normal = 0.d0
      NA          = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)


          ! Clifford Gaussian ! 

          gamma_x = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          xp_C    = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))


          !  Real Gaussian    !  
          
          albe = alpha + beta 
          inv_albe = 1.d0/albe 
          mu = alpha * beta * inv_albe

          yp_R = (alpha * y1 + beta * y2) * inv_albe
          zp_R = (alpha * z1 + beta * z2) * inv_albe

          gamma_y  = alpha+beta
          gamma_z  = alpha+beta

          kc       = c1 * c2  * 2.d0 * dsqrt(pi) / ax * Lx * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))
  
          do k = 1 , number_of_atoms

            xc = geometry(k,1) 
            yc = geometry(k,2)
            zc = geometry(k,3)
                
            charge_atom = (-1)*atoms(k)%charge

            if (AO2%orbital=="px") then 

              kc       = c1 * c2  * 2.d0 * dsqrt(pi) / ax * Lx * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))
                  
              call integrate_NA_spx_Toroidal(gamma_x,xp_C,xc,x2,yp_R-yc,zp_R-zc,alpha+beta,NA)
              
            end if 

            if (AO2%orbital=="py") then 
                  
              kc       = c1 * c2  * 2.d0 * dsqrt(pi) * Lx * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))

              call integrate_NA_spy_Toroidal(gamma_x,xp_C-xc,yp_R-yc,yp_R-y2,zp_R-zc,albe, NA)
              
            end if

            if (AO2%orbital=="pz") then 
                  
              kc       = c1 * c2  * 2.d0 * dsqrt(pi) * Lx * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))

              call integrate_NA_spz_Toroidal(gamma_x,xp_C-xc,yp_R-yc,zp_R-zc,zp_R-z2,albe, NA)
              
            end if

            S_sp_normal =  S_sp_normal + charge_atom * kc * NA

          end do

            

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_sp_toroidal

subroutine integrate_NA_spx_Toroidal(gamma_x,xP,xc,xb,yPC,zPC,albe, result)
    
      use quadpack, only : dqagi
      use torus_init
      use gsl_bessel_mod
      use HeavisideModule
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , xP , albe , xc , xb 
      double precision, intent(in)  :: yPC , zPC
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)

        double precision, intent(in) :: t
        double precision             :: ft

        double precision             :: I_1_x , xD , nu_x , t2 

        t2   = t * t 

        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5*Lx * Heaviside(-gamma_x*cos(ax*xp)-t2*cos(ax*xc))

        nu_x = dsqrt(gamma_x**2 + t2 * t2  + 2.d0 * gamma_x * t2 * dcos(ax*(xp-xc)))

        I_1_x = bessi_scaled(1, 2.d0*nu_x/ax**2)

        ft  = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe-nu_x)/ax**2) * dsin(ax*(xD-xB)) * I_1_x * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) + (zPC*zPC) )

      end function f_decay
    
end subroutine integrate_NA_spx_Toroidal

subroutine integrate_NA_spy_Toroidal(gamma_x,xPc,yPC,yPB,zPC,albe, result)
    
      use quadpack, only : dqagi
      use torus_init
      use gsl_bessel_mod
      use HeavisideModule
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , xPc , albe
      double precision, intent(in)  :: yPC , zPC 
      double precision, intent(in)  :: yPB
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)

        double precision, intent(in) :: t
        double precision             :: ft
        double precision             :: I_0_x ,dx 

        double precision             :: t2 ,t4 

        t2 = t  * t 
        t4 = t2 * t2 

        dx = 2.d0*dsqrt(dabs(gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XPC)))) / ax**2

        I_0_x = bessi_scaled(0, dx)
        
        ft  = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe)/ax**2 + dx)  * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) + (zPC*zPC) )  *  I_0_x * (yPB - t2/(albe+t2) * (yPC))

      end function f_decay
    
end subroutine integrate_NA_spy_Toroidal


subroutine integrate_NA_spz_Toroidal(gamma_x,xPc,yPC,zPC,zPB,albe, result)
    
      use quadpack, only : dqagi
      use torus_init
      use gsl_bessel_mod
      use HeavisideModule
      use bessel_functions
      use, intrinsic :: ieee_arithmetic
    
      implicit none
    
      ! Input parameters
      double precision, intent(in)  :: gamma_x , xPc , albe
      double precision, intent(in)  :: yPC , zPC 
      double precision, intent(in)  :: zPB
      
      ! Output parameters
      double precision, intent(out) :: result
    
      ! Local variables

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from sp NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)

        double precision, intent(in) :: t
        double precision             :: ft
        double precision             :: I_0_x ,dx 

        double precision             :: t2 ,t4 

        t2 = t  * t 
        t4 = t2 * t2 

        dx = 2.d0*dsqrt(dabs(gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XPC)))) / ax**2

        I_0_x = bessi_scaled(0, dx)
        
        ft  = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe)/ax**2 + dx)  * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) + (zPC*zPC) )  *  I_0_x * (zPB - t2/(albe+t2) * (zPC)) 

      end function f_decay
    
end subroutine integrate_NA_spz_Toroidal

!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !
!////////////////////////////////////////////////////////////////////// !

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
      double precision              :: xc , yc , zc
      double precision              :: x1 , x2 , X
      double precision              :: y1 , y2 , Y 
      double precision              :: z1 , z2 , Z 
      double precision              :: xp_C , yp_R , zp_R 
      double precision              :: gamma_x , gamma_y , gamma_z
      double precision              :: NA

      double precision              :: albe 
      double precision              :: inv_albe
      double precision              :: mu
      
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

          !  Clifford Gaussian ! 
            
          gamma_x  = dsqrt(dabs(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(X))))
          xp_C     = datan((alpha*dsin(ax*x1)+beta*dsin(ax*x2))/(alpha*dcos(ax*x1)+beta*dcos(ax*x2)))/ax + 0.5*Lx * Heaviside(-alpha*cos(ax*x1)-beta*cos(ax*x2))

          !   Real Gaussian    !

          albe      = alpha + beta 
          inv_albe  = 1.d0/albe 
          mu        = alpha * beta * inv_albe
          yp_R      = (alpha * y1 + beta * y2) * inv_albe
          zp_R      = (alpha * z1 + beta * z2) * inv_albe

          ! ----------------------------------------------------------- !


          kc       = c1 * c2 * 2.d0 / dsqrt(pi) * Lx  * dexp(-mu * (Y*Y)) * dexp(-mu * (Z*Z))


            if (AO1%orbital == "px" .and. AO2%orbital == "px") then 

              do k = 1 , number_of_atoms

                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_px_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
                                                 
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax / ax 
              
              end do

            end if

            if (AO1%orbital == "px" .and. AO2%orbital == "py") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_py_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax
              
              end do

            end if

            if (AO1%orbital == "px" .and. AO2%orbital == "pz") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_px_pz_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax
              
              end do

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "px") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_px_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA / ax
              
              end do

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "py") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_py_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA
              
              end do

            end if

            if (AO1%orbital == "py" .and. AO2%orbital == "pz") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_py_pz_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA
              
              end do

            end if

            if (AO1%orbital == "pz" .and. AO2%orbital == "px") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_pz_px_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA /ax 
              
              end do

            end if

          
            if (AO1%orbital == "pz" .and. AO2%orbital == "py") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_pz_py_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA
              
              end do

            end if

            if (AO1%orbital == "pz" .and. AO2%orbital == "pz") then 

              do k = 1 , number_of_atoms
                
                xc = geometry(k,1) 
                yc = geometry(k,2)
                zc = geometry(k,3)
    
                charge_atom = (-1)*atoms(k)%charge
  
                call integrate_NA_pz_pz_Toroidal(xp_C-xc,yp_R-yc,zp_R-zc, &
                                                 gamma_x, gamma_y, gamma_z, &
                                                 xP_C , xc  , x1  , x2 , &
                                                 yp_R , yc  , y1  , y2 , &
                                                 zp_R , zc  , z1  , z2 , &
                                                 albe , NA)
              
                S_pp_normal =  S_pp_normal +   charge_atom * kc * NA
              
              end do

            end if

            

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_pp_toroidal


subroutine integrate_NA_px_px_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
      integer         ,parameter    :: lenw  = limit*4
      integer                       :: ier, iwork(limit), last, neval
      double precision              :: abserr, work(lenw)

      call dqagi(f_decay, bound, inf, epsabs, epsrel, result,abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(/,A,/)')  'Error code from ppx NA :'
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(t) result(ft)

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: ft

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 
        double precision             :: t2    , t4
        double precision,parameter    :: pi = 3.14159265358979323846D00

        t2 = t  * t 
        t4 = t2 * t2


        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t2*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = 0.5d0 * dexp(-2.d0*(t2+albe)/ax**2 + dx) * ( dcos(ax*(xb-xa)) * I_0_x - dcos(ax*(2.d0*xD-xa-xb)) *  I_2_x  )
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) )
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC) )

        ft    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_px_Toroidal

subroutine integrate_NA_px_py_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 

        double precision             :: t2   , t4
        double precision             :: yPB
        double precision,parameter   :: pi = 3.14159265358979323846D00

        t2 = t  * t
        t4 = t2 * t2

        yPB = yp-yb 

        xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xA))
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC)) * (yPB - t2/(albe+t2) * (yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_py_Toroidal


subroutine integrate_NA_px_pz_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 

        double precision             :: t2   , t4
        double precision             :: zPB
        double precision,parameter   :: pi = 3.14159265358979323846D00

        t2 = t  * t
        t4 = t2 * t2

        zPB = zp-zb 

        xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xA))
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC)) * (zPB - t2/(albe+t2) * (zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_px_pz_Toroidal


subroutine integrate_NA_py_px_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 

        double precision             :: t2   , t4
        double precision             :: yPA
        double precision,parameter   :: pi = 3.14159265358979323846D00

        t2 = t  * t
        t4 = t2 * t2

        yPA = yp-ya

        xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xB))
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC)) * (yPA - t2/(albe+t2) * (yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_px_Toroidal


subroutine integrate_NA_py_py_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 
        double precision             :: t2    , t4
        double precision             :: yPA , yPB
        double precision,parameter   :: pi = 3.14159265358979323846D00

        yPA = yp-ya
        yPB = yp-yb

        t2 = t  * t 
        t4 = t2 * t2


        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t2*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) ) * ( (yPA * yPB) - t2/(albe+t2) * yPC * (yPA+yPB) + 0.5d0 * 1.d0/(albe+t2) + (t2/(albe+t2)) * (t2/(albe+t2)) * (yPC*yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC) )

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_py_Toroidal

subroutine integrate_NA_py_pz_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 
        double precision             :: t2    , t4
        double precision             :: yPA , zPB
        double precision,parameter   :: pi = 3.14159265358979323846D00


        yPA = yp-ya
        zPB = zp-zb 

        t2 = t  * t 
        t4 = t2 * t2


        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t2*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe)/ax**2 + dx) *  I_0_x
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) ) * (yPA - t2/(albe+t2) * (yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC) ) * (zPB - t2/(albe+t2) * (zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_py_pz_Toroidal

subroutine integrate_NA_pz_px_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 

        double precision             :: t2   , t4
        double precision             :: zPA
        double precision,parameter   :: pi = 3.14159265358979323846D00

        t2 = t  * t
        t4 = t2 * t2

        zPA = zp-za

        xD   = datan((gamma_x*dsin(ax*xp)+t**2.d0*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t**2.d0*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t**2.d0*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t**4 + 2.d0 * gamma_x * t**2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t2+albe)/ax**2 + dx) * I_1_x * dsin(ax*(xD-xB))
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC)) 
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC)) * (zPA - t2/(albe+t2) * (zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_px_Toroidal


subroutine integrate_NA_pz_py_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 
        double precision             :: t2    , t4
        double precision             :: yPB , zPA
        double precision,parameter   :: pi = 3.14159265358979323846D00


        yPB = yp-yb
        zPA = zp-za

        t2 = t  * t 
        t4 = t2 * t2


        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t2*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = 1.d0/(albe+t2) * dexp(-2.d0*(t2+albe)/ax**2 + dx) *  I_0_x
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) ) * (yPB - t2/(albe+t2) * (yPC))
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC) ) * (zPA - t2/(albe+t2) * (zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_py_Toroidal

subroutine integrate_NA_pz_pz_Toroidal(xpC,ypC,zpC, &
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

      double precision,parameter    :: epsabs = 1.0e-16 , epsrel = 1.0e-10
      integer         ,parameter    :: inf   = 1 
      double precision,parameter    :: bound = 0.0d0
      integer         ,parameter    :: limit = 100
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

        use gsl_bessel_mod
        double precision, intent(in) :: t
        double precision             :: fx

        double precision             :: I_0_x
        double precision             :: I_1_x
        double precision             :: I_2_x

        double precision             :: dx
        double precision             :: xD
        double precision             :: Nax   , Nay   , NAz 
        double precision             :: t2    , t4
        double precision             :: zPA , zPB
        double precision,parameter   :: pi = 3.14159265358979323846D00

        zPA = zp-za
        zPB = zp-zb

        t2 = t  * t 
        t4 = t2 * t2


        xD   = datan((gamma_x*dsin(ax*xp)+t2*dsin(ax*xc))/(gamma_x*dcos(ax*xp)+t2*dcos(ax*xc)))/ax + 0.5d0 * Lx * Heaviside(-gamma_x*dcos(ax*xp)-t2*dcos(ax*xc))

        dx = 2.d0*dsqrt( gamma_x**2 + t4 + 2.d0 * gamma_x * t2 * dcos(ax*(XpC))) / ax**2
        
        I_0_x = bessi_scaled(0, dx)
        I_1_x = bessi_scaled(1, dx)
        I_2_x = bessi_scaled(2, dx)

        Nax   = dexp(-2.d0*(t**2+albe)/ax**2 + dx) * I_0_x
        Nay   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (yPC*yPC) ) 
        Naz   = dsqrt(pi/(albe+t2)) * dexp(-(albe*t2)/(albe+t2) * (zPC*zPC) ) * ( (zPA * zPB) - t2/(albe+t2) * yPC * (zPA+zPB) + 0.5d0 * 1.d0/(albe+t2) + (t2/(albe+t2)) * (t2/(albe+t2)) * (zPC*zPC))

        fx    = Nax * Nay * Naz

      end function f_decay
    
end subroutine integrate_NA_pz_pz_Toroidal