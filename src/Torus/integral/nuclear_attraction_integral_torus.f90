subroutine nuclear_attraction_integral_ss_torus(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_ss_normal)

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
      double precision  ,intent(out) :: S_ss_normal

      integer                        :: i , j , k 
      integer                        :: charge_atom
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: p,mu
      double precision               :: Two_PIP , P_R
      double precision               :: x1 , x2 , X
      double precision               :: y1 , y2 , Y 
      double precision               :: z1 , z2 , Z
      double precision               :: xp , yp , zp 
      double precision               :: D_normal 
      double precision               :: xPC , yPC ,  zPC
      double precision               :: Boys_func , R2PC

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      call PBC(x1,x2,X)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!
 
      S_ss_normal = 0.d0

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

            do k = 1 , number_of_atoms
              
              xPC = xp - geometry(k,1) 
              yPC = yp - geometry(k,2)
              zPC = zp - geometry(k,3)

              call euc(xp,geometry(k,1),xPC)
              
              R2PC = xPC*xPC + yPC*yPC + zPC*zPC
            
              charge_atom = (-1)*atoms(k)%charge
            
              S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * Boys_func(0,p*R2PC)
            
            end do 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_ss_torus

subroutine nuclear_attraction_integral_sp_torus(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_sp_normal)

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

end subroutine nuclear_attraction_integral_sp_torus


subroutine nuclear_attraction_integral_pp_torus(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: AO1 , AO2
      type(atom)        ,intent(in)   :: atoms(number_of_atoms)
      integer           ,intent(in)   :: number_of_atoms
      double precision  ,intent(in)   :: r1(3) , r2(3)
      double precision  ,intent(in)   :: geometry(number_of_atoms,3)
      double precision  ,intent(out)  :: S_pp_normal


      integer                       :: i , j , k 
      integer                       :: charge_atom
      double precision,parameter    :: pi = 3.14159265358979323846D00
      double precision              :: alpha , beta
      double precision              :: c1    , c2 
      double precision              :: p,mu
      double precision              :: Two_PIP , P_R
      double precision              :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision              :: xp , yp , zp 
      double precision              :: X , Y , Z
      double precision              :: x_PB_normal , y_PB_normal , z_PB_normal
      double precision              :: x_PA_normal , y_PA_normal , z_PA_normal
      double precision              :: D_normal 
      double precision              :: x_PC , y_PC , z_PC 
      double precision              :: term1 , term2 , term3
      double precision              :: const 
      double precision              :: Boys_func , R2PC
      

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!
 
      S_pp_normal = 0.d0

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

            X_PB_normal  =  (alpha*p_R)*(X)
            Y_PB_normal  =  (alpha*p_R)*(Y)
            Z_PB_normal  =  (alpha*p_R)*(Z)

            X_PA_normal  = -(beta*p_R)*(X)
            Y_PA_normal  = -(beta*p_R)*(Y)
            Z_PA_normal  = -(beta*p_R)*(Z)

            do k = 1 , number_of_atoms
              
              x_PC = xp - geometry(k,1) 
              y_PC = yp - geometry(k,2)
              z_PC = zp - geometry(k,3)

              call euc(xp,geometry(k,1),x_PC)
              
              R2PC = x_PC*x_PC + y_PC*y_PC + z_PC*z_PC
            
              charge_atom = (-1)*atoms(k)%charge

              const = c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal)
            
              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  

                term1 =  ( X_PB_normal *  X_PA_normal + 0.50d0 * P_R )            *   Boys_func(0,p*R2PC)
                term2 = -( X_PC * (X_PA_normal + X_PB_normal ) + 0.50d0 * P_R )   *   Boys_func(1,p*R2PC)
                term3 =  ( X_PC **2 )                                             *   Boys_func(2,p*R2PC)
              
                if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  const * (term1 + term2 + term3) 

                if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  const * ( X_PA_normal * Y_PB_normal * Boys_func(0,p*R2PC) - ( X_PA_normal * Y_PC  + Y_PB_normal * X_PC  ) * Boys_func(1,p*R2PC) +  X_PC * Y_PC *  Boys_func(2,p*R2PC) )      

                if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  const * ( X_PA_normal * Z_PB_normal * Boys_func(0,p*R2PC) - ( X_PA_normal * Z_PC  + Z_PB_normal * X_PC  ) * Boys_func(1,p*R2PC) +  X_PC * Z_PC *  Boys_func(2,p*R2PC) )      

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  

                if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  const * ( Y_PA_normal * X_PB_normal * Boys_func(0,p*R2PC) - ( Y_PA_normal * X_PC  + X_PB_normal * Y_PC  ) * Boys_func(1,p*R2PC) +  Y_PC * X_PC *  Boys_func(2,p*R2PC) )      ! Py-Px

                term1 =  ( Y_PB_normal *  Y_PA_normal + 0.50d0 * P_R )           *   Boys_func(0,p*R2PC)
                term2 = -( Y_PC * (Y_PA_normal + Y_PB_normal ) + 0.50d0 * P_R )  *   Boys_func(1,p*R2PC)
                term3 =  ( Y_PC **2 )                                            *   Boys_func(2,p*R2PC)

                if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  const * (term1 + term2 + term3 )

                if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  const * ( Y_PA_normal * Z_PB_normal * Boys_func(0,p*R2PC) - ( Y_PA_normal * Z_PC  + Z_PB_normal * Y_PC  ) * Boys_func(1,p*R2PC) +  Y_PC * Z_PC *  Boys_func(2,p*R2PC) )      ! Py-Pz

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  
             
                if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  const * ( Z_PA_normal * X_PB_normal * Boys_func(0,p*R2PC) - ( Z_PA_normal * X_PC  + X_PB_normal * Z_PC  ) * Boys_func(1,p*R2PC) +  Z_PC * X_PC *  Boys_func(2,p*R2PC) )      ! Pz-Px

                if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  const * ( Z_PA_normal * Y_PB_normal * Boys_func(0,p*R2PC) - ( Z_PA_normal * Y_PC  + Y_PB_normal * Z_PC  ) * Boys_func(1,p*R2PC) +  Z_PC * Y_PC *  Boys_func(2,p*R2PC) )     ! Pz-Py

                term1 =  ( Z_PB_normal *  Z_PA_normal + 0.50d0 * P_R )            *   Boys_func(0,p*R2PC)
                term2 = -( Z_PC  * (Z_PA_normal + Z_PB_normal ) + 0.50d0 * P_R )  *   Boys_func(1,p*R2PC)
                term3 =  ( Z_PC  **2 )                                            *   Boys_func(2,p*R2PC)

                if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  const * (term1 + term2 + term3 )

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  

            end do 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine nuclear_attraction_integral_pp_torus