subroutine nuclear_attraction_integral_ss(n_atoms,geometry,r1,r2,atom1,atom2,atoms,index1,index2,S_ss_normal)

      use torus_init
      use atom_basis

      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      type(atom)                   :: atoms(n_atoms)
      integer,intent(in)           :: n_atoms
      double precision,intent(in)  :: geometry(n_atoms,3)
      integer                      :: index1 , index2 

      integer                      :: i , j  , k 
      integer                      :: charge_atom
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: Two_PIP , P_R
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: xp , yp , zp 
      double precision             :: X , Y , Z
      double precision             :: D_normal 
      double precision             :: xPC , yPC , zPC
      double precision             :: Boys_func , R2PC
      double precision,intent(out) :: S_ss_normal

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

  !-----------------------------------------------------------------!

      S_ss_normal = 0.d0

      do i = 1 , atom1%num_exponent_s
        alpha = atom1%exponent_s(i)
        c1    = atom1%coefficient_s(i,index1)
        do j = 1 , atom2%num_exponent_s
          beta = atom2%exponent_s(j)
          c2   = atom2%coefficient_s(j,index2)
              p       = alpha + beta 
              P_R     = 1.d0 / p
              mu      = alpha*beta*p_R
              Two_PIP = 2.0d0*pi*p_R

              xp = (alpha*x1+beta*x2)*p_R
              if (torus) call bary_center(alpha,x1,beta,x2,p_R,xp)
              yp = (alpha*y1+beta*y2)*p_R
              zp = (alpha*z1+beta*z2)*p_R
        
              do k = 1 , n_atoms
              
                xPC = xp - geometry(k,1) 
                yPC = yp - geometry(k,2)
                zPC = zp - geometry(k,3)

                if (torus) call euc(xp,geometry(k,1),xPC)
                
                R2PC = xPC*xPC + yPC*yPC + zPC*zPC
              
                charge_atom = (-1)*atoms(k)%charge
              
                S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * Boys_func(0,p*R2PC)
              
              end do 
            
        end do 
      end do

!-----------------------------------------------------------------!

end subroutine

subroutine nuclear_attraction_integral_sp(n_atoms,geometry,r1,r2,atom1,atom2,atoms,index1,index2,S_sp_normal)

      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      type(atom)                   :: atoms(n_atoms)
      integer,intent(in)           :: n_atoms
      double precision,intent(in)  :: geometry(n_atoms,3)
      integer                      :: index1 , index2 

      integer                      :: i , j  , k 
      integer                      :: charge_atom
      double precision,parameter   :: pi = dacos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: Two_PIP , P_R
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: xp , yp , zp 
      double precision             :: X , Y , Z
      double precision             :: X_PB_normal , Y_PB_normal , Z_PB_normal
      double precision             :: D_normal 
      double precision             :: X_PC_normal , Y_PC_normal , Z_PC_normal
      double precision             :: Boys_func , R2PC
      double precision,intent(out) :: S_sp_normal(3)

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!

      S_sp_normal(:) = 0.d0

      do i = 1 , atom1%num_exponent_s
        alpha = atom1%exponent_s(i)
        c1    = atom1%coefficient_s(i,index1)
        do j = 1 , atom2%num_exponent_p
          beta = atom2%exponent_p(j)
          c2   = atom2%coefficient_p(j,index2)
              p       = alpha + beta 
              P_R     = 1.d0 / p
              mu      = alpha*beta*p_R
              Two_PIP = 2.0d0*pi*p_R

              xp = (alpha*x1+beta*x2)*p_R
              if (torus) call bary_center(alpha,x1,beta,x2,p_R,xp)
              yp = (alpha*y1+beta*y2)*p_R
              zp = (alpha*z1+beta*z2)*p_R
        
              X = x1 - x2 
              if (torus) call SSD(x1,x2,X)
              X_PB_normal  =  (alpha/p)*(X)
              Y_PB_normal  =  (alpha/p)*(Y)
              Z_PB_normal  =  (alpha/p)*(Z)
        
              do k = 1 , n_atoms
              
                X_PC_normal = xp - geometry(k,1) 
                Y_PC_normal = yp - geometry(k,2)
                Z_PC_normal = zp - geometry(k,3)

                if (torus) call euc(xp,geometry(k,1),X_PC_normal)

                R2PC = X_PC_normal*X_PC_normal + Y_PC_normal*Y_PC_normal + Z_PC_normal*Z_PC_normal
              
                charge_atom = (-1)*atoms(k)%charge
              
                S_sp_normal(1) =  S_sp_normal(1) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (X_PB_normal * Boys_func(0,p*R2PC) - X_PC_normal * Boys_func(1,p*R2PC) ) 
                S_sp_normal(2) =  S_sp_normal(2) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Y_PB_normal * Boys_func(0,p*R2PC) - Y_PC_normal * Boys_func(1,p*R2PC) ) 
                S_sp_normal(3) =  S_sp_normal(3) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Z_PB_normal * Boys_func(0,p*R2PC) - Z_PC_normal * Boys_func(1,p*R2PC) ) 

              end do 
            
        end do 
      end do


!-----------------------------------------------------------------!

end subroutine

subroutine nuclear_attraction_integral_ps(n_atoms,geometry,r1,r2,atom1,atom2,atoms,index1,index2,S_ps_normal)
  
      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      type(atom)                   :: atoms(n_atoms)
      integer,intent(in)           :: n_atoms
      double precision,intent(in)  :: geometry(n_atoms,3)
      integer                      :: index1 , index2 

      integer                      :: i , j  , k 
      integer                      :: charge_atom
      double precision,parameter   :: pi = dacos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: Two_PIP , P_R
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: xp , yp , zp 
      double precision             :: X , Y , Z
      double precision             :: X_PA_normal , Y_PA_normal , Z_PA_normal
      double precision             :: D_normal 
      double precision             :: X_PC_normal , Y_PC_normal , Z_PC_normal 
      double precision             :: Boys_func , R2PC
      double precision,intent(out) :: S_ps_normal(3)

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!

      S_ps_normal(:) = 0.d0

      do i = 1 , atom1%num_exponent_p
        alpha = atom1%exponent_p(i)
        c1    = atom1%coefficient_p(i,index1)
        do j = 1 , atom2%num_exponent_s
          beta = atom2%exponent_s(j)
          c2   = atom2%coefficient_s(j,index2)
              p       = alpha + beta 
              P_R     = 1.d0 / p
              mu      = alpha*beta*p_R
              Two_PIP = 2.0d0*pi*p_R
        
              xp = (alpha*x1+beta*x2)*p_R
              if (torus) call bary_center(alpha,x1,beta,x2,p_R,xp)
              yp = (alpha*y1+beta*y2)*p_R
              zp = (alpha*z1+beta*z2)*p_R
        
              X = x1 - x2 
              if (torus) call SSD(x1,x2,X)
              X_PA_normal  = -(beta*p_R)*(X)
              Y_PA_normal  = -(beta*p_R)*(Y)
              Z_PA_normal  = -(beta*p_R)*(Z)
        
              do k = 1 , n_atoms
              
                X_PC_normal = xp - geometry(k,1) 
                Y_PC_normal = yp - geometry(k,2)
                Z_PC_normal = zp - geometry(k,3)

                if (torus) call euc(xp,geometry(k,1),X_PC_normal)
              
                R2PC = X_PC_normal*X_PC_normal + Y_PC_normal*Y_PC_normal + Z_PC_normal*Z_PC_normal
              
                charge_atom = (-1)*atoms(k)%charge
              
                S_ps_normal(1) =  S_ps_normal(1) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (X_PA_normal * Boys_func(0,p*R2PC) - X_PC_normal * Boys_func(1,p*R2PC) ) 
                S_ps_normal(2) =  S_ps_normal(2) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Y_PA_normal * Boys_func(0,p*R2PC) - Y_PC_normal * Boys_func(1,p*R2PC) ) 
                S_ps_normal(3) =  S_ps_normal(3) +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * (Z_PA_normal * Boys_func(0,p*R2PC) - Z_PC_normal * Boys_func(1,p*R2PC) ) 
              
              end do 
            
        end do 
      end do


!-----------------------------------------------------------------!

end subroutine


subroutine nuclear_attraction_integral_pp(n_atoms,geometry,r1,r2,atom1,atom2,atoms,index1,index2,S_pp_normal)
  
      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      type(atom)                   :: atoms(n_atoms)
      integer,intent(in)           :: n_atoms
      double precision,intent(in)  :: geometry(n_atoms,3)
      integer                      :: index1 , index2 

      integer                      :: i , j  , k 
      integer                      :: charge_atom
      double precision,parameter   :: pi = dacos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: Two_PIP , P_R
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: xp , yp , zp 
      double precision             :: X , Y , Z
      double precision             :: X_PA_normal , Y_PA_normal , Z_PA_normal
      double precision             :: X_PB_normal , Y_PB_normal , Z_PB_normal
      double precision             :: D_normal 
      double precision             :: X_PC_normal , Y_PC_normal , Z_PC_normal 
      double precision             :: term1 , term2 , term3
      double precision             :: const 
      double precision             :: Boys_func , R2PC
      double precision,intent(out) :: S_pp_normal(3,3)

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!

      S_pp_normal(:,:) = 0.d0

      do i = 1 , atom1%num_exponent_p
        alpha = atom1%exponent_p(i)
        c1    = atom1%coefficient_p(i,index1)
        do j = 1 , atom2%num_exponent_p
          beta = atom2%exponent_p(j)
          c2   = atom2%coefficient_p(j,index2)
              p       = alpha + beta 
              P_R     = 1.d0 / p
              mu      = alpha*beta*p_R
              Two_PIP = 2.0d0*pi*p_R
        
              xp = (alpha*x1+beta*x2)*p_R
              if (torus) call bary_center(alpha,x1,beta,x2,p_R,xp)
              yp = (alpha*y1+beta*y2)*p_R
              zp = (alpha*z1+beta*z2)*p_R
        
              X = x1 - x2 
              if (torus) call SSD(x1,x2,X)
              X_PB_normal  =  (alpha*p_R)*(X)
              Y_PB_normal  =  (alpha*p_R)*(Y)
              Z_PB_normal  =  (alpha*p_R)*(Z)

              X_PA_normal  = -(beta*p_R)*(X)
              Y_PA_normal  = -(beta*p_R)*(Y)
              Z_PA_normal  = -(beta*p_R)*(Z)
        
              do k = 1 , n_atoms
              
                X_PC_normal = xp - geometry(k,1) 
                Y_PC_normal = yp - geometry(k,2)
                Z_PC_normal = zp - geometry(k,3)

                if (torus) call euc(xp,geometry(k,1),X_PC_normal)
              
                R2PC = X_PC_normal*X_PC_normal + Y_PC_normal*Y_PC_normal + Z_PC_normal*Z_PC_normal
              
                charge_atom = (-1)*atoms(k)%charge

                const = c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal)

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  

                term1 =  ( X_PB_normal *  X_PA_normal + 0.50d0 * P_R )                  *   Boys_func(0,p*R2PC)
                term2 = -( X_PC_normal * (X_PA_normal + X_PB_normal ) + 0.50d0 * P_R )  *   Boys_func(1,p*R2PC)
                term3 =  ( X_PC_normal**2 )                                             *   Boys_func(2,p*R2PC)
              
                S_pp_normal(1,1) =  S_pp_normal(1,1) +  const * (term1 + term2 + term3) 

                S_pp_normal(1,2) =  S_pp_normal(1,2) +  const * ( X_PA_normal * Y_PB_normal * Boys_func(0,p*R2PC) - ( X_PA_normal * Y_PC_normal  + Y_PB_normal * X_PC_normal  ) * Boys_func(1,p*R2PC) +  X_PC_normal * Y_PC_normal *  Boys_func(2,p*R2PC) )      ! Px-Py

                S_pp_normal(1,3) =  S_pp_normal(1,3) +  const * ( X_PA_normal * Z_PB_normal * Boys_func(0,p*R2PC) - ( X_PA_normal * Z_PC_normal  + Z_PB_normal * X_PC_normal  ) * Boys_func(1,p*R2PC) +  X_PC_normal * Z_PC_normal *  Boys_func(2,p*R2PC) )      ! Px-Pz

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  

                S_pp_normal(2,1) =  S_pp_normal(2,1) +  const * ( Y_PA_normal * X_PB_normal * Boys_func(0,p*R2PC) - ( Y_PA_normal * X_PC_normal  + X_PB_normal * Y_PC_normal  ) * Boys_func(1,p*R2PC) +  Y_PC_normal * X_PC_normal *  Boys_func(2,p*R2PC) )      ! Py-Px

                term1 =  ( Y_PB_normal *  Y_PA_normal + 0.50d0 * P_R )                  *   Boys_func(0,p*R2PC)
                term2 = -( Y_PC_normal * (Y_PA_normal + Y_PB_normal ) + 0.50d0 * P_R )  *   Boys_func(1,p*R2PC)
                term3 =  ( Y_PC_normal**2 )                                             *   Boys_func(2,p*R2PC)

                S_pp_normal(2,2) =  S_pp_normal(2,2) +  const * (term1 + term2 + term3 )

                S_pp_normal(2,3) =  S_pp_normal(2,3) +  const * ( Y_PA_normal * Z_PB_normal * Boys_func(0,p*R2PC) - ( Y_PA_normal * Z_PC_normal  + Z_PB_normal * Y_PC_normal  ) * Boys_func(1,p*R2PC) +  Y_PC_normal * Z_PC_normal *  Boys_func(2,p*R2PC) )      ! Py-Pz

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  
             
                S_pp_normal(3,1) =  S_pp_normal(3,1) +  const * ( Z_PA_normal * X_PB_normal * Boys_func(0,p*R2PC) - ( Z_PA_normal * X_PC_normal  + X_PB_normal * Z_PC_normal  ) * Boys_func(1,p*R2PC) +  Z_PC_normal * X_PC_normal *  Boys_func(2,p*R2PC) )      ! Pz-Px

                S_pp_normal(3,2) =  S_pp_normal(3,2) +  const * ( Z_PA_normal * Y_PB_normal * Boys_func(0,p*R2PC) - ( Z_PA_normal * Y_PC_normal  + Y_PB_normal * Z_PC_normal  ) * Boys_func(1,p*R2PC) +  Z_PC_normal * Y_PC_normal *  Boys_func(2,p*R2PC) )     ! Pz-Py

                term1 =  ( Z_PB_normal *  Z_PA_normal + 0.50d0 * P_R )                  *   Boys_func(0,p*R2PC)
                term2 = -( Z_PC_normal * (Z_PA_normal + Z_PB_normal ) + 0.50d0 * P_R )  *   Boys_func(1,p*R2PC)
                term3 =  ( Z_PC_normal**2 )                                             *   Boys_func(2,p*R2PC)

                S_pp_normal(3,3) =  S_pp_normal(3,3) +  const * (term1 + term2 + term3 )

              ! / ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- / !  


              end do 
            
        end do 
      end do


!-----------------------------------------------------------------!

end subroutine