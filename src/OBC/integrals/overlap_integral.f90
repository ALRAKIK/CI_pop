subroutine overlap_integral_ss(r1,r2,atom1,atom2,index1,index2,S_ss_normal)

      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = 3.14159265358979323846D00
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: D_normal 
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
            if (c1*c2 == 0.d0) cycle  
              p  = alpha + beta 
              mu = alpha*beta/p 
              S_ss_normal =  S_ss_normal +  c1*c2*(dsqrt(pi/p)**3.d0)*exp(-mu*D_normal)
        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine



subroutine overlap_integral_sp(r1,r2,atom1,atom2,index1,index2,S_sp_normal)

      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = dacos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: D_normal  
      double precision             :: X_PB_normal    , Y_PB_normal , Z_PB_normal 
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
            if (c1*c2 == 0.d0) cycle  
              p  = alpha + beta 
              mu = alpha*beta/p 

              X            =  (x1 - x2)
              if (torus) call SSD(x1,x2,X)
              X_PB_normal  =  (alpha/p)*(X)

              Y            =  (y1 - y2)
              Y_PB_normal  =  (alpha/p)*(Y)

              Z            =  (z1 - z2)
              Z_PB_normal  =  (alpha/p)*(Z)


              S_sp_normal(1) =  S_sp_normal(1) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal
              S_sp_normal(2) =  S_sp_normal(2) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal
              S_sp_normal(3) =  S_sp_normal(3) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal


        end do 
      end do

!-----------------------------------------------------------------!

end subroutine


subroutine overlap_integral_ps(r1,r2,atom1,atom2,index1,index2,S_ps_normal)

      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = dacos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: D_normal  
      double precision             :: X_PA_normal    , Y_PA_normal , Z_PA_normal 
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
            if (c1*c2 == 0.d0) cycle  
              p  = alpha + beta 
              mu = alpha*beta/p 
        
              X            =  (x1 - x2)
              if (torus) call SSD(x1,x2,X)
              X_PA_normal  = -(beta/p)*(X)
        
              Y            =  (y1 - y2)
              Y_PA_normal  = -(beta/p)*(Y)
        
              Z            =  (z1 - z2)
              Z_PA_normal  = -(beta/p)*(Z)
        
        
              S_ps_normal(1) =  S_ps_normal(1) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PA_normal
              S_ps_normal(2) =  S_ps_normal(2) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PA_normal
              S_ps_normal(3) =  S_ps_normal(3) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PA_normal
        
        
        end do 
      end do

!-----------------------------------------------------------------!

end subroutine


subroutine overlap_integral_pp(r1,r2,atom1,atom2,index1,index2,S_pp_normal)

      use torus_init
      use atom_basis
      implicit none 

      double precision,intent(in)  :: r1(3) , r2(3)
      type(atom),intent(in)        :: atom1 , atom2 
      integer                      :: index1 , index2 

      integer                      :: i , j 
      double precision,parameter   :: pi = Acos(-1.d0)
      double precision             :: alpha , beta
      double precision             :: c1    , c2 
      double precision             :: p,mu
      double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision             :: X , Y , Z
      double precision             :: D_normal  
      double precision             :: X_PB_normal , Y_PB_normal , Z_PB_normal 
      double precision             :: X_PA_normal , Y_PA_normal , Z_PA_normal 
      double precision             :: C_X_normal  , C_Y_normal  , C_Z_normal 
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
            if (c1*c2 == 0.d0) cycle  
              p  = alpha + beta 
              mu = alpha*beta/p 
        
              X            =  (x1 - x2)
              if (torus) call SSD(x1,x2,X)
              X_PB_normal  =  (alpha/p)*(X)
              X_PA_normal  = -(beta/p) *(X)
              C_X_normal   = X_PB_normal*X_PA_normal+(1/(2.d0*p))
        
              Y            =  (y1 - y2)
              Y_PB_normal  =  (alpha/p)*(Y)
              Y_PA_normal  = -(beta/p) *(Y)
              C_Y_normal   = Y_PB_normal*Y_PA_normal+(1/(2.d0*p))
        
              Z            =  (z1 - z2)
              Z_PB_normal  =  (alpha/p)*(Z)
              Z_PA_normal  = -(beta/p) *(Z)
              C_Z_normal   = Z_PB_normal*Z_PA_normal+(1/(2.d0*p))
        
        
              S_pp_normal(1,1) =  S_pp_normal(1,1) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_X_normal                    ! Px-Px      
              S_pp_normal(1,2) =  S_pp_normal(1,2) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal*Y_PA_normal       ! Px-Py
              S_pp_normal(1,3) =  S_pp_normal(1,3) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal*Z_PA_normal       ! Px-Pz

              S_pp_normal(2,1) =  S_pp_normal(2,1) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal*X_PA_normal       ! Py-Px
              S_pp_normal(2,2) =  S_pp_normal(2,2) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_Y_normal                    ! Py-Py
              S_pp_normal(2,3) =  S_pp_normal(2,3) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal*Z_PA_normal       ! Py-Pz

              S_pp_normal(3,1) =  S_pp_normal(3,1) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal*X_PA_normal       ! Pz-Px
              S_pp_normal(3,2) =  S_pp_normal(3,2) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal*Y_PA_normal       ! Pz-Py
              S_pp_normal(3,3) =  S_pp_normal(3,3) +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_Z_normal                    ! Pz-Pz
        
        
        end do 
      end do

!-----------------------------------------------------------------!

end subroutine