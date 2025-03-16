subroutine overlap_integral_ss_alt(r1,r2,AO1,AO2,S_ss_normal)

      use torus_init
      use classification_ERI

      implicit none 

      double precision,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2

      integer                        :: i , j 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: p,mu
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: D_normal 
      double precision,intent(out)   :: S_ss_normal

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
      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)
          p  = alpha + beta 
          mu = alpha*beta/p 
          S_ss_normal =  S_ss_normal +  c1*c2*(dsqrt(pi/p)**3.d0)*exp(-mu*D_normal)
        end do 
      end do
    
      !-----------------------------------------------------------------!
    
end subroutine overlap_integral_ss_alt


subroutine overlap_integral_sp_alt(r1,r2,AO1,AO2,S_sp_normal)

      use torus_init
      use classification_ERI

      implicit none 

      double precision,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2

      integer                        :: i , j 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: p,mu
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: D_normal 
      double precision               :: X_PB_normal    , Y_PB_normal , Z_PB_normal
      double precision,intent(out)   :: S_sp_normal

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
      Y            = (y1 - y2)
      Z            = (z1 - z2)

      D_normal     = (X*X+Y*Y+Z*Z)

      !-----------------------------------------------------------------!

      S_sp_normal = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha =  AO1%exponent(i)
        c1    =  AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta =   AO2%exponent(j)
          c2   =   AO2%coefficient(j)

          p  = alpha + beta 
          mu = alpha*beta/p 

          X            =  (x1 - x2)
          if (torus) call SSD(x1,x2,X)
          X_PB_normal  =  (alpha/p)*(X)

          Y            =  (y1 - y2)
          Y_PB_normal  =  (alpha/p)*(Y)

          Z            =  (z1 - z2)
          Z_PB_normal  =  (alpha/p)*(Z)
          
          if (AO2%orbital=="px") S_sp_normal =  S_sp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal
          if (AO2%orbital=="py") S_sp_normal =  S_sp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal
          if (AO2%orbital=="pz") S_sp_normal =  S_sp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_sp_alt

subroutine overlap_integral_pp_alt(r1,r2,AO1,AO2,S_pp_normal)

      use torus_init
      use classification_ERI

      implicit none 

      double precision,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)  :: AO1 , AO2

      integer                        :: i , j 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      double precision               :: alpha , beta
      double precision               :: c1    , c2 
      double precision               :: p,mu
      double precision               :: x1 , x2 , y1 , y2 , z1 , z2 
      double precision               :: X , Y , Z
      double precision               :: D_normal 
      double precision               :: X_PB_normal    , Y_PB_normal , Z_PB_normal
      double precision               :: X_PA_normal    , Y_PA_normal , Z_PA_normal
      double precision               :: C_X_normal     , C_Y_normal  , C_Z_normal
      double precision,intent(out)   :: S_pp_normal

      x1 = r1(1) ; x2 = r2(1) 
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)

      X            = (x1 - x2)
      if (torus) call PBC(x1,x2,X)
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

          if (AO1%orbital == "px" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_X_normal                   
          if (AO1%orbital == "px" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal*Y_PA_normal
          if (AO1%orbital == "px" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * X_PB_normal*Z_PA_normal

          if (AO1%orbital == "py" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal*X_PA_normal
          if (AO1%orbital == "py" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_Y_normal             
          if (AO1%orbital == "py" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Y_PB_normal*Z_PA_normal

          if (AO1%orbital == "pz" .and. AO2%orbital == "px") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal*X_PA_normal
          if (AO1%orbital == "pz" .and. AO2%orbital == "py") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * Z_PB_normal*Y_PA_normal
          if (AO1%orbital == "pz" .and. AO2%orbital == "pz") S_pp_normal =  S_pp_normal +  c1*c2*(dsqrt(pi/p)**3)*exp(-mu*D_normal) * C_Z_normal             

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_pp_alt