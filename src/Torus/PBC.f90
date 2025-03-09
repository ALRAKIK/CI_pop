subroutine PBC(x1,x2,X)

      use torus_init
      implicit none 

      double precision,intent(in)   :: x1 , x2 
      double precision,intent(out)  :: X

      X = x1 - x2                       
              
      if (dabs(X) > 0.5d0*Lx) then
        X = Lx - dabs(X)
      end if

end subroutine PBC

subroutine bary_center(alpha,x1,beta,x2,p_R,xp)

      use torus_init
      implicit none

      double precision,intent(in)  :: alpha , x1 , beta , x2 , p_R
      double precision,intent(out) :: xp

      double precision             :: X 

      X = x1 - x2 

      if (dabs(X) > 0.50d0*Lx) then 

        if (x1 > x2) then 
          xp = alpha * x1 + beta * (x2 + Lx) 
        else if (x2 > x1) then 
          xp = alpha * (x1 + Lx) + beta * x2 
        else 
          xp = alpha * x1 + beta * x2  
        end if 

      else
        xp = alpha * x1 + beta * x2
      end if 

      xp = xp * p_R

end subroutine bary_center

subroutine SSD(x1,x2,X)

      use torus_init
      implicit none

      double precision, intent(in)  :: x1 , x2         ! the coordinates on X 
      double precision, intent(out) :: X               ! the distance

      ! local ! 

      double precision      :: sign = 0.d0

      X = x1 - x2 

      if (X >= 0.d0) then 
        sign = +1.d0 
      else
        sign = -1.d0  
      end if 
  
      if (dabs(X) > 0.5d0*Lx) then 
        X = (-1.d0)*sign*(Lx-dabs(X))
      end if 

end subroutine SSD

subroutine euc(x1,x2,X)

      use torus_init
      implicit none 

      double precision,intent(in)   :: x1 , x2 
      double precision,intent(out)  :: X

      double precision              :: sign = 0.d0 
              

      call SSD(x1,x2,X)
    
      if (X >= 0.d0) then 
        sign = +1.d0 
      else
        sign = -1.d0  
      end if 

      X = sign*dabs(X)

      X = sign * (dsqrt(2.d0 * (1.d0 - dcos(ax*x)) / (ax * ax)))

end subroutine euc