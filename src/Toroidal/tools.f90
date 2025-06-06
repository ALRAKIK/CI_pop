subroutine bary_center_toroidal(e1,e2,r1,r2,rp)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: r1, r2
      double precision, intent(out) :: rp

      ! ----------------------------------------------------------------!

      rp = datan((e1*dsin(ax*r1)+e2*dsin(ax*r2))/(e1*dcos(ax*r1)+e2*dcos(ax*r2)))/ax + 0.5d0 * Lx * Heaviside(-e1*dcos(ax*r1)-e2*dcos(ax*r2)) 

      if (abs(e1 - e2) < 1.d-10) then 

        if (dabs(r1-r2) > 0.50d0*Lx) then 
          
          if (r1 > r2) then 
            rp =  r1 + (r2 + Lx)
          elseif (r2 > r1) then 
            rp = (r1 + Lx) + r2 
          else 
            rp = r1 + r2  
          end if

          rp = 0.5d0 * rp
        
        else 

          rp = 0.5d0 * ( r1 + r2 )

        end if 

      end if

      if (abs(r1 - r2) < 1e-9) then 
        rp = 0.5d0 * abs((r1 + r2))
      end if







end subroutine bary_center_toroidal