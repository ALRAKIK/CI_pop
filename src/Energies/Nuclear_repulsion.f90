subroutine NRE(number_of_atoms,geometry,atoms,E)

      use torus_init
      use atom_basis 

      implicit none 

      integer         ,intent(in)    :: number_of_atoms
      double precision,intent(in)    :: geometry(number_of_atoms,3)
      type(atom)      ,intent(in)    :: atoms(number_of_atoms)

      double precision,intent(out)   :: E 

      integer                        :: i , j 
      double precision               :: x , y , z , dist 


      E = 0.d0 

        do i = 1 , number_of_atoms-1 
          do j = i+1 , number_of_atoms
            x  = geometry(i,1) - geometry(j,1)
            if (torus) x  = (dsqrt(2.d0 * (1.d0 - dcos(ax*x)) / (ax * ax)))
            y  = geometry(i,2) - geometry(j,2)
            z  = geometry(i,3) - geometry(j,3)
            dist =  x*x+y*y+z*z
            dist = dsqrt(dist)
            E = E + atoms(i)%charge*atoms(j)%charge/dist
          end do 
        end do 

        if (E > 1e18) then 
          call header("Error",-1)
          write(*,'(a80)')'*******************************************************************************************'
          write(*,'(a80)')           "* Nuclear repulsion energy is too large, Two atoms are too close to each other *"
          write(*,'(a,a)')        "*","                                                                              *"
          write(*,'(a,f13.6,a)')             "* E_NUC = ",E,"                                                        *"
          write(*,'(a,a)')        "*","                                                                              *"
          write(*,'(a80)')'*******************************************************************************************'
          stop 
        end if

end subroutine NRE 