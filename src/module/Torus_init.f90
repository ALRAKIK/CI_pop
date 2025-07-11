module torus_init

      implicit none

      Logical                      :: Torus = .FALSE. 
      double precision             :: Lx = 0.d0 , Ly = 0.d0 , Lz = 0.d0
      double precision             :: ax = 0.d0 , ay = 0.d0 , az = 0.d0
      integer                      :: number_of_atom_in_unitcell

      contains 

      subroutine Torus_def()

        implicit none

        double precision,parameter   :: pi = 3.14159265358979323846D00

        open(50,file="torus_parameters.inp")
          read(50,*) Lx 
          read(50,*)    number_of_atom_in_unitcell
        close(50)
        
        Torus = .TRUE. 

        Ly = Lx ; Lz = Lx 

        ax    = 2*pi/Lx
        ay    = 2*pi/Ly
        az    = 2*pi/Lz

      end subroutine Torus_def 
      

end module torus_init