module torus_init

      implicit none

      Logical                      :: Torus = .FALSE. 
      double precision             :: Lx = 0.d0 , Ly = 0.d0 , Lz = 0.d0
      double precision             :: ax = 0.d0 , ay = 0.d0 , az = 0.d0

      contains 

      subroutine Torus_def()

        implicit none

        double precision,parameter   :: pi = 3.14159265358979323846D00

        Lx    = 24.0d0
        Ly    = 24.0d0
        Lz    = 24.0d0

        Torus = .TRUE. 

        ax    = 2*pi/Lx
        ay    = 2*pi/Ly
        az    = 2*pi/Lz

      end subroutine Torus_def 
      

end module torus_init