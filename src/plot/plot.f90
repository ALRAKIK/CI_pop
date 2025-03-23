subroutine  plot(n_atoms,geometry)

      use atom_basis
      use torus_init
      implicit none 

      integer,intent(in)               :: n_atoms 
      double precision,intent(in)      :: geometry(n_atoms,3)

      call Torus_def()

      
      open(30,file="./tmp/parameter.dat")
        write(30,"(f16.8)")  Lx
        write(30,"(10f16.8)")  geometry(1,1), geometry(2,1) , geometry((n_atoms/2)+1,1) 
      close(30)

      call system("python src/plot/plot_gaussian.py")


end subroutine