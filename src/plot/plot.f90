subroutine  plot(n_atoms,geometry)

      use atom_basis
      use torus_init
      implicit none 

      integer,intent(in)               :: n_atoms 
      double precision,intent(in)      :: geometry(n_atoms,3)



      open(1,file="./tmp/parameter.dat")
        write(1,"(f16.8)")  Lx
        write(1,"(10f16.8)")  geometry(1,1), geometry(2,1) , geometry((n_atoms/2)+1,1) 
      close(1)

      call system("python src/plot/plot_gaussian.py")


end subroutine