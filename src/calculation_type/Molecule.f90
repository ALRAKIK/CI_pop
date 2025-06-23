subroutine molecule(n_atoms,geometry,atoms)

      use files
      use torus_init
      use atom_basis

      implicit none

      ! - input - !
       
      integer          ,intent(in) :: n_atoms
      double precision ,intent(in) :: geometry(n_atoms,3)
      type(atom)       ,intent(in) :: atoms(n_atoms)

      ! - local - !

      double precision             :: start,end,time

      ! - output - ! 

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header(" OBC Calculation ",-1)
      call header_under(" Calculate the integerals ",-1)


      call cpu_time(start)

        call overlap_matrix             (n_atoms,geometry,atoms)
        call kinetic_matrix             (n_atoms,geometry,atoms)
        call nuclear_attraction_matrix  (n_atoms,geometry,atoms)  
        call ERI_integral               (n_atoms,geometry,atoms)

      call cpu_time(end)


      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',time,' seconds'
      write(outfile,*)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine molecule 