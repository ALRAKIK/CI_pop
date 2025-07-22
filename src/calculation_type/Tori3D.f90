subroutine Tori3D(n_atoms,number_of_functions,atoms,AO,geometry)

      use files
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none

      ! - input - !
       
      integer           ,intent(in)  :: n_atoms , number_of_functions
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)
      type(ERI_function),intent(in)  :: AO (number_of_functions)

      ! - local - !

      double precision             :: start,end,time

      ! - output - ! 

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


      call header("Toroidal  3D Real Gaussian",-1)
      call header_under("Calculate the integerals",-1)
      call Torus_def()
      write(outfile,*) ""
      write(outfile,'(a,f12.8,a,f12.8,a,f12.8,a)') "The length of the box:  Lx , Ly , Lz  = (", Lx ," , ", Ly ," , ",Lz, " )"  
      write(outfile,*) ""


      call cpu_time(start)

        call overlap_matrix_toroidal_3D(n_atoms,number_of_functions,atoms,AO)
        call kinetic_matrix_toroidal_3D(n_atoms,number_of_functions,atoms,AO)
        call nuclear_attraction_matrix_toroidal_3D(n_atoms,number_of_functions,geometry,atoms,AO)
        call ERI_integral_toroidal_3D(n_atoms,geometry,atoms)
                
      call cpu_time(end)

      write(outfile,*) ""

      time = end - start

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',time,' seconds'
      write(outfile,*)

      call system("rm torus_parameters.inp")

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !


end subroutine Tori3D