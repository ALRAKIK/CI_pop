subroutine basis_info(number_of_functions,number_of_primitives,number_of_shells,n_atom_unitcell)

      use files

      implicit none 

      integer     :: number_of_functions,number_of_primitives
      integer     :: number_of_shells,n_atom_unitcell


      ! --------------------------------------------------------------- !

      write(outfile,'(A,I10)')'The number of AO functions            :',&
      &                        number_of_functions 

      write(outfile,'(A)')

      write(outfile,'(A,I10)')'The number of the Gaussian primitives :',&
      &                        number_of_primitives

      write(outfile,'(A)')

      write(outfile,'(A,I10)')'The number of shells                  :',&
      &                        number_of_shells

      write(outfile,'(A)')

      write(outfile,'(A,I10)')'The number of atoms in the unitcell   :',&
      &                     n_atom_unitcell

      write(outfile,'(A)')

      flush(outfile)

      ! --------------------------------------------------------------- !

end subroutine basis_info