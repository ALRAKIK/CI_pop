subroutine ERI_integral(number_of_atoms,geometry,atoms)

      use atom_basis
      use classification_ERI

      implicit none 

      integer                        :: i , j , k , l 
      integer                        :: number_of_atoms
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision               :: value 
      integer                        :: number_of_functions


      number_of_functions = 0 
      do i = 1 , number_of_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

!      call print_orbital_table(ERI,number_of_functions)

      open(1,file="./tmp/ERI.dat")

      do i = 1, number_of_functions
        do j = 1 , number_of_functions
          do k = 1 , number_of_functions
            do l = 1 , number_of_functions
              call ERI_integral_4_function(ERI(i),ERI(j),ERI(k),ERI(l),value)
              if (abs(value) > 1e-10 ) write(1,"(I5,I5,I5,I5,f16.10)") i , j , k , l , value 
            end do 
          end do 
        end do 
      end do 

      close(1)

      deallocate(ERI)

      


end subroutine