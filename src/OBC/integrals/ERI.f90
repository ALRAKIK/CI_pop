subroutine ERI_integral(number_of_atoms,geometry,atoms)

      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      type(atom),intent(in)          :: atoms(number_of_atoms)
      integer,intent(in)             :: number_of_atoms


      type(ERI_function),allocatable :: ERI  (:)
      integer                        :: i , j , k , l 
      integer                        :: number_of_functions
      double precision               :: geometry(number_of_atoms,3)
      double precision               :: value 
      
      !-----------------------------------------------------------------!

      number_of_functions = 0 
      do i = 1 , number_of_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

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