subroutine ERI_integral(number_of_atoms,geometry,atoms)

      use atom_basis
      use classification_ERI

      implicit none 

      integer                        :: i 
      integer                        :: number_of_atoms
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)

      integer                        :: number_of_functions


      number_of_functions = 0 
      do i = 1 , number_of_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)


      write(*,*) ""
      write(*,*) ""
       
      do i = 1 , number_of_functions
        write(*,'(I10,3f12.8,2x,a,2x,I0)') i , ERI(i)%x , ERI(i)%y , ERI(i)%z , ERI(i)%orbital , ERI(i)%n_contract  
      end do 






end subroutine