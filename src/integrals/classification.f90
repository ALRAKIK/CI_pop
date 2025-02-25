module classification_ERI
      
      implicit none

      type :: ERI_function

        integer                        :: n_contract 
        double precision               :: x , y , z 
        double precision , allocatable :: coefficient(:,:)
        double precision , allocatable :: exponent(:)
        character(LEN=2)               :: orbital 

      end type ERI_function


      contains

      subroutine classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)
    
      use atom_basis

      implicit none 

      integer             :: i , j , k 
      integer             :: index 
      integer             :: index_s , index_p
      integer             :: number_of_functions
      integer             :: number_of_atoms 
      double precision    :: geometry(number_of_atoms,3)
      type(atom)          :: atoms(number_of_atoms)
      type(ERI_function)  :: ERI  (number_of_functions)
      

      j = 0 

      do i = 1 , number_of_atoms

        index_s =     atoms(i)%num_s_function
        index_p = 3 * atoms(i)%num_p_function   
        index   =     index_s + index_p
    
        ! Assign 's' orbitals first

        do k = 1 , index_s 
            j = j + 1 
            ERI(j)%x = geometry(i,1) 
            ERI(j)%y = geometry(i,2) 
            ERI(j)%z = geometry(i,3)
            ERI(j)%orbital = 's'
            ERI(j)%n_contract   = atoms(i)%num_s_function
            ERI(j)%exponent     = atoms(i)%exponent_s
            ERI(j)%coefficient  = atoms(i)%coefficient_s
        end do 
    
        ! Assign 'p' orbitals in order p_x, p_y, p_z

        do k = 1 , index_p / 3

            j = j + 1 
            ERI(j)%x = geometry(i,1) 
            ERI(j)%y = geometry(i,2) 
            ERI(j)%z = geometry(i,3)
            ERI(j)%orbital = 'px'
            ERI(j)%n_contract   = atoms(i)%num_p_function
            ERI(j)%exponent     = atoms(i)%exponent_p
            ERI(j)%coefficient  = atoms(i)%coefficient_p
    
            j = j + 1 
            ERI(j)%x = geometry(i,1) 
            ERI(j)%y = geometry(i,2) 
            ERI(j)%z = geometry(i,3)
            ERI(j)%orbital = 'py'
            ERI(j)%n_contract   = atoms(i)%num_p_function
            ERI(j)%exponent     = atoms(i)%exponent_p
            ERI(j)%coefficient  = atoms(i)%coefficient_p
    
            j = j + 1 
            ERI(j)%x = geometry(i,1) 
            ERI(j)%y = geometry(i,2) 
            ERI(j)%z = geometry(i,3)
            ERI(j)%orbital = 'pz'
            ERI(j)%n_contract   = atoms(i)%num_p_function
            ERI(j)%exponent     = atoms(i)%exponent_p
            ERI(j)%coefficient  = atoms(i)%coefficient_p

        end do 
      end do 

      end subroutine


end module 