module classification_ERI
      
      implicit none

      type :: ERI_function

        double precision               :: x , y , z 
        double precision , allocatable :: coefficient(:)
        double precision , allocatable :: exponent(:)
        character(LEN=2)               :: orbital 

      end type ERI_function


      contains

      subroutine classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)
    
      use atom_basis

      implicit none 

      integer             :: i , j , k
      integer             :: number_of_functions
      integer             :: number_of_atoms 
      double precision    :: geometry(number_of_atoms,3)
      type(atom)          :: atoms(number_of_atoms)
      type(ERI_function)  :: ERI  (number_of_functions)
      
      j = 0 
      
      do k = 1 , number_of_atoms

        do i = 1 , atoms(k)%num_s_function
              
          j= j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "s"
          ERI(j)%exponent    = atoms(k)%exponent_s
          ERI(j)%coefficient = atoms(k)%coefficient_s(:,i)

        end do 

        do i = 1 , atoms(k)%num_p_function

          j= j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "px"
          ERI(j)%exponent    = atoms(k)%exponent_p
          ERI(j)%coefficient = atoms(k)%coefficient_p(:,i)

          j= j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "py"
          ERI(j)%exponent    = atoms(k)%exponent_p
          ERI(j)%coefficient = atoms(k)%coefficient_p(:,i)

          j= j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "pz"
          ERI(j)%exponent    = atoms(k)%exponent_p
          ERI(j)%coefficient = atoms(k)%coefficient_p(:,i)

          end do 

      end do 


      end subroutine


end module 