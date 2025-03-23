subroutine ERI_integral_torus(number_of_atoms,geometry,atoms)

      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l
      integer                        :: number_of_atoms
      integer                        :: number_of_atoms_per_unitcell 
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision,allocatable   ::      two_eri(:,:,:,:)
      double precision               :: value
      integer                        :: number_of_functions
      integer                        :: number_of_functions_per_unitcell 

      !-----------------------------------------------------------------!

      number_of_functions = 0 
      do i = 1 , number_of_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      number_of_functions_per_unitcell = 0 
      do i = 1 , number_of_atom_in_unitcell
        number_of_functions_per_unitcell = number_of_functions_per_unitcell + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(number_of_functions))
      allocate(two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions))
      allocate(     two_eri(number_of_functions,number_of_functions,number_of_functions,number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)


      open(1,file="./tmp/ERI.dat")

      ! 2-fold symmetry implementation (k,l permutation only)

      do i = 1, number_of_functions_per_unitcell
!      do i = 1, number_of_functions
        do j = 1, number_of_functions
            do k = 1, number_of_functions
                ! Only calculate for k ≤ l to avoid duplicates
                do l = k, number_of_functions
                    ! Calculate integral once
                    call ERI_integral_4_function_torus(ERI(i),ERI(j),ERI(k),ERI(l), value)
                    
                    ! Store in original position
                    two_electron(i,j,k,l) = value
                    
                    ! Apply symmetry for k!=l
                    two_electron(i,j,l,k) = value
                end do
            end do
        end do
      end do
      
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      call shift_integrals(two_electron,two_eri,number_of_functions,number_of_functions_per_unitcell)

      do i = 1, number_of_functions
        do j = 1 , number_of_functions
          do k = 1 , number_of_functions
            do l = 1 , number_of_functions
              if (abs(two_eri(i,j,k,l)) > 1e-10 ) write(1,"(I5,I5,I5,I5,f16.10)") i , j , k , l , two_eri(i,j,k,l)
            end do 
          end do 
        end do 
      end do 

      close(1)

      deallocate(ERI)
      deallocate(two_electron)
      deallocate(two_eri)

end subroutine ERI_integral_torus