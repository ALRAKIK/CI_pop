subroutine nuclear_attraction_matrix_toroidal(number_of_atoms,number_of_functions,geometry,atoms,AO,NA)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use omp_lib

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l
      integer                        :: index_atom1 , index_unitcell 
      integer                        :: number_of_atoms
      integer                        :: number_of_functions 

      type(atom)                     :: atoms(number_of_atoms)

      type(ERI_function)             :: AO (number_of_functions)
      type(ERI_function)             :: AO1 , AO2

      double precision               :: geometry(number_of_atoms,3)

      double precision               :: r1(3) , r2(3)

      double precision,parameter     :: pi = dacos(-1.d0)

      ! output ! 

      double precision,intent(out)   :: NA(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!
      ! Parallelization variables
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      NA(:,:) = 0.d0 
    
      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function

      index_unitcell = 0 

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!
      ! Setup OpenMP
      !-----------------------------------------------------------------!
      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Nuclear Attraction: Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel
      
      !$omp parallel
        !$omp single
          num_threads = omp_get_num_threads()
        !$omp end single
      !$omp end parallel
      
      ! Calculate optimal chunk size based on available threads
      if (num_threads <= 16) then
        optimal_chunk_size = 16
      else if (num_threads <= 64) then
        optimal_chunk_size = 8
      else 
        optimal_chunk_size = 1
      end if

      !-----------------------------------------------------------------!
      ! Precompute all i-j pairs
      !-----------------------------------------------------------------!
      total_ij_pairs = index_unitcell * number_of_functions
      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      ij_index = 0
      do i = 1, index_unitcell
        do j = 1, number_of_functions
          ij_index = ij_index + 1
          i_index(ij_index) = i
          j_index(ij_index) = j
        end do
      end do
      
      !-----------------------------------------------------------------!
      ! Parallel computation
      !-----------------------------------------------------------------!
      !$omp parallel do private(ij_index,i,j,k,l,AO1,AO2,r1,r2) &
      !$omp shared(NA, AO, i_index, j_index, number_of_atoms, geometry, atoms) &
      !$omp schedule(dynamic,optimal_chunk_size)
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        
        AO1 = AO(i)
        AO2 = AO(j)
        
        r1(1) = AO1%x ; r2(1) = AO2%x
        r1(2) = AO1%y ; r2(2) = AO2%y
        r1(3) = AO1%z ; r2(3) = AO2%z
        
        if (AO1%orbital =="s" .and. AO2%orbital == "s") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_ss_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA(i,j))
            end do 
          end do 
        end if 
        
        if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_sp_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA(i,j))
            end do 
          end do
        end if
        
        if (AO1%orbital(:1) =="p" .and. AO2%orbital == "s") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_sp_toroidal(number_of_atoms,geometry,atoms,r2,r1,AO2,AO1,NA(i,j))
            end do 
          end do
        end if
        
        if (AO1%orbital(:1) =="p" .and. AO2%orbital(:1) == "p") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_pp_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA(i,j))
            end do 
          end do
        end if
        
      end do
      !$omp end parallel do
      
      deallocate(i_index, j_index)

      

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      

      do i = 1 , index_unitcell
        do j = 1 , number_of_functions
          if (abs(NA(i,j)) < 1e-15) NA(i,j) = 0.d0 
        end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
        do j = index_unitcell + 1 , number_of_functions
          NA(i,j) = NA(i-index_unitcell,j-index_unitcell)
        end do 
      end do 

      do i = 1 , number_of_functions - 1 
        do j = i , number_of_functions
          NA(j,i) = NA(i,j)
        end do 
      end do 


end subroutine nuclear_attraction_matrix_toroidal 