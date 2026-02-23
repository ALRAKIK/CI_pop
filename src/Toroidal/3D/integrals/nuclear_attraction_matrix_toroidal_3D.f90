subroutine nuclear_attraction_matrix_toroidal_3D(number_of_atoms,number_of_functions,geometry,atoms,AO,NA)

      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use omp_lib

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l
      integer                        :: number_of_atoms
      integer                        :: number_of_functions
      integer                        :: fpuc

      type(atom)                     :: atoms(number_of_atoms)

      type(ERI_function)             :: AO (number_of_functions)
      type(ERI_function)             :: AO1 , AO2

      double precision               :: geometry(number_of_atoms,3)

      double precision               :: r1(3) , r2(3)

      double precision,parameter     :: pi = dacos(-1.d0)

      ! output ! 
      double precision               :: NA_tmp(number_of_functions,number_of_functions)
      double precision,intent(out)   :: NA(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!
      ! Parallelization variables
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      
      NA(:,:) = 0.d0 

      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!

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
      total_ij_pairs = fpuc * number_of_functions
      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      ij_index = 0
      do i = 1, fpuc
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
              call nuclear_attraction_integral_ss_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
            end do 
          end do 
        end if 
        
        if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
            end do 
          end do
        end if
        
        if (AO1%orbital(:1) =="p" .and. AO2%orbital == "s") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r2,r1,AO2,AO1,NA_tmp(i,j))
            end do 
          end do
        end if
        
        if (AO1%orbital(:1) =="p" .and. AO2%orbital(:1) == "p") then
          do k = 1 , size(AO1%exponent)
            do l = 1 , size(AO2%exponent)
              call nuclear_attraction_integral_pp_toroidal(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
            end do 
          end do
        end if
        
      end do
      !$omp end parallel do
      
      deallocate(i_index, j_index)


  
      ! do i = 1 , fpuc
      !   do j = 1 , number_of_functions
        
      !     AO1 = AO(i)
      !     AO2 = AO(j)

      !     r1(1) = AO1%x ; r2(1) = AO2%x
      !     r1(2) = AO1%y ; r2(2) = AO2%y
      !     r1(3) = AO1%z ; r2(3) = AO2%z
      
      !     if (AO1%orbital =="s" .and. AO2%orbital == "s") then
            
      !       do k = 1 , size  (AO1%exponent)
      !         do l = 1 , size  (AO2%exponent)
      !           call nuclear_attraction_integral_ss_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
      !         end do 
      !       end do 

      !     end if

      !     if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
      !         do k = 1 , size  (AO1%exponent)
      !           do l = 1 , size  (AO2%exponent)
      !             call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
      !           end do 
      !         end do

      !     end if

      !     if (AO1%orbital(:1) =="p" .and. AO2%orbital == "s") then
            
      !       do k = 1 , size  (AO1%exponent)
      !         do l = 1 , size  (AO2%exponent)
      !           call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r2,r1,AO2,AO1,NA_tmp(i,j))
      !         end do 
      !       end do
      !     end if

      !     if (AO1%orbital(:1) =="p" .and. AO2%orbital(:1) == "p") then
          
            
      !       do k = 1 , size  (AO1%exponent)
      !         do l = 1 , size  (AO2%exponent)
      !           call nuclear_attraction_integral_pp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
      !         end do 
      !       end do
          
      !     end if
          
      !   end do 
      ! end do 

!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
!      !                    symmetry of the integrals                    !
!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

        call symmetry_of_integrals(number_of_functions,fpuc,NA_tmp,NA)

end subroutine nuclear_attraction_matrix_toroidal_3D