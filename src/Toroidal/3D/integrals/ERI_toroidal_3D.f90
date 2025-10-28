subroutine ERI_integral_toroidal_3D(number_of_atoms,geometry,number_of_functions,atoms,two_electron_integrals)

      use files
      use omp_lib
      use files
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l , num_total_int
      integer                        :: number_of_atoms
      integer                        :: number_of_functions
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision,allocatable   :: two_eri(:,:,:,:)
      double precision               :: value
      double precision               :: start_time, end_time
      integer                        :: number_of_functions_per_unitcell 
      integer                        :: index_sym

      integer                        :: days, hours, minutes, seconds , t 

      double precision,intent(out)   :: two_electron_integrals(number_of_functions,number_of_functions,number_of_functions,number_of_functions)

      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size

      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      number_of_functions_per_unitcell = 0 
      do i = 1 , number_of_atom_in_unitcell
        number_of_functions_per_unitcell = number_of_functions_per_unitcell + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(number_of_functions))
      allocate(two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions))
      allocate(two_eri(number_of_functions,number_of_functions,number_of_functions,number_of_functions))

      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

      index_sym   = 0 
    
      do i = 1 , number_of_atoms/2 + 1 
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel


      !$omp parallel
        !$omp single
          num_threads = omp_get_num_threads()
        !$omp end single
      !$omp end parallel


      ! Calculate optimal chunk size based on available threads
      if (num_threads <= 16) then
        optimal_chunk_size = 16    ! Good for laptops (8-16 cores)
      else if (num_threads <= 64) then
        optimal_chunk_size = 8     ! Good for medium clusters
      else 
        optimal_chunk_size = 1     ! Good for large clusters (128+ cores)
      end if

      start_time = omp_get_wtime()

      ! Precompute all i-j pairs for manual collapse


      total_ij_pairs = 0
      !do i = 1, number_of_functions_per_unitcell
      do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
        end do
      end do


      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))


      total_ij_pairs = 0
      !do i = 1, number_of_functions_per_unitcell
      do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
          i_index(total_ij_pairs) = i
          j_index(total_ij_pairs) = j
        end do
      end do


      num_total_int = 0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, number_of_functions
          do l = k, number_of_functions    
            if (i <= k .or. (i == k .and. j <= l)) then
              num_total_int = num_total_int + 1
            end if
          end do
        end do
      end do

      write(*,*) 'Will compute ', num_total_int, ' unique integrals'

      !$omp parallel do private(ij_index,i,j,k,l,value) &
      !$omp shared(two_electron, ERI, i_index, j_index) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do ij_index = 1, total_ij_pairs
          i = i_index(ij_index)
          j = j_index(ij_index)


        do k = 1, number_of_functions
          do l = k, number_of_functions


         if (i <= k .or. (i == k .and. j <= l)) then

           call ERI_integral_4_function_toroidal_3D(ERI(i),ERI(j),ERI(k),ERI(l), value)
               
           two_electron(i,j,k,l) = value
           two_electron(i,j,l,k) = value
           two_electron(j,i,k,l) = value
           two_electron(j,i,l,k) = value
           two_electron(k,l,i,j) = value
           two_electron(k,l,j,i) = value  
           two_electron(l,k,i,j) = value
           two_electron(l,k,j,i) = value

         end if 
       end do
      end do

      end do
      
      !$omp end parallel do


      deallocate(i_index, j_index)

      end_time = omp_get_wtime()
      
      write(outfile,"(a)") ""
      write(outfile,"(a)") "8-fold symmetry with translation applied to integrals"
      write(outfile,"(a)") "" 

      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec     "
      write(outfile,"(a)") "" 
      FLUSH(outfile)

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      !call shift_integrals(two_electron,two_eri,number_of_functions,number_of_functions_per_unitcell)

      two_electron_integrals = 0.d0

      !do i = 1, number_of_functions_per_unitcell
      do i = 1, number_of_functions
        do j = 1 , number_of_functions
          do k = 1 , number_of_functions
            do l = 1 , number_of_functions
              if (abs(two_electron(i,j,k,l)) > 1e-30 )  two_electron_integrals(i,j,k,l) = two_electron(i,j,k,l) 
            end do 
          end do 
        end do 
      end do

      deallocate(ERI)
      deallocate(two_electron)
      deallocate(two_eri)

end subroutine ERI_integral_toroidal_3D