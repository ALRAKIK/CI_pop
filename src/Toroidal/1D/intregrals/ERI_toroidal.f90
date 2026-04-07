subroutine ERI_integral_toroidal(number_of_atoms,geometry,number_of_functions,atoms,two_electron_integrals)

      use files
      use omp_lib
      use files
      use torus_init
      use atom_basis
      use classification_ERI
      use unitcell_module
      use precomputed_bessel

      implicit none 

      !-----------------------------------------------------------------!

      integer                        :: i , j , k , l
      integer                        :: number_of_atoms
      integer                        :: number_of_functions
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      double precision               :: geometry(number_of_atoms,3)
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision               :: num_total_int
      integer                        :: fpuc

      ! --------------------------------------------------------------- !
      
      double precision,intent(out)   :: two_electron_integrals(number_of_functions,number_of_functions,number_of_functions,number_of_functions)

      ! ----------------------    Time     ---------------------------- !
      integer                        :: days, hours, minutes, seconds , t , time 
      double precision               :: start,end
      double precision               :: start_time, end_time
      !-----------------------------------------------------------------!
      integer                        :: total_ij_pairs, ij_index
      integer, allocatable           :: i_index(:), j_index(:)
      integer                        :: num_threads, optimal_chunk_size
      !-----------------------------------------------------------------!

      integer                        :: current_pct 
      integer                        :: last_pct = -1
      double precision               :: integrals_done = 0.d0
      
      !-----------------------------------------------------------------!

      !                   Precompute Bessel functions                   !

      call initialize_bessel_table_64_Lx()

      !-----------------------------------------------------------------!


      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function + 3 *  atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())


      allocate(ERI(number_of_functions))
      allocate(two_electron(fpuc,number_of_functions,number_of_functions,number_of_functions))

      call classification_tor(number_of_atoms,number_of_functions,geometry,atoms,ERI)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel

      !$omp    parallel
          !$omp     single
            num_threads = omp_get_num_threads()
          !$omp end single
      !$omp end parallel

      ! Calculate optimal chunk size based on available threads

      if (num_threads <= 16) then
        optimal_chunk_size = 16           ! Good for laptops (8-16 cores)
      else if (num_threads <= 64) then
        optimal_chunk_size = 8            ! Good for medium clusters
      else 
        optimal_chunk_size = 1            ! Good for large clusters (128+ cores)
      end if

      start_time = omp_get_wtime()

      ! Precompute all i-j pairs for manual collapse

      total_ij_pairs = 0
      do i = 1, fpuc
      !do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
        end do
      end do

      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))

      total_ij_pairs = 0
      do i = 1, fpuc
      !do i = 1, number_of_functions
        do j = i, number_of_functions
          total_ij_pairs = total_ij_pairs + 1
          i_index(total_ij_pairs) = i
          j_index(total_ij_pairs) = j
        end do
      end do

      num_total_int = 0.d0
      do ij_index = 1, total_ij_pairs
        i = i_index(ij_index)
        j = j_index(ij_index)
        do k = 1, number_of_functions
          do l = k, number_of_functions    
            if (i <= k .or. (i == k .and. j <= l)) then
              num_total_int = num_total_int + 1.d0
            end if
          end do
        end do
      end do

      write(outfile,'(a,I0,a)') 'Need to  compute ', int(num_total_int,16), ' unique integrals (Two electron Integrals)'
      write(*      ,'(a,I0,a)') 'Need to  compute ', int(num_total_int,16), ' unique integrals (Two electron Integrals)'
      write(outfile,*) ''
      flush(outfile)

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

      !$omp parallel do private(ij_index,i,j,k,l) &
      !$omp shared(two_electron, ERI, i_index, j_index) &
      !$omp schedule(dynamic,optimal_chunk_size)

      do ij_index = 1, total_ij_pairs
          i = i_index(ij_index)
          j = j_index(ij_index)
        do k = 1, number_of_functions  
          do l = k, number_of_functions

            if (i <= k .or. (i == k .and. j <= l)) then

              call ERI_integral_4_function_toroidal(ERI(i),ERI(j),ERI(k),ERI(l), two_electron(i,j,k,l))

              !$omp critical
              integrals_done = integrals_done + 1.d0
              current_pct = int((integrals_done * 100.0d0) / num_total_int)
              if (current_pct > last_pct .and. mod(current_pct, 5) == 0) then
                  write(*,'(I3,"% done")', advance='no') current_pct
                  write(*,'(A)') repeat(char(8), 10)
                  flush(6)
                  last_pct = current_pct
              endif
              !$omp end critical

            end if
          end do
        end do

      end do
      
      !$omp end parallel do

      write(*,'(a)') ""
      write(*,'(a)') "*************************************************"
      write(*,'(a)') "* Two-electron integrals calculation completed  *"
      write(*,'(a)') "*************************************************"
      flush(6)
      write(*,'(a)') ""

      deallocate(i_index, j_index)

      end_time = omp_get_wtime()
    
      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,"(a)") "" 
      FLUSH(outfile)

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      call cpu_time(start)
        call symmetry_of_integrals_ERI(number_of_functions,fpuc,two_electron,two_electron_integrals)
      call cpu_time(end)

      time = int(end - start)
      days = (time/86400)
      hours=mod(time,86400)/3600
      minutes=mod(mod(time,86400),3600)/60
      seconds=mod(mod(mod(time,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') 'CPU time for Translational and 8 fold symmetry = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec"
      write(outfile,"(a)") ""
      write(outfile,"(a)") "8-fold symmetry with translation applied to integrals"
      write(outfile,"(a)") ""

      deallocate(ERI)
      deallocate(two_electron)

end subroutine ERI_integral_toroidal