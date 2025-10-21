subroutine ERI_integral(number_of_atoms,number_of_functions,geometry,atoms,two_electron)

      use files 
      use omp_lib
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      type(atom),intent(in)          :: atoms(number_of_atoms)
      integer,intent(in)             :: number_of_atoms


      type(ERI_function),allocatable :: ERI  (:)
      integer                        :: i , j , k , l  , p , q  , actual_total_int , num_int
      integer                        :: number_of_functions
      double precision               :: geometry(number_of_atoms,3)
      double precision               :: two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions)
      double precision               :: value 

      double precision               :: start_time, end_time
      integer                        :: days, hours, minutes, seconds , t 
      
      !-----------------------------------------------------------------!

      call omp_set_dynamic(.false.)
      call omp_set_num_threads(omp_get_max_threads())

      
      allocate(ERI(number_of_functions))
      
      call classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)

      !$omp parallel
      if (omp_get_thread_num() == 0) then
        print *, "Running with", omp_get_num_threads(), "threads"
      endif
      !$omp end parallel

      start_time = omp_get_wtime()

      ! 8-fold symmetry implementation ((i,j),(k,l) permutation only)

      actual_total_int = 0 

      do i = 1, number_of_functions
        do j = 1, i
            do k = 1, number_of_functions
                do l = 1, k
                    p = i*(i-1)/2 + j
                    q = k*(k-1)/2 + l
                    if (p >= q) actual_total_int = actual_total_int + 1
                end do
            end do
        end do
      end do

      write(*,*) 'Will compute ', actual_total_int, ' unique integrals'

!$omp parallel do private(j,k,l,value,p,q) shared(two_electron,ERI)
      
      do i = 1, number_of_functions
        do j = 1, i
            do k = 1, number_of_functions
                do l = 1, k
                    p = i*(i-1)/2 + j
                    q = k*(k-1)/2 + l

                    if (p >= q) then
                      call ERI_integral_4_function(ERI(i),ERI(j),ERI(k),ERI(l),value)

                        ! Store only in the canonical position
                        !$omp atomic
                        num_int = num_int + 1
                        two_electron(i,j,k,l) = value
                        two_electron(j,i,k,l) = value
                        two_electron(i,j,l,k) = value     
                        two_electron(j,i,l,k) = value
                        two_electron(k,l,i,j) = value     
                        two_electron(l,k,i,j) = value
                        two_electron(k,l,j,i) = value
                        two_electron(l,k,j,i) = value


                    end if
                end do
            end do
        end do
      end do
!$omp end parallel do

      end_time = omp_get_wtime()

      write(outfile,"(a)") ""
      write(outfile,"(a)") " (8-fold) symmetry applied to integrals"
      write(outfile,"(a)") "" 

      t = int(end_time - start_time)
      days= (t/86400)
      hours=mod(t,86400)/3600
      minutes=mod(mod(t,86400),3600)/60
      seconds=mod(mod(mod(t,86400),3600),60)

      write(outfile,'(A65,5X,I0,a,I0,a,I0,a,I0,4x,a)') '2 el-integrals calculation time = ',days,":",hours,":",minutes,":",seconds, "days:hour:min:sec     "
      write(outfile,"(a)") "" 

      open(1,file=trim(tmp_file_name)//"/ERI.dat")
        do i = 1, number_of_functions
          do j = 1 , number_of_functions
            do k = 1 , number_of_functions
              do l = 1 , number_of_functions
                if (abs(two_electron(i,j,k,l)) > 1e-24 ) write(1,*) i , j , k , l , two_electron(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine