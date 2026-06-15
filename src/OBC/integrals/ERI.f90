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
      integer                        :: i , j , k , l  , p , q  , actual_total_int
      integer                        :: number_of_functions
      double precision               :: geometry(number_of_atoms,3)
      double precision               :: two_electron(number_of_functions,number_of_functions,number_of_functions,number_of_functions)
      double precision               :: value 

      double precision               :: start_time, end_time
      integer                        :: days, hours, minutes, seconds , t 

      !-----------------------------------------------------------------!
      !                      Schwarz screening                          !
      
      double precision, allocatable :: Q_schwarz(:,:)
      double precision              :: val_iijj
      double precision, parameter   :: schwarz_thresh = 1.d-14
      integer                       :: n_screened, n_computed
      
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

      !-----------------------------------------------------------------!
      !              Schwarz screening precomputation                   !
      !   Q_schwarz(i,j) = sqrt(|(ij|ij)|)  for all unique pairs        !
      !   Done ONCE here, reused for every screening check below        !
      !-----------------------------------------------------------------!

      allocate(Q_schwarz(number_of_functions, number_of_functions))
      Q_schwarz = 0.d0
      write(*,'(A)') 'Precomputing Schwarz screening matrix...'
      write(outfile,'(A)') 'Precomputing Schwarz screening matrix...'


      do i = 1, number_of_functions
        do j = i, number_of_functions
          call ERI_integral_4_function(ERI(i), ERI(j), ERI(i), ERI(j), val_iijj)
          Q_schwarz(i,j) = dsqrt(dabs(val_iijj))
          Q_schwarz(j,i) = Q_schwarz(i,j)
        end do
      end do

      write(*,'(A)') 'Schwarz matrix done.'
      write(outfile,'(A)') 'Schwarz matrix done.'
      flush(outfile)

      write(*,'(A)') 'Sample Q_schwarz values:'
      write(*,'(A,ES12.4)') 'Max Q value = ', maxval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Min Q value = ', minval(Q_schwarz)
      write(*,'(A,ES12.4)') 'Threshold   = ', schwarz_thresh

      write(*,*) 'Will compute ', actual_total_int, ' unique integrals'

!$omp parallel do private(j,k,l,value,p,q) shared(two_electron,ERI)
      
       do i = 1, number_of_functions
         do j = 1, i
             do k = 1, number_of_functions
                 do l = 1, k

                    p = i*(i-1)/2 + j
                    q = k*(k-1)/2 + l

                    if (p >= q) then

                        ! ---- Schwarz screen: ONE multiply + ONE compare ----
                        if (Q_schwarz(i,j) * Q_schwarz(k,l) < schwarz_thresh) then
                          two_electron(i,j,k,l) = 0.d0   ! screened integrals are zero
                          cycle
                        end if

                        call ERI_integral_4_function(ERI(i),ERI(j),ERI(k),ERI(l),value)
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


end subroutine