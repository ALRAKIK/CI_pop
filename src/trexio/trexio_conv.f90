subroutine trexio_conv_init(calculation_type, n_atoms)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 

      character(len=10),intent(in)    :: calculation_type
      integer          ,intent(in)    :: n_atoms

      ! local ! 

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  
      character(len=100)              :: trexio_file_name


      ! open the trexio file for writing 

      write(trexio_file_name,'(A,A,I0,A)') trim(calculation_type),"_",n_atoms,".trexio"

      call system("rm -r  "// trim(trexio_file_name))

      trexio_file = trexio_open (trim(trexio_file_name), 'w', TREXIO_HDF5, rc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

end subroutine trexio_conv_init

subroutine trexio_conv_close()
   
      use files 
      use trexio
  
      implicit none

      ! local ! 

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  


      ! close the trexio file !

      rc = trexio_close(trexio_file)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

end subroutine trexio_conv_close


subroutine trexio_conv_global(n_atoms,label,geometry,charge,E_nuc,n_electron,&
                              number_of_functions)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 

      integer          ,intent(in)    :: n_atoms
      integer          ,intent(in)    :: charge(n_atoms)
      integer          ,intent(in)    :: n_electron
      integer          ,intent(in)    :: number_of_functions
      character(len=2) ,intent(in)    :: label(n_atoms)
      double precision ,intent(in)    :: geometry(n_atoms,3)
      double precision ,intent(in)    :: E_nuc
      

      ! local ! 

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  


      ! writing into the trexio file ! 

      rc = trexio_write_nucleus_num (trexio_file, n_atoms)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_nucleus_label (trexio_file, label,2)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_nucleus_coord (trexio_file, geometry)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
        print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_nucleus_charge (trexio_file, real(charge,8))
      if (rc /= TREXIO_SUCCESS) then
         call trexio_string_of_error(rc, err_msg)
         print *, 'Error: '//trim(err_msg)
         call exit(-1)
       end if


      rc = trexio_write_nucleus_repulsion (trexio_file, E_nuc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if


      rc = trexio_write_electron_num (trexio_file, n_electron)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_basis_type(trexio_file,"Gaussian",8)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_basis_prim_num(trexio_file, number_of_functions)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_num(trexio_file, number_of_functions)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_cartesian(trexio_file, 1)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

end subroutine trexio_conv_global 

subroutine trexio_conv_integrals(nBas,S,T,V,Hc,ERI)
   
      use files 
      use trexio
  
      implicit none

      ! input ! 
      
      integer          ,intent(in)    :: nBas
      double precision ,intent(in)    :: S(nBas,nBas)
      double precision ,intent(in)    :: T(nBas,nBas)
      double precision ,intent(in)    :: V(nBas,nBas)
      double precision ,intent(in)    :: Hc(nBas,nBas)
      double precision ,intent(in)    :: ERI(nBas,nBas,nBas,nBas)

      ! local ! 

      integer                         :: i , j , k , l 

      double precision                :: Eri_p(nBas,nBas,nBas,nBas)

      integer(trexio_exit_code)       :: rc       
      character*(128)                 :: err_msg  
      integer(8)                      :: offset, icount
      integer(8), parameter           :: BUFSIZE = 10000_8
      integer                         :: buffer_index(4,BUFSIZE)
      double precision                :: buffer_values(BUFSIZE)


      ! writing into the trexio file ! 

      rc = trexio_write_ao_1e_int_overlap(trexio_file, S)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_kinetic(trexio_file, T)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_potential_n_e(trexio_file, V)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if

      rc = trexio_write_ao_1e_int_core_hamiltonian(trexio_file, Hc)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if


      do i = 1 , nBas
        do j = 1 , nBas
          do k = 1 , nBas
            do l = 1 , nBas
              Eri_p(i,j,k,l) = ERI(i,k,j,l)
            end do 
          end do
        end do 
      end do

      icount = 1
      offset = 0_8
      do i = 1 , nBas
        do j = 1 , nBas
          do k = 1 , nBas
            do l = 1 , nBas
              buffer_index(1,icount) = i
              buffer_index(2,icount) = j
              buffer_index(3,icount) = k
              buffer_index(4,icount) = l
              buffer_values(icount) = ERI_p(i,j,k,l)
              icount = icount + 1
            end do 
          end do
        end do 
      end do

      rc = trexio_write_ao_2e_int_eri(trexio_file,offset,BUFSIZE,buffer_index,buffer_values)
      if (rc /= TREXIO_SUCCESS) then
        call trexio_string_of_error(rc, err_msg)
          print *, 'Error: '//trim(err_msg)
        call exit(-1)
      end if



end subroutine trexio_conv_integrals