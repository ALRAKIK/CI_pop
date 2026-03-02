subroutine read_basis_class(atom_type)
        
      use files
      use atom_basis

      implicit none

      integer             :: i , j 
      type(atom)          :: atom_type
      character(len=100)  :: lines
      character(len=10)   :: atom_type_charge
      logical             :: found_p , found_d  
      
      write(atom_type_charge,'(A,I0)') "A ", atom_type%charge 

      atom_type%num_s_function = 0 ; atom_type%num_exponent_s = 0
      atom_type%num_p_function = 0 ; atom_type%num_exponent_p = 0
      atom_type%num_d_function = 0 ; atom_type%num_exponent_d = 0
      
      found_p                  = .false.
      found_d                  = .false.

    
      open(1,file=trim(tmp_file_name)//"/Basis_normalized")

        do 

          read(1,'(A)',end=3) lines 

          if (lines == atom_type_charge) then 
            do 
              read(1,'(A)',end=3) lines 
              if (lines == "$ S-TYPE FUNCTIONS") then
                read(1,*) atom_type%num_exponent_s , atom_type%num_s_function
                allocate(atom_type%exponent_s(atom_type%num_exponent_s))
                allocate(atom_type%coefficient_s(atom_type%num_exponent_s,atom_type%num_s_function))
               do i = 1 , atom_type%num_exponent_s
                 read(1,*) atom_type%exponent_s(i) , (atom_type%coefficient_s(i,j),j=1,atom_type%num_s_function)
               end do 
              end if 
              if (lines == "$ P-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%exponent_p(atom_type%num_exponent_p))
                allocate(atom_type%coefficient_p(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_p(i,j),j=1,atom_type%num_p_function)
               end do 
              end if
              if (lines == "$ D-TYPE FUNCTIONS") then
                found_d = .true.
                read(1,*) atom_type%num_exponent_d , atom_type%num_d_function
                allocate(atom_type%exponent_d(atom_type%num_exponent_d))
                allocate(atom_type%coefficient_d(atom_type%num_exponent_d,atom_type%num_d_function))
               do i = 1 , atom_type%num_exponent_d
                 read(1,*) atom_type%exponent_d(i) , (atom_type%coefficient_d(i,j),j=1,atom_type%num_d_function)
               end do 
              end if
              if (adjustl(lines(1:1))=="A" .or. adjustl(lines(1:1))=="a") exit
            end do 

            if (.not. found_p) then
              atom_type%num_exponent_p = 0
              atom_type%num_p_function = 0
            end if

            if (.not. found_d) then
              atom_type%num_exponent_d = 0
              atom_type%num_d_function = 0
            end if

            exit

          end if 

        end do 

3      close(1)

end subroutine read_basis_class

subroutine read_basis_class_tor(atom_type)
        
      use files
      use atom_basis
      
      implicit none

      type(atom)          :: atom_type
      integer             :: i , j
      character(len=100)  :: lines 
      character(len=10)   :: atom_type_charge
      logical             :: found_p , found_d
      
      write(atom_type_charge,'(A,I0)') "A ", atom_type%charge 

      atom_type%num_s_function = 0 ; atom_type%num_exponent_s = 0
      atom_type%num_p_function = 0 ; atom_type%num_exponent_p = 0
      atom_type%num_d_function = 0 ; atom_type%num_exponent_d = 0
      
      found_p                  = .false.
      found_d                  = .false.

        open(1,file=trim(tmp_file_name)//"/Basis_normalized")

        do 

          read(1,'(A)',end=3) lines 

          !------------------------------------------------
          ! S-TYPE
          !------------------------------------------------

          if (lines == atom_type_charge) then 
            do 
              read(1,'(A)',end=3) lines 
              if (lines == "$ S-TYPE FUNCTIONS") then
                read(1,*) atom_type%num_exponent_s , atom_type%num_s_function
                allocate(atom_type%exponent_s(atom_type%num_exponent_s))
                allocate(atom_type%coefficient_s(atom_type%num_exponent_s,atom_type%num_s_function))
               do i = 1 , atom_type%num_exponent_s
                 read(1,*) atom_type%exponent_s(i) , (atom_type%coefficient_s(i,j),j=1,atom_type%num_s_function)
               end do 
              end if 
              !------------------------------------------------
              ! PX-TYPE
              !------------------------------------------------
              if (lines == "$ PX-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%exponent_p(atom_type%num_exponent_p))
                allocate(atom_type%coefficient_px(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_px(i,j),j=1,atom_type%num_p_function)
               end do
              end if
              !------------------------------------------------
              ! PY-TYPE
              !------------------------------------------------
              if (lines == "$ PY-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%coefficient_py(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_py(i,j),j=1,atom_type%num_p_function)
               end do
              end if
              !------------------------------------------------
              ! PZ-TYPE
              !------------------------------------------------
              if (lines == "$ PZ-TYPE FUNCTIONS") then
                found_p = .true.
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%coefficient_pz(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_pz(i,j),j=1,atom_type%num_p_function)
               end do
              end if

              if (lines == "$ D-TYPE FUNCTIONS") then
                found_d = .true.
                read(1,*) atom_type%num_exponent_d , atom_type%num_d_function
                allocate(atom_type%exponent_d(atom_type%num_exponent_d))
                allocate(atom_type%coefficient_d(atom_type%num_exponent_d,atom_type%num_d_function))
               do i = 1 , atom_type%num_exponent_d
                 read(1,*) atom_type%exponent_d(i) , (atom_type%coefficient_d(i,j),j=1,atom_type%num_d_function)
               end do
              end if
              if (adjustl(lines(1:1))=="A" .or. adjustl(lines(1:1))=="a") exit
            end do 

            if (.not. found_p) then
              atom_type%num_exponent_p = 0
              atom_type%num_p_function = 0
            end if

            if (.not. found_d) then
              atom_type%num_exponent_d = 0
              atom_type%num_d_function = 0
            end if


            exit 

          end if 

        end do 

3      close(1)

end subroutine read_basis_class_tor