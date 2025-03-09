subroutine read_basis_class(atom_type)
        
      use atom_basis
      implicit none

      integer             :: i , j 
      type(atom)          :: atom_type
      character(len=100)  :: lines 

      character(len=10)    :: atom_type_charge
      
      write(atom_type_charge,'(A,I0)') "A ", atom_type%charge 

      atom_type%num_s_function = 0
      atom_type%num_p_function = 0

    
      open(1, file="./tmp/Basis_normalized")

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
                read(1,*) atom_type%num_exponent_p , atom_type%num_p_function
                allocate(atom_type%exponent_p(atom_type%num_exponent_p))
                allocate(atom_type%coefficient_p(atom_type%num_exponent_p,atom_type%num_p_function))
               do i = 1 , atom_type%num_exponent_p
                 read(1,*) atom_type%exponent_p(i) , (atom_type%coefficient_p(i,j),j=1,atom_type%num_p_function)
               end do 
              end if 
              if (adjustl(lines(1:1))=="A" .or. adjustl(lines(1:1))=="a") exit
            end do 
          end if 

        end do 

3      close(1)

end subroutine read_basis_class