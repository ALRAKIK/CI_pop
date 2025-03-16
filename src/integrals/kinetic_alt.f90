subroutine kinetic_matrix_alt(number_of_atoms,number_of_functions,geometry,atoms,AO)

      use atom_basis
      use classification_ERI

      implicit none 

      integer                         :: i , j , k , l 
      integer                         :: index_atom1 , index_sym
      integer                         :: group , j_base , j_offset , j_orig
      integer                         :: number_of_atoms
      integer                         :: number_of_functions

      double precision                :: geometry(number_of_atoms,3)

      type(atom)                      :: atoms(number_of_atoms)

      type(ERI_function)              :: AO (number_of_functions)
      type(ERI_function)              :: AO1 , AO2


      double precision,allocatable :: kinetic(:,:)
      double precision             :: r1(3) , r2(3)

      double precision             :: SS 
      double precision             :: SP(3) , PS(3)
      double precision             :: PP(3,3)


      allocate(kinetic(number_of_functions,number_of_functions))

      kinetic(:,:) = 0.d0 

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function

      index_sym   = 0 
    
      do i = 1 , number_of_atoms/2 + 1 
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      index_sym = index_sym + 1 


      do i = 1 , index_atom1
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
        
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_ss_alt(r1,r2,AO1,AO2,kinetic(i,j))
              end do 
            end do 

          end if 

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_sp_alt(r1,r2,AO1,AO2,kinetic(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_sp_alt(r2, r1, AO2, AO1, kinetic(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_pp_alt(r1, r2, AO1, AO2, kinetic(i,j))
              end do 
            end do

          end if
        
        
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      do i = index_atom1 + 1   , number_of_functions
        do j = index_atom1 + 1 , number_of_functions
          kinetic(i,j) = kinetic(i-index_atom1,j-index_atom1)
        end do 
      end do 

      do i = 1 , number_of_functions - 1 
        do j = i , number_of_functions
          kinetic(j,i) = kinetic(i,j)
        end do 
      end do 

      open(1,file="./tmp/KI.dat")
      do i = 1 , size(kinetic,1)
        do j = i , size(kinetic,1)
          write(1,'(I5,I5,f16.8)') i, j , kinetic(i,j)
        end do 
      end do 
      close(1)

      open(1,file="./tmp/KI_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(kinetic,1))
      do i = 1 , size(kinetic,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,   (kinetic(i,j),j=1,size(kinetic,1))
      end do 
      close(1)

      deallocate(kinetic)


end subroutine
