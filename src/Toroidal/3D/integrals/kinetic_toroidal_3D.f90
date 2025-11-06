subroutine kinetic_matrix_toroidal_3D(number_of_atoms,number_of_functions,atoms,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: number_of_atoms
      integer                      :: number_of_functions
      integer                      :: fpuc

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      ! output ! 

      double precision             :: kinetic_tmp(number_of_functions,number_of_functions)
      double precision,intent(out) :: kinetic(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!
      
      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function +  3 * atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!


      kinetic(:,:) = 0.d0 

      do i = 1 , fpuc
      !do i = 1 , number_of_functions
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
        
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_ss_toroidal_3D(r1,r2,AO1,AO2,kinetic_tmp(i,j))
              end do 
            end do 

          end if        

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_sp_toroidal_3D(r1,r2,AO1,AO2,kinetic_tmp(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_sp_toroidal_3D(r2, r1, AO2, AO1, kinetic_tmp(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_pp_toroidal_3D(r1, r2, AO1, AO2, kinetic_tmp(i,j))
              end do 
            end do

          end if

        
        end do 
      end do 

!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
!      !                    symmetry of the integrals                    !
!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      call symmetry_of_integrals(number_of_functions,fpuc,kinetic_tmp,kinetic)

      ! do i = 1 , number_of_functions - 1 
      !   do j = 1 , number_of_functions
      !     if (abs(kinetic(i,j)) < 1e-15) kinetic(i,j) = 0.d0 
      !     kinetic(j,i) = kinetic(i,j)
      !   end do 
      ! end do


end subroutine kinetic_matrix_toroidal_3D