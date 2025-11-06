subroutine nuclear_attraction_matrix_toroidal_3D(number_of_atoms,number_of_functions,geometry,atoms,AO,NA)

      use files
      use torus_init
      use atom_basis
      use classification_ERI

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
      
      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!

      NA(:,:) = 0.d0 
  
      !do i = 1 , number_of_functions

      do i = 1 , fpuc
        do j = 1 , number_of_functions
        
          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z
      
          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call nuclear_attraction_integral_ss_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
              end do 
            end do 

          end if

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
              do k = 1 , size  (AO1%exponent)
                do l = 1 , size  (AO2%exponent)
                  call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
                end do 
              end do

          end if

          if (AO1%orbital(:1) =="p" .and. AO2%orbital == "s") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call nuclear_attraction_integral_sp_toroidal_3D(number_of_atoms,geometry,atoms,r2,r1,AO2,AO1,NA_tmp(i,j))
              end do 
            end do
          end if

          if (AO1%orbital(:1) =="p" .and. AO2%orbital(:1) == "p") then
          
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call nuclear_attraction_integral_pp_toroidal_3D(number_of_atoms,geometry,atoms,r1,r2,AO1,AO2,NA_tmp(i,j))
              end do 
            end do
          
          end if
          
        end do 
      end do 

!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
!      !                    symmetry of the integrals                    !
!      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

        call symmetry_of_integrals(number_of_functions,fpuc,NA_tmp,NA)

      ! do i = 1 , number_of_functions - 1
      !   do j = 1 , number_of_functions
      !     if (abs(NA(i,j)) < 1e-15) NA(i,j) = 0.d0 
      !     NA(j,i) = NA(i,j)
      !   end do 
      ! end do


end subroutine nuclear_attraction_matrix_toroidal_3D