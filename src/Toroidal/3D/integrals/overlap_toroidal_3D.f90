subroutine overlap_matrix_toroidal_3D(number_of_atoms,number_of_functions,atoms,AO,overlap)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use keywords


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

      double precision             :: diag_save(number_of_functions)

      

      ! output !

      double precision             :: overlap_tmp(number_of_functions,number_of_functions)
      double precision,intent(out) :: overlap(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      ! functions_per_unitcell ! 

      fpuc = 0 

      do i = 1 , number_of_atom_in_unitcell 
        fpuc = fpuc + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0

      ! --------------------------------------------------------------- !

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
                call overlap_integral_ss_toroidal_3D(r1,r2,AO1,AO2,overlap_tmp(i,j))
              end do 
            end do 

          end if

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call overlap_integral_sp_toroidal_3D(r1,r2,AO1,AO2,overlap_tmp(i,j))
              end do 
            end do 

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call overlap_integral_sp_toroidal_3D(r2,r1,AO2,AO1,overlap_tmp(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call overlap_integral_pp_toroidal_3D(r1, r2, AO1, AO2, overlap_tmp(i,j))
              end do 
            end do

          end if
          

        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      call symmetry_of_integrals(number_of_functions,fpuc,overlap_tmp,overlap)

      if (c_OV) then 

        open(1,file=trim(tmp_file_name)//"/OV.dat ")
          do i = 1 , size(overlap,1)
            do j = i , size(overlap,1)
              if (abs(overlap(i,j)) > 1e-15 ) write(1,*) i , j , overlap(i,j)
            end do 
          end do 
        close(1)

        open(1,file=trim(tmp_file_name)//"/OV_matrix.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(overlap,1))
        do i = 1 , size(overlap,1)
          write(1,'(i3,6x,1000(f16.12,2x))') i ,  (overlap(i,j),j=1,size(overlap,1))
        end do 
        close(1)

      end if 

end subroutine overlap_matrix_toroidal_3D