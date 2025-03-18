subroutine nuclear_attraction_matrix_SAO(n_atoms,nbas,geometry,atoms,AO)

      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer           ,intent(in)  :: n_atoms
      integer           ,intent(in)  :: nbas 
      type(ERI_function),intent(in)  :: AO (nbas)
      type(atom)        ,intent(in)  :: atoms(n_atoms)
      double precision  ,intent(in)  :: geometry(n_atoms,3)


      
      type(ERI_function)  :: AO1 , AO2

      integer             :: mu , k , k_prime 
      integer             :: i , j , l  , index_atom1
      
      double precision    :: S(nbas,nbas) , SAO(nbas,nbas)

      double precision    :: r1(3) , r2(3)


      double precision             :: Kronecker_delta
      double precision,parameter   :: pi = 3.14159265358979323846D00

      !-----------------------------------------------------------------!


      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function

      S(:,:) = 0d0 

      do i = 1 , index_atom1
        do j = 1 , nbas 

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call nuclear_attraction_integral_ss_SAO(n_atoms,geometry,atoms,r1,r2,AO1,AO2,S(i,j))
              end do 
            end do 

          end if

        end do 
      end do 

      do k = 1 , nbas 
        do k_prime  = 1 , nbas 
          do mu = 1 , nbas 
            SAO(k,k_prime) = SAO(k,k_prime) + Kronecker_delta(k,k_prime) * dcos(((2*pi)/nbas)*k_prime*(mu-1)) * S(1,mu)
          end do 
        end do 
      end do 

      open(1,file="./tmp/NA.dat")
      do i = 1 , size(SAO,1)
        do j = i , size(SAO,1)
          write(1,'(I5,I5,f16.8)') i , j , SAO(i,j)
        end do 
      end do 
      close(1)

end subroutine nuclear_attraction_matrix_SAO