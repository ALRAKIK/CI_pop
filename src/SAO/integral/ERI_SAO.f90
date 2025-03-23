subroutine ERI_integral_SAO(number_of_atoms,nbas,geometry,atoms,S_SAO_diag)

      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      
      integer,intent(in)             :: number_of_atoms
      integer,intent(in)             :: nbas
      double precision,intent(in)    :: geometry(number_of_atoms,3)
      double precision,intent(in)    :: S_SAO_diag(nbas,nbas)
      integer                        :: number_of_atoms_per_unitcell 
      type(atom)                     :: atoms(number_of_atoms)
      type(ERI_function),allocatable :: ERI  (:)

      integer                        :: i , j , k , l
      
      double precision,allocatable   :: two_electron(:,:,:,:)
      double precision,allocatable   ::      two_eri(:,:,:,:)
      double precision               :: value , Kronecker_delta , phase , norm 
      double precision,parameter     :: pi = 3.14159265358979323846D00
      integer                        :: number_of_functions_per_unitcell 
      integer                        :: kp , kpp , kppp
      integer                        :: mu , rho , sigma  

      !-----------------------------------------------------------------!

      number_of_atoms_per_unitcell = 1 

      number_of_functions_per_unitcell = 0 
      do i = 1 , number_of_atoms_per_unitcell 
        number_of_functions_per_unitcell = number_of_functions_per_unitcell + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(ERI(nbas))
      allocate(two_electron(nbas,nbas,nbas,nbas))
      allocate(     two_eri(nbas,nbas,nbas,nbas))

      call classification(number_of_atoms,nbas,geometry,atoms,ERI)

      two_electron(:,:,:,:) = 0.d0 


      open(1,file="./tmp/ERI.dat")
      open(2,file="./tmp/ERI_SAO.dat")


      ! 2-fold symmetry implementation (k,l permutation only)

      do i = 1, number_of_functions_per_unitcell
        do j = 1, nbas
            do k = 1, nbas
                ! Only calculate for k ≤ l to avoid duplicates
                do l = k, nbas
                    ! Calculate integral once
                    call ERI_integral_4_function_SAO(ERI(i),ERI(j),ERI(k),ERI(l), value)
                    
                    ! Store in original position
                    two_electron(i,j,k,l) = value
                    
                    ! Apply symmetry for k!=l
                    two_electron(i,j,l,k) = value
                end do
            end do
        end do
      end do

      do i = 1, nbas
        do j = 1 , nbas
          do k = 1 , nbas
            do l = 1 , nbas
              if (abs(two_electron(i,j,k,l)) > 1e-10 ) write(2,"(I5,I5,I5,I5,f16.10)") i , j , k , l , two_electron(i,j,k,l)
            end do 
          end do 
        end do 
      end do 

      
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!

      two_eri(:,:,:,:) = 0d0

      do k = 1 , nbas
        do kp = 1 , nbas
          do kpp = 1 , nbas
            do kppp = 1 , nbas

              if (Kronecker_delta(k+kp-kpp-kppp) == 1) then

                norm = (S_SAO_diag(k,k) * S_SAO_diag(kp,kp) * &
                        S_SAO_diag(kpp,kpp) * S_SAO_diag(kppp,kppp))**(-1.d0/2.d0)

              do mu = 1 , nbas 
                do rho =1 , nbas 
                  do sigma = 1 , nbas 
                    phase = ((2*pi)/nbas)*((kp-1)*(mu-1)-(kpp-1)*(rho-1)-(kppp-1)*(sigma-1))
                    two_eri(k,kp,kpp,kppp) =  two_eri(k,kp,kpp,kppp) + &
                                              dcos(phase) * two_electron(1,mu,rho,sigma) * norm * 1/nbas 
                  end do 
                end do 
              end do 

              end if 

            end do 
          end do 
        end do 
      end do 

      do i = 1, nbas
        do j = 1 , nbas
          do k = 1 , nbas
            do l = 1 , nbas
              if (abs(two_eri(i,j,k,l)) > 1e-10 ) write(1,"(I5,I5,I5,I5,f16.10)") i , j , k , l , two_eri(i,j,k,l)
            end do 
          end do 
        end do 
      end do 

      close(1)
      close(2)

      deallocate(ERI)
      deallocate(two_electron)
      deallocate(two_eri)

end subroutine ERI_integral_SAO