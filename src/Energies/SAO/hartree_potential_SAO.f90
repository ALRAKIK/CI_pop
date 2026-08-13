subroutine hartree_potential_SAO(nBas,P,ERI,J)

      ! Compute Coulomb matrix

      use keywords
      use unitcell_module

      implicit none

      ! Input variables

      integer,intent(in)            ::                    nBas
      double precision,intent(in)   ::             P(nBas,nBas)
      double precision,intent(in)   :: ERI(nfuc,nBas,nBas,nBas)


      ! Local variables

      integer                       :: mu,nu,la,si


      ! Output variables

      double precision,intent(out)  :: J(nBas,nBas)

      J(:,:) = 0d0


      do si=1,nBas
        do la=1,nBas
          do nu=1,nBas
            do mu=1,nfuc
              J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,nu,la,si)
            end do
          end do
        end do
      end do

      do mu = 1 , nfuc 
        do nu = 1 , nBas 
          if (dabs(J(mu,nu)) < 1.d-14) J(mu,nu) = 0.d0
        end do 
      end do 


end subroutine hartree_potential_SAO