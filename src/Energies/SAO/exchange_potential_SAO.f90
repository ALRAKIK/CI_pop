subroutine exchange_potential_SAO(nBas,P,ERI,K)

      ! Compute exchange matrix in the AO basis

      use keywords
      use unitcell_module

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)
      double precision,intent(in)   :: ERI(nfuc,nBas,nBas,nBas)

      ! Local variables

      integer                       :: mu,nu,la,si

      ! Output variables

      double precision,intent(out)  :: K(nBas,nBas)
      
      K(:,:) = 0d0
      
      do nu=1,nBas
        do si=1,nBas
          do la=1,nBas
            do mu=1,nfuc
              K(mu,nu) = K(mu,nu) - 0.5d0 * P(la,si)*ERI(mu,la,si,nu)
            end do
          end do
        end do
      end do

      do mu = 1 , nfuc 
        do nu = 1 , nBas 
          if (dabs(K(mu,nu)) < 1.d-14) K(mu,nu) = 0.d0
        end do 
      end do 

end subroutine exchange_potential_SAO