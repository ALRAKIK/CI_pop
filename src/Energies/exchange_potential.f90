subroutine exchange_potential(nBas,P,ERI,K)

      ! Compute exchange matrix in the AO basis

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

      ! Local variables

      integer                       :: mu,nu,la,si

      ! Output variables

      double precision,intent(out)  :: K(nBas,nBas)
      
      K(:,:) = 0d0


      do nu=1,nBas
        do mu=1,nBas
          do si=1,nBas
            do la=1,nBas
              K(mu,nu) = K(mu,nu) - 0.5d0 * P(la,si)*ERI(mu,la,si,nu) 
            end do
          end do
        end do
      end do

end subroutine exchange_potential

subroutine exchange_potential_spare(nBas,P,K,non_zero,ERI_spare,ERI_index)

      ! Compute exchange matrix in the AO basis

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)

      integer,intent(in)            :: non_zero
      double precision,intent(in)   :: ERI_spare(non_zero)
      integer         ,intent(in)   :: ERI_index(4,non_zero)

      ! Local variables

      integer                       :: mu,nu,la,si
      integer                       :: idx

      ! Output variables

      double precision,intent(out)  :: K(nBas,nBas)
      
      K(:,:) = 0.d0

      do idx = 1, non_zero
        mu = ERI_index(1,idx)
        nu = ERI_index(4,idx)  
        la = ERI_index(2,idx)
        si = ERI_index(3,idx)
        K(mu,nu) = K(mu,nu) - 0.5d0 *  P(la,si) * ERI_spare(idx)
      end do

end subroutine exchange_potential_spare