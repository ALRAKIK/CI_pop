subroutine hartree_potential(nBas,P,ERI,J)

      ! Compute Coulomb matrix

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
      double precision,intent(in)   :: P(nBas,nBas)
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)


      ! Local variables

      integer                       :: mu,nu,la,si


      ! Output variables

      double precision,intent(out)  :: J(nBas,nBas)

      J(:,:) = 0d0

      do nu=1,nBas
        do mu=1,nBas
          do si=1,nBas
            do la=1,nBas
              J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,nu,la,si)
            end do
          end do
        end do
      end do


end subroutine hartree_potential


subroutine hartree_potential_spare(nBas,P,J,non_zero,ERI_spare,ERI_index)

      ! Compute Coulomb matrix

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

      double precision,intent(out)  :: J(nBas,nBas)

      J(:,:) = 0.d0

      do idx = 1, non_zero
        mu = ERI_index(1,idx)
        nu = ERI_index(2,idx)  
        la = ERI_index(3,idx)
        si = ERI_index(4,idx)
        J(mu,nu) = J(mu,nu) + P(la,si) * ERI_spare(idx)
      end do



end subroutine hartree_potential_spare