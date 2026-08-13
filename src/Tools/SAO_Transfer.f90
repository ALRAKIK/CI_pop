subroutine transform_AO_to_SAO(nBas, nfuc, Nc, AO_matrix, SAO_matrix)

      use constants_module

      implicit none

      ! input ! 

      integer, intent(in)       :: nBas, nfuc, Nc
      real(dp), intent(in)      :: AO_matrix(nBas, nBas)

      ! local ! 

      integer                   :: k, a, b, m, mu_sao, nu_sao, j
      real(dp)                  :: theta
      complex(dpc)              :: val

      ! output ! 

      complex(dpc), intent(out) :: SAO_matrix(nBas, nBas)
      
      ! / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - / - !


      SAO_matrix = (0.0_dp, 0.0_dp)

      do k = 1, Nc                                                       ! number of cells 
        do a = 1, nfuc                                                   ! number of basis functions in the first cell
          mu_sao = (k-1)*nfuc + a
          do b = 1, nfuc                                                 ! number of basis functions in the first cell
            nu_sao = (k-1)*nfuc + b
            val = (0.0_dp, 0.0_dp)
            do m = 1, Nc                                                 ! number of cells
              j = b + (m-1)*nfuc                                         ! index of the basis function in the m-th cell
              theta = (2.0_dp * pi / dble(Nc)) * dble(k-1) * dble(m-1)
              val = val + exp(dcmplx(0.0_dp, theta)) * AO_matrix(a, j)
            end do
            SAO_matrix(mu_sao, nu_sao) = val
          end do
        end do
      end do

      call clean_and_symmetrize(nBas, SAO_matrix)

end subroutine transform_AO_to_SAO

subroutine transform_SAO_to_AO(nBas, nloc, Nc, M_SAO, M_AO)

      ! Given a matrix M_SAO in the SAO basis, compute its representation in the AO basis.
      ! M_SAO is block-diagonal with blocks of size nloc for each k-point.
      ! M_AO is the full nBas x nBas matrix in the AO basis.

      use constants_module
      implicit none
      integer, intent(in)       :: nBas, nloc, Nc
      complex(dpc), intent(in)  :: M_SAO(nBas, nBas)
      complex(dpc), intent(out) :: M_AO(nBas, nBas)
      integer                   :: k, kp, a, b, m, mp
      integer                   :: mu, nu , mu_sao, nu_sao
      real(dp)                  :: theta, thetap
      complex(dpc)              :: phase, phasep

      M_AO = (0.0_dp, 0.0_dp)

      do k = 1, Nc
        do kp = 1, Nc
          do a = 1, nloc
            do b = 1, nloc
              ! SAO indices
              mu_sao = (k-1)*nloc + a
              nu_sao = (kp-1)*nloc + b
            
              ! For each combination of cells m (for a) and mp (for b)
              do m = 1, Nc
                do mp = 1, Nc
                  mu = a + (m-1)*nloc
                  nu = b + (mp-1)*nloc
                  theta  = (2.0_dp * pi / dble(Nc)) * dble(k-1) * dble(m-1)
                  thetap = (2.0_dp * pi / dble(Nc)) * dble(kp-1) * dble(mp-1)
                  phase  = exp(dcmplx(0.0_dp,  theta))
                  phasep = exp(dcmplx(0.0_dp,  thetap))
                  M_AO(mu, nu) = M_AO(mu, nu) + phase * M_SAO(mu_sao, nu_sao) * conjg(phasep)
                end do
              end do
            end do
          end do
        end do
      end do

      M_AO = M_AO / dble(Nc)   ! because the SAO basis is not normalised

      call clean_and_symmetrize(nBas, M_AO)

end subroutine transform_SAO_to_AO

subroutine clean_and_symmetrize(N, A)

      use constants_module

      implicit none

      integer, intent(in)         :: N
      complex(dpc), intent(inout) :: A(N,N)
      integer                     :: i, j
      real(dp), parameter         :: tol = 1.0d-14

      do i = 1, N
        do j = 1, i-1
          ! Zero out tiny elements (both real and imag)
          if (abs(A(i,j)) < tol) A(i,j) = (0.0_dp, 0.0_dp)
          ! Enforce Hermitian: upper = conjg(lower)
          A(j,i) = conjg(A(i,j))
        end do
        ! Diagonal must be real
        A(i,i) = dble(A(i,i))
      end do
      
end subroutine clean_and_symmetrize


subroutine sort_array(arr, n)

      use constants_module
      implicit none
      integer, intent(in)     :: n
      real(dp), intent(inout) :: arr(n)
      integer                 :: i, j
      real(dp)                :: tmp

      do i = 1, n-1
        do j = i+1, n
          if (arr(i) > arr(j)) then
            tmp = arr(i)
            arr(i) = arr(j)
            arr(j) = tmp
          end if
        end do
      end do

end subroutine sort_array