subroutine AO_to_MO_ERI_SAO(nBas,c,ERI_AO,ERI_MO)

      use constants_module
      use unitcell_module
      use files
      implicit none

      integer, intent(in)             :: nBas
      complex(dpc), intent(in)        :: c(nBas,nBas)
      double precision, intent(in)    :: ERI_AO(nBas,nBas,nBas,nBas) ! Still real
      complex(dpc), intent(out)       :: ERI_MO(nBas,nBas,nBas,nBas)

      complex(dpc), allocatable       :: scr(:,:,:,:)
      integer                         :: mu,nu,la,si, p,q,r,s

      allocate(scr(nBas,nBas,nBas,nBas))


      ! 1. Transform 4th index (ket) – no conjugate
      scr = (0d0,0d0)
      do mu=1,nBas
        do nu=1,nBas
          do la=1,nBas
            do si=1,nBas
              do s=1,nBas
                scr(mu,nu,la,s) = scr(mu,nu,la,s) + c(si,s) * ERI_AO(mu,nu,la,si)
              enddo
            enddo
          enddo
        enddo
      enddo
    
      ! 2. Transform 3rd index (ket) – USE CONJG
      ERI_MO = (0d0,0d0)
      do mu=1,nBas
        do nu=1,nBas
          do la=1,nBas
            do r=1,nBas
              do s=1,nBas
                ERI_MO(mu,nu,r,s) = ERI_MO(mu,nu,r,s) + conjg(c(la,r)) * scr(mu,nu,la,s)
              enddo
            enddo
          enddo
        enddo
      enddo
    
      ! 3. Transform 2nd index (bra) – no conjugate
      scr = (0d0,0d0)
      do mu=1,nBas
        do nu=1,nBas
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
                scr(mu,q,r,s) = scr(mu,q,r,s) + c(nu,q) * ERI_MO(mu,nu,r,s)
              enddo
            enddo
          enddo
        enddo
      enddo
    
      ! 4. Transform 1st index (bra) – USE CONJG
      ERI_MO = (0d0,0d0)
      do mu=1,nBas
        do p=1,nBas
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
                ERI_MO(p,q,r,s) = ERI_MO(p,q,r,s) + conjg(c(mu,p)) * scr(mu,q,r,s)
              enddo
            enddo
          enddo
        enddo
      enddo
    
      deallocate(scr)

      call check_eri_hermiticity(nBas,ERI_MO, 1.d-10)
      call check_coulomb_reality(nBas,ERI_MO, 1.d-10)
      call check_eri_symmetry_complex(nBas, ERI_MO, 1.d-12)

end subroutine AO_to_MO_ERI_SAO


subroutine check_eri_hermiticity(nBas, ERI_MO, tol)

      use constants_module
      implicit none
      integer, intent(in)        :: nBas
      complex(dpc), intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
      real(8), intent(in)        :: tol
      integer :: p,q,r,s
      real(8) :: max_abs, max_rel, denom, err_abs, err_rel
      complex(dpc) :: diff

      max_abs = 0d0
      max_rel = 0d0

      do p=1,nBas
        do q=1,nBas
          do r=1,nBas
            do s=1,nBas
              diff = ERI_MO(p,q,r,s) - conjg(ERI_MO(r,s,p,q))
              err_abs = abs(diff)
              denom = max( abs(ERI_MO(p,q,r,s)), abs(ERI_MO(r,s,p,q)), 1d-14 )
              err_rel = err_abs / denom
              if (err_abs > max_abs) max_abs = err_abs
              if (err_rel > max_rel) max_rel = err_rel
            end do
          end do
        end do
      end do

      write(*,'(a,1x,es12.4)') "CHECK ERI Hermiticity max abs err:", max_abs
      write(*,'(a,1x,es12.4)') "CHECK ERI Hermiticity max rel err:", max_rel
      if (max_abs > tol .and. max_rel > tol) then
        write(*,*) "WARNING: ERI Hermiticity check failed (likely wrong conjugation/index order)."
      end if
end subroutine



subroutine check_coulomb_reality(nBas, ERI_MO, tol)
      use constants_module
      implicit none
      integer, intent(in)        :: nBas
      complex(dpc), intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
      real(8), intent(in)        :: tol
      integer :: p,q
      real(8) :: max_im

      max_im = 0d0
      do p=1,nBas
        do q=1,nBas
          max_im = max(max_im, abs(aimag(ERI_MO(p,p,q,q))))
        end do
      end do
    
      write(*,'(a,1x,es12.4)') "CHECK (pp|qq) max |Im|:", max_im
      if (max_im > tol) write(*,*) "WARNING: Large imaginary part in Coulomb terms."
end subroutine


subroutine check_eri_symmetry_complex(nBas, ERI_MO, tol)

  use constants_module
  implicit none
  integer, intent(in)        :: nBas
  complex(dpc), intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  real(8), intent(in)        :: tol
  integer :: p,q,r,s
  real(8) :: max_abs, max_rel, denom, err_abs, err_rel
  complex(dpc) :: diff

  max_abs = 0d0
  max_rel = 0d0

  ! Correct symmetry for complex orbitals: (pq|rs) == (rs|pq)
  do p=1,nBas
    do q=1,nBas
      do r=1,nBas
        do s=1,nBas
          diff = ERI_MO(p,q,r,s) - ERI_MO(r,s,p,q)  ! NO conjugate here!
          err_abs = abs(diff)
          denom = max( abs(ERI_MO(p,q,r,s)), abs(ERI_MO(r,s,p,q)), 1d-14 )
          err_rel = err_abs / denom
          if (err_abs > max_abs) max_abs = err_abs
          if (err_rel > max_rel) max_rel = err_rel
        end do
      end do
    end do
  end do

  write(*,'(a,1x,es12.4)') "CHECK ERI Pair-Swap max abs err:", max_abs
  write(*,'(a,1x,es12.4)') "CHECK ERI Pair-Swap max rel err:", max_rel
  if (max_abs > tol .and. max_rel > tol) then
    write(*,*) "WARNING: ERI Pair-Swap symmetry failed."
  end if

end subroutine check_eri_symmetry_complex