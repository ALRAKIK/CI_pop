subroutine bary_center_toroidal(e1,e2,r1,r2,rp)

      use torus_init
      use HeavisideModule

      implicit none
      double precision, intent(in)  :: e1, e2
      double precision, intent(in)  :: r1, r2
      double precision, intent(out) :: rp

      ! local !

      double precision,parameter    :: epsilon = 1.d0-10

      ! ----------------------------------------------------------------!

      rp = datan((e1*dsin(ax*r1)+e2*dsin(ax*r2))/(e1*dcos(ax*r1)+e2*dcos(ax*r2)+epsilon))/ax + 0.5d0 * Lx * Heaviside(-e1*dcos(ax*r1)-e2*dcos(ax*r2)) 

      if (dabs(r1-r2) < 0.5*Lx + epsilon .and. dabs(r1-r2) > 0.5*Lx - epsilon ) then

        if (dabs(e1-e2) < 1.d-10) then 
          rp = 0.5d0 * ( r1 + r2 )
        else 
          if (e1 > e2) then 
            rp = r1 
          else if (e2 > e1) then 
            rp = r2 
          else 
            rp = 0.5d0 * ( r1 + r2 )
          end if 
        end if 

      end if

end subroutine bary_center_toroidal

      ! ----------------------------------------------------------------!
      ! ----------------------------------------------------------------!

subroutine symmetry_of_integrals(nf,fpuc,mat_tmp,mat)

      use unitcell_module

      implicit none 

      integer                      :: i , j
      integer                      :: i_cell  , j_cell  , k_cell 
      integer                      :: i_shift , j_shift , k_shift 
      integer                      :: m , n , l
      integer                      :: nf
      integer                      :: fpuc
      integer                      :: func_i , func_j
      integer                      :: equivalent_j
      
      double precision             :: mat_tmp(nf,nf)
      double precision,intent(out) :: mat(nf,nf)

      !-----------------------------------------------------------------!

      mat(:,:) = 0.d0

      ! First, copy what you already computed (first fpuc rows)
      do i = 1, fpuc
        do j = 1, nf
          mat(i,j) = mat_tmp(i,j)
        end do
      end do

      ! Now fill the rest using symmetry
      
      do i = fpuc + 1, nf
        ! Find which unit cell and function this corresponds to
        call find_cell_and_function(i, fpuc, i_cell, j_cell, k_cell, func_i)
        
        do j = 1, nf
          ! Find which unit cell and function j corresponds to
          call find_cell_and_function(j, fpuc, i_shift, j_shift, k_shift, func_j)

          ! The overlap(i,j) should equal overlap(func_i, equivalent_j)
          ! where equivalent_j is the function in the reference unit cell that corresponds to j
          ! relative to i's unit cell
          
          ! Calculate the relative shift
          m = modulo(i_shift - i_cell, nx)
          n = modulo(j_shift - j_cell, ny)  
          l = modulo(k_shift - k_cell, nz)
          
          ! Find the equivalent function in the first fpuc rows
          equivalent_j = m * ny * nz * fpuc + n * nz * fpuc + l * fpuc + func_j
          
          ! Use the symmetry: S(i,j) = S(func_i, equivalent_j)

          if (equivalent_j <= nf) then
            mat(i,j) = mat_tmp(func_i, equivalent_j)
          endif
        end do
      end do

      !-----------------------------------------------------------------!

end subroutine 

subroutine symmetry_of_integrals_ERI(nf, fpuc, eri_tmp, eri)

      use unitcell_module

      implicit none 

      integer                      :: i, j, k, l
      integer                      :: i_cell, j_cell, k_cell
      integer                      :: i_shift, j_shift, k_shift
      integer                      :: m, n, p
      integer                      :: nf
      integer                      :: fpuc
      integer                      :: func_i, func_j, func_k, func_l
      integer                      :: equiv_j, equiv_k, equiv_l

      double precision             :: eri_tmp(fpuc, nf, nf, nf)
      double precision,intent(out) ::     eri(nf, nf, nf, nf)


      !-----------------------------------------------------------------!

      eri(:,:,:,:) = 0.d0
    
      ! First, copy what you already computed (first index up to fpuc)

      do i = 1, fpuc
        do j = 1, nf
          do k = 1, nf
            do l = 1, nf
              eri(i,j,k,l) = eri_tmp(i,j,k,l)
            end do
          end do
        end do
      end do

      do i = 1, nf
        do j = 1, nf
          do k = 1, nf
            do l = 1, nf
              eri(i,j,l,k) = eri(i,j,k,l)
              eri(j,i,k,l) = eri(i,j,k,l)
              eri(j,i,l,k) = eri(i,j,k,l)
              eri(k,l,i,j) = eri(i,j,k,l)
              eri(k,l,j,i) = eri(i,j,k,l)
              eri(l,k,i,j) = eri(i,j,k,l)
              eri(l,k,j,i) = eri(i,j,k,l)
            end do
          end do
        end do
      end do

      ! Now fill the rest using symmetry
    
      do i = fpuc + 1, nf
        ! Find which unit cell and function i corresponds to
        call find_cell_and_function(i, fpuc, i_cell, j_cell, k_cell, func_i)
        
        do j = 1, nf
          ! Find which unit cell and function j corresponds to  
          call find_cell_and_function(j, fpuc, i_shift, j_shift, k_shift, func_j)
          
          ! Calculate relative shift for j
          m = modulo(i_shift - i_cell, nx)
          n = modulo(j_shift - j_cell, ny)  
          p = modulo(k_shift - k_cell, nz)
          
          ! Find equivalent j in reference frame
          equiv_j = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_j
          
          do k = 1, nf
            ! Find which unit cell and function k corresponds to  
            call find_cell_and_function(k, fpuc, i_shift, j_shift, k_shift, func_k)
            
            ! Calculate relative shift for k
            m = modulo(i_shift - i_cell, nx)
            n = modulo(j_shift - j_cell, ny)  
            p = modulo(k_shift - k_cell, nz)
            
            ! Find equivalent k in reference frame
            equiv_k = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_k
            
            do l = 1, nf

              ! Find which unit cell and function l corresponds to  

              call find_cell_and_function(l, fpuc, i_shift, j_shift, k_shift, func_l)
              
              ! Calculate relative shift for l
              m = modulo(i_shift - i_cell, nx)
              n = modulo(j_shift - j_cell, ny)  
              p = modulo(k_shift - k_cell, nz)
              
              ! Find equivalent l in reference frame

              equiv_l = m * ny * nz * fpuc + n * nz * fpuc + p * fpuc + func_l
              
              ! Use the symmetry: (i,j,k,l) = (func_i, equiv_j, equiv_k, equiv_l)

              if (equiv_j <= nf .and. &
                  equiv_k <= nf .and. &
                  equiv_l <= nf) then
                  eri(i,j,k,l) = eri(func_i, equiv_j, equiv_k, equiv_l)
              endif
            end do
          end do
        end do
      end do

      !-----------------------------------------------------------------!

end subroutine symmetry_of_integrals_ERI

subroutine find_cell_and_function(index, fpuc, i_cell, j_cell, k_cell, func_index)

      use unitcell_module
      implicit none

      integer, intent(in) :: index, fpuc
      integer, intent(out) :: i_cell, j_cell, k_cell, func_index
  
      integer :: temp
  
      ! Function index within unit cell (1 to fpuc)
      func_index = mod(index - 1, fpuc) + 1
  
      ! Remove function part and get cell index
      temp = (index - 1) / fpuc
  
      ! Decompose cell index into i,j,k coordinates
      k_cell = mod(temp, nz)
      temp = temp / nz
      j_cell = mod(temp, ny) 
      i_cell = temp / ny
  
end subroutine