subroutine overlap_matrix_toroidal(number_of_atoms,number_of_functions,atoms,AO,overlap)

      use files
      use atom_basis
      use torus_init
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: index_atom1 , index_sym
      integer                      :: index_unitcell
      integer                      :: number_of_atoms
      integer                      :: number_of_functions

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2

      
      double precision             :: r1(3) , r2(3)


      double precision,parameter   :: pi = 3.14159265358979323846D00


      ! output !

      double precision,intent(out) :: overlap(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function
      
      index_unitcell = 0

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      index_sym   = 0

      do i = 1 , number_of_atoms/2 + 1
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      ! --------------------------------------------------------------- !

      do i = 1 , index_unitcell
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z
          
          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call overlap_integral_ss_toroidal(r1,r2,AO1,AO2,overlap(i,j))
              end do 
            end do 

          end if

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call overlap_integral_sp_toroidal(r1,r2,AO1,AO2,overlap(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call overlap_integral_sp_toroidal(r2, r1, AO2, AO1, overlap(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call overlap_integral_pp_toroidal(r1, r2, AO1, AO2, overlap(i,j))
              end do 
            end do

          end if

        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      do i = 1 , index_unitcell
       do j = 1 , number_of_functions
         if (abs(overlap(i,j)) < 1e-15) overlap(i,j) = 0.d0 
       end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
       do j = index_unitcell + 1 , number_of_functions
         overlap(i,j) = overlap(i-index_unitcell,j-index_unitcell)
       end do 
      end do 
 
      do i = 1 , number_of_functions - 1 
       do j = i , number_of_functions
         overlap(j,i) = overlap(i,j)
       end do 
      end do

end subroutine overlap_matrix_toroidal

      !-----------------------------------------------------------------!


subroutine overlap_matrix_toroidal_n(number_of_atoms,number_of_functions,atoms,AO,overlap)

      use files
      use atom_basis
      use torus_init
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j
      integer                      :: index_atom1 , index_sym
      integer                      :: index_unitcell
      integer                      :: number_of_atoms
      integer                      :: number_of_functions

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2

      
      double precision             :: r1(3) , r2(3)

      integer                      :: pattern_id, encode_orbital_pattern_AO

      ! output !

      double precision,intent(out) :: overlap(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      overlap(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function
      
      index_unitcell = 0

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      index_sym   = 0

      do i = 1 , number_of_atoms/2 + 1
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do

      ! --------------------------------------------------------------- !

      do i = 1 , index_unitcell
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call overlap_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,overlap(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      do i = 1 , index_unitcell
       do j = 1 , number_of_functions
         if (abs(overlap(i,j)) < 1e-15) overlap(i,j) = 0.d0 
       end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
       do j = index_unitcell + 1 , number_of_functions
         overlap(i,j) = overlap(i-index_unitcell,j-index_unitcell)
       end do 
      end do 
 
      do i = 1 , number_of_functions - 1 
       do j = i , number_of_functions
         overlap(j,i) = overlap(i,j)
       end do 
      end do

end subroutine overlap_matrix_toroidal_n

      !-----------------------------------------------------------------!

subroutine overlap_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,S_prime)

      use torus_init
      use classification_ERI
      use bessel_functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: S_prime

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha  , beta
      double precision                 :: c1     , c2     , const
      double precision                 :: x1     , x2     , X  
      double precision                 :: y1     , y2     , Y
      double precision                 :: z1     , z2     , Z 
      double precision                 :: kx     , ky     , kz
      double precision                 :: sx     , sy     , sz  
      double precision                 :: sx_int , sy_int , sz_int  
      double precision                 :: s 
       

      ! Clifford ! 

      double precision                 :: ax2 , inv_ax , inv_ax2 
      double precision                 :: px
      double precision                 :: I_0_A , I_1_A , I_2_A
      double precision                 :: A 
      double precision                 :: xp
      double precision                 :: spa  , spb
      double precision                 :: cpa  , cpb
      double precision                 :: svp  , cvp
      double precision                 :: svp2 , cvp2

      !   Real   ! 

      double precision                 :: albe, inv_albe , mu 
      double precision                 :: yp   , zp
      double precision                 :: ypa  , ypb 
      double precision                 :: zpa  , zpb 
      double precision                 :: yvp  , zvp
      double precision                 :: yvp2 , zvp2

      
      integer                          :: pattern_id

      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!

      x1  = r1(1) ; x2  = r2(1) 
      y1  = r1(2) ; y2  = r2(2)
      z1  = r1(3) ; z2  = r2(3)

      X   = (x1 - x2)
      Y   = (y1 - y2)
      Z   = (z1 - z2)

      ax2 = ax * ax 

      inv_ax  = 1.d0 / ax 
      inv_ax2 = inv_ax * inv_ax 

      !-----------------------------------------------------------------!
 
      S       = 0.d0
      s_prime = 0.d0 

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        c1    = AO1%coefficient(i)
        do j = 1 , size(AO2%exponent)
          beta = AO2%exponent(j)
          c2   = AO2%coefficient(j)

          const   =  c1 * c2

          ! ----------------- !
          ! Clifford Gaussian !
          ! ----------------- ! 

          call bary_exponent_x         (alpha,beta,X,px)
          call bary_center_toroidal_x  (alpha,beta,x1,x2,xp)

          A      = 2.d0*px/ax2

          I_0_A  = iv_scaled(0, A)
          I_1_A  = iv_scaled(1, A)
          I_2_A  = iv_scaled(2, A)

          kx     = dexp(-2.d0 * ( alpha + beta - px ) / ax2)

          sx_int = Lx * I_0_A

          spa    = dsin(ax*(xp-x1)) ; cpa    = dcos(ax*(xp-x1))
          spb    = dsin(ax*(xp-x2)) ; cpb    = dcos(ax*(xp-x2)) 

          ! ************ !
          ! The integral !
          ! ************ ! 

          svp    = 0.d0
          svp2   = Lx * I_1_A / A 

          cvp    = Lx * I_1_A 
          cvp2   = Lx * (I_1_A + A * I_2_A) / A 


          ! ----------------- !
          !   Real Gaussian   !
          ! ----------------- !

          albe     = alpha + beta
          inv_albe = 1.d0 / albe
          mu       = alpha * beta * inv_albe
          yp       = (alpha * y1 + beta * y2) * inv_albe
          zp       = (alpha * z1 + beta * z2) * inv_albe

          ky       = dexp(- mu * ( Y * Y ))
          kz       = dexp(- mu * ( Z * Z ))

          sy_int       = dsqrt(pi * inv_albe)
          sz_int       = dsqrt(pi * inv_albe)

          ypa = yp - y1
          ypb = yp - y2

          zpa = zp - z1 
          zpb = zp - z2

          ! ************ !
          ! the integral !
          ! ************ !

          yvp  = 0.d0 
          yvp2 = 0.5d0 * inv_albe * sy_int
          zvp  = 0.d0 
          zvp2 = 0.5d0 * inv_albe * sz_int
          

      select case(pattern_id)

      ! case (00) ! | s   s    ( 1 )

      ! zeta_x    = sx
      ! zeta_y    = sy
      ! zeta_z    = sz

      ! case (01) ! | s   px   ( 2 )
      
      ! zeta_x    = (cpb*svp + cvp*spb) * inv_ax
      ! zeta_y    = (sy)
      ! zeta_z    = (sz)
      
      ! case (02) ! | s   py   ( 3 )
      
      ! zeta_x    = (sx)
      ! zeta_y    = (ypb + yvp) * sy
      ! zeta_z    = (sz)
      
      ! case (03) ! | s   pz   ( 4 )
      
      ! zeta_x    = (sx)
      ! zeta_y    = (sy)
      ! zeta_z    = (zpb + zvp) * sz
      
      ! case (10) ! | px  s    ( 5 )
      
      ! zeta_x    = (cpa*svp + cvp*spa) * inv_ax
      ! zeta_y    = (sy)
      ! zeta_z    = (sz)
      
      ! case (11) ! | px  px   ( 6 )
      
      ! zeta_x    = (cpa*cpb*svp2 + cpa*cvp*spb*svp + cpb*cvp*spa*svp + cvp2*spa*spb) * inv_ax2
      ! zeta_y    = sy
      ! zeta_z    = sz
      
      ! case (12) ! | px  py   ( 7 )
      
      ! zeta_x    = (cpa*svp + cvp*spa) * inv_ax
      ! zeta_y    = (ypb + yvp) * sy
      ! zeta_z    = sz
      
      ! case (13) ! | px  pz   ( 8 )
      
      ! zeta_x    = (cpa*svp + cvp*spa) * inv_ax
      ! zeta_y    = sy
      ! zeta_z    = (zpb + zvp) * sz
      
      ! case (20) ! | py  s    ( 9 )
      
      ! zeta_x    = (sx)
      ! zeta_y    = (ypa + yvp) * sy
      ! zeta_z    = (sz)
      
      ! case (21) ! | py  px   ( 10)
      
      ! zeta_x    = (cpb*svp + cvp*spb) * inv_ax
      ! zeta_y    = (ypa + yvp) * sy
      ! zeta_z    = sz
      
      ! case (22) ! | py  py   ( 11)
      
      ! zeta_x    = sx
      ! zeta_y    = (ypa*ypb + ypa*yvp + ypb*yvp + yvp2) * sy
      ! zeta_z    = sz
      
      ! case (23) ! | py  pz   ( 12)
      
      ! zeta_x    = sx
      ! zeta_y    = (ypa + yvp) * sy
      ! zeta_z    = (zpb + zvp) * sz
      
      ! case (30) ! | pz  s    ( 13)
      
      ! zeta_x    = (sx)
      ! zeta_y    = (sy)
      ! zeta_z    = (zpa + zvp) * sz
      
      ! case (31) ! | pz  px   ( 14)
      
      ! zeta_x    = (cpb*svp + cvp*spb) * inv_ax
      ! zeta_y    = sy
      ! zeta_z    = (zpa + zvp) * sz
      
      ! case (32) ! | pz  py   ( 15)
      
      ! zeta_x    = sx
      ! zeta_y    = (ypb + yvp) * sy
      ! zeta_z    = (zpa + zvp) * sz
      
      ! case (33) ! | pz  pz   ( 16)
      
      ! zeta_x    = sx
      ! zeta_y    = sy
      ! zeta_z    = (zpa*zpb + zpa*zvp + zpb*zvp + zvp2) * sz

      case (00) ! | s     s     ( 1 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx = sx_int
        sy = sy_int
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (01) ! | s     px    ( 2 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx = (svp*cpb + cvp*spb) * inv_ax
        sy = sy_int
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (02) ! | s     py    ( 3 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx = sx_int
        sy = sy_int*ypb + yvp
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (03) ! | s     pz    ( 4 )
      
        ! G1 (i=0, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx = sx_int
        sy = sy_int
        sz = sz_int*zpb + zvp
      
        s  = sx * sy * sz
      
      case (10) ! | px    s     ( 5 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx = (svp*cpa + cvp*spa) * inv_ax
        sy = sy_int
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (11) ! | px    px    ( 6 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx = (svp2*cpa*cpb + svp*cvp*cpa*spb + svp*cvp*cpb*spa + cvp2*spa*spb) * inv_ax2
        sy = sy_int
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (12) ! | px    py    ( 7 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx = (svp*cpa + cvp*spa) * inv_ax
        sy = sy_int*ypb + yvp
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (13) ! | px    pz    ( 8 )
      
        ! G1 (i=1, j=0, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx = (svp*cpa + cvp*spa) * inv_ax
        sy = sy_int
        sz = sz_int*zpb + zvp
      
        s  = sx * sy * sz
      
      case (20) ! | py    s     ( 9 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=0)
      
        sx = sx_int
        sy = sy_int*ypa + yvp
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (21) ! | py    px    ( 10 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=1, m=0, n=0)
      
        sx = (svp*cpb + cvp*spb) * inv_ax
        sy = sy_int*ypa + yvp
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (22) ! | py    py    ( 11 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=1, n=0)
      
        sx = sx_int
        sy = sy_int*ypa*ypb + yvp*ypa + yvp*ypb + yvp2
        sz = sz_int
      
        s  = sx * sy * sz
      
      case (23) ! | py    pz    ( 12 )
      
        ! G1 (i=0, j=1, k=0)
        ! G2 (l=0, m=0, n=1)
      
        sx = sx_int
        sy = sy_int*ypa + yvp
        sz = sz_int*zpb + zvp
      
        s  = sx * sy * sz
      
      case (30) ! | pz    s     ( 13 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=0)
      
        sx = sx_int
        sy = sy_int
        sz = sz_int*zpa + zvp
      
        s  = sx * sy * sz
      
      case (31) ! | pz    px    ( 14 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=1, m=0, n=0)
      
        sx = (svp*cpb + cvp*spb) * inv_ax
        sy = sy_int
        sz = sz_int*zpa + zvp
      
        s  = sx * sy * sz
      
      case (32) ! | pz    py    ( 15 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=1, n=0)
      
        sx = sx_int
        sy = sy_int*ypb + yvp
        sz = sz_int*zpa + zvp
      
        s  = sx * sy * sz
      
      case (33) ! | pz    pz    ( 16 )
      
        ! G1 (i=0, j=0, k=1)
        ! G2 (l=0, m=0, n=1)
      
        sx = sx_int
        sy = sy_int
        sz = sz_int*zpa*zpb + zvp*zpa + zvp*zpb + zvp2
      
        s  = sx * sy * sz

      case default
        s = 0.0d0

      end select

        s_prime  =  s_prime + const * kx * ky * kz * s 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine overlap_integral_toroidal_1D