subroutine kinetic_matrix_toroidal(number_of_atoms,number_of_functions,atoms,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: index_atom1 , index_sym , index_unitcell 
      integer                      :: number_of_atoms
      integer                      :: number_of_functions

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      ! output ! 

      double precision,intent(out) :: kinetic(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      kinetic(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function

      index_sym   = 0 
    
      do i = 1 , number_of_atoms/2 + 1 
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      index_sym = index_sym + 1 

      index_unitcell = 0 

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 
      
      do i = 1 ,  index_unitcell
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z

          if (AO1%orbital =="s" .and. AO2%orbital == "s") then
        
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_ss_toroidal(r1,r2,AO1,AO2,kinetic(i,j))
              end do 
            end do 

          end if 

          if (AO1%orbital =="s" .and. AO2%orbital(:1) == "p") then
            
            do k = 1 , size  (AO1%exponent)
              do l = 1 , size  (AO2%exponent)
                call kinetic_integral_sp_toroidal(r1,r2,AO1,AO2,kinetic(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital == "s") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_sp_toroidal(r2, r1, AO2, AO1, kinetic(i,j))
              end do 
            end do

          end if

          if (AO1%orbital(:1) == "p" .and. AO2%orbital(:1) == "p") then
                
            do k = 1, size(AO1%exponent)
              do l = 1, size(AO2%exponent)
                call kinetic_integral_pp_toroidal(r1, r2, AO1, AO2, kinetic(i,j))
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
          if (abs(kinetic(i,j)) < 1e-15) kinetic(i,j) = 0.d0 
        end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
        do j = index_unitcell + 1 , number_of_functions
          kinetic(i,j) = kinetic(i-index_unitcell,j-index_unitcell)
        end do 
      end do 

      do i = 1 , number_of_functions - 1 
        do j = i , number_of_functions
          kinetic(j,i) = kinetic(i,j)
        end do 
      end do



end subroutine kinetic_matrix_toroidal


subroutine kinetic_matrix_toroidal_n(number_of_atoms,number_of_functions,atoms,AO,kinetic)

      use files 
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!

      integer                      :: i , j , k , l 
      integer                      :: index_atom1 , index_sym , index_unitcell 
      integer                      :: number_of_atoms
      integer                      :: number_of_functions

      type(atom)                   :: atoms(number_of_atoms)

      type(ERI_function)           :: AO (number_of_functions)
      type(ERI_function)           :: AO1 , AO2


      double precision             :: r1(3) , r2(3)

      integer                      :: pattern_id, encode_orbital_pattern_AO

      ! output ! 

      double precision,intent(out) :: kinetic(number_of_functions,number_of_functions)

      !-----------------------------------------------------------------!

      kinetic(:,:) = 0.d0

      index_atom1 = atoms(1)%num_s_function + 3*atoms(1)%num_p_function

      index_sym   = 0 
    
      do i = 1 , number_of_atoms/2 + 1 
        index_sym = index_sym + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 

      index_sym = index_sym + 1 

      index_unitcell = 0 

      do i = 1 , number_of_atom_in_unitcell
        index_unitcell = index_unitcell  + atoms(i)%num_s_function + 3*atoms(i)%num_p_function
      end do 
      
      do i = 1 , index_unitcell
        do j = 1 , number_of_functions

          AO1 = AO(i)
          AO2 = AO(j)

          r1(1) = AO1%x ; r2(1) = AO2%x
          r1(2) = AO1%y ; r2(2) = AO2%y
          r1(3) = AO1%z ; r2(3) = AO2%z


          pattern_id = encode_orbital_pattern_AO(AO1%orbital, AO2%orbital)

          call Kinetic_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,kinetic(i,j))
          
        end do 
      end do 

      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!
      !                    symmetry of the integrals                    !
      !-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-!


      do i = 1 , index_unitcell
        do j = 1 , number_of_functions
          if (abs(kinetic(i,j)) < 1e-15) kinetic(i,j) = 0.d0 
        end do 
      end do 

      do i = index_unitcell + 1   , number_of_functions
        do j = index_unitcell + 1 , number_of_functions
          kinetic(i,j) = kinetic(i-index_unitcell,j-index_unitcell)
        end do 
      end do 

      do i = 1 , number_of_functions - 1 
        do j = i , number_of_functions
          kinetic(j,i) = kinetic(i,j)
        end do 
      end do



end subroutine kinetic_matrix_toroidal_n


subroutine kinetic_integral_toroidal_1D(pattern_id,r1,r2,AO1,AO2,K_t)

      use torus_init
      use classification_ERI
      use bessel_functions

      implicit none 

      !-----------------------------------------------------------------!

      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: K_t

      integer                          :: i , j 
      double precision,parameter       :: pi = 3.14159265358979323846D00
      double precision                 :: alpha  , beta   , beta2 
      double precision                 :: c1     , c2     , const
      double precision                 :: x1     , x2     , X  
      double precision                 :: y1     , y2     , Y
      double precision                 :: z1     , z2     , Z 
      double precision                 :: kx     , ky     , kz 
      double precision                 :: sx_int , sy_int , sz_int  
      double precision                 :: K 
       

      ! Clifford ! 

      double precision                 :: ax2 , inv_ax , inv_ax2 
      double precision                 :: px
      double precision                 :: I_0_A , I_1_A , I_2_A , I_3_A
      double precision                 :: A , A2 , B , B2 
      double precision                 :: xp
      double precision                 :: spa  , spb
      double precision                 :: spa2 , spb2
      double precision                 :: spa3 , spb3
      double precision                 :: cpa  , cpb
      double precision                 :: cpa2 , cpb2
      double precision                 :: cpa3 , cpb3


      double precision                 :: svp  , cvp
      double precision                 :: svp2 , cvp2
      double precision                 :: svp3 , cvp3
      double precision                 :: svp4 , cvp4


      double precision                 :: scvp11 , scvp21 , scvp12
      double precision                 :: scvp22 , scvp31 , scvp13

      !   Real   ! 

      double precision                 :: albe, inv_albe , mu 
      double precision                 :: yp   , zp
      double precision                 :: ypa  , ypb
      double precision                 :: ypa2 , ypb2
      double precision                 :: ypa3 , ypb3
      double precision                 :: zpa  , zpb
      double precision                 :: zpa2 , zpb2
      double precision                 :: zpa3 , zpb3
      double precision                 :: yvp  , zvp
      double precision                 :: yvp2 , zvp2
      double precision                 :: yvp3 , zvp3
      double precision                 :: yvp4 , zvp4

      
      integer                          :: pattern_id

      double precision                 :: D_00_x , D_00_y , D_00_z
      double precision                 :: D_01_x , D_01_y , D_01_z
      double precision                 :: D_10_x , D_10_y , D_10_z
      double precision                 :: D_11_x , D_11_y , D_11_z
      double precision                 :: S_00_x , S_00_y , S_00_z
      double precision                 :: S_01_x , S_01_y , S_01_z
      double precision                 :: S_10_x , S_10_y , S_10_z
      double precision                 :: S_11_x , S_11_y , S_11_z

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
 
      K   = 0.d0
      K_t = 0.d0 

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

          beta2  = beta * beta

          A      = 2.d0*px/ax2
          B      = 2.d0*beta/ax2
          B2     = B * B 

          I_0_A  = iv_scaled(0, A)
          I_1_A  = iv_scaled(1, A)
          I_2_A  = iv_scaled(2, A)
          I_3_A  = iv_scaled(3, A)

          kx     = dexp(-2.d0 * ( alpha + beta - px ) / ax2)

          sx_int = Lx * I_0_A

          spa    = dsin(ax*(xp-x1)) ; cpa    = dcos(ax*(xp-x1))
          spb    = dsin(ax*(xp-x2)) ; cpb    = dcos(ax*(xp-x2)) 

          spa2 = spa * spa 
          cpa2 = cpa * cpa

          spb2 = spb  * spb
          spb3 = spb2 * spb
          cpb2 = cpb  * cpb 
          cpb3 = cpb2 * cpb 


          ! ************ !
          ! The integral !
          ! ************ ! 

          svp    = 0.d0
          svp2   = Lx * I_1_A / A 
          svp3   = 0.d0
          svp4   = 3 * Lx * I_2_A / (A * A)

          cvp    = Lx * I_1_A 
          cvp2   = Lx * (I_1_A + A * I_2_A) / A 
          cvp3   = Lx * (3.d0 * I_2_A / A + I_3_A)
          cvp4   = Lx * (3.d0  + A * A ) * I_2_A / (A * A)


          scvp11 = 0.d0
          scvp21 = Lx * I_2_A / A 
          scvp12 = 0.d0
          scvp22 = Lx * ( A * I_1_A - 3.d0 * I_2_A ) / (A * A)
          scvp31 = 0.d0 
          scvp13 = 0.d0 
          


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

          sy_int   = dsqrt(pi * inv_albe)
          sz_int   = dsqrt(pi * inv_albe)

          ypa  = yp - y1
          ypa2 = ypa  * ypa 
          
          ypb  = yp - y2
          ypb2 = ypb  * ypb
          ypb3 = ypb2 * ypb

          zpa  = zp - z1 
          zpa2 = zpa * zpa

          zpb  = zp - z2
          zpb2 = zpb * zpb
          zpb3 = zpb2 * zpb

          ! ************ !
          ! the integral !
          ! ************ !

          yvp  = 0.d0 
          yvp2 = 0.5d0  * sy_int * inv_albe
          yvp4 = 0.75d0 * sy_int * inv_albe * inv_albe
          zvp  = 0.d0 
          zvp2 = 0.5d0  * sz_int * inv_albe
          zvp4 = 0.75d0 * sz_int * inv_albe * inv_albe
          

      select case(pattern_id)

      case (00) ! | s     s     ( 1 )

  ! G1 (i=0, j=0, k=0)
  ! G2 (l=0, m=0, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_00_y = sy_int
  S_00_z = sz_int
  S_00_x = sx_int
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_00_x * S_00_y * S_00_z + S_00_x * D_00_y * S_00_z + S_00_x * S_00_y * D_00_z )

case (01) ! | s     px    ( 2 )

  ! G1 (i=0, j=0, k=0)
  ! G2 (l=1, m=0, n=0)

  D_01_x = svp3*B2*ax*cpb3 + 3*scvp21*B2*ax*cpb2*spb + 3*scvp12*B2*ax*cpb*spb2 + cvp3*B2*ax*spb3 - 3*scvp11*B*ax*cpb2 - 3*cvp2*B*ax*cpb*spb + 3*svp2*B*ax*cpb*spb + 3*scvp11*B*ax*spb2 - svp*ax*cpb - cvp*ax*spb
  S_00_y = sy_int
  S_00_z = sz_int
  S_01_x = (svp*cpb + cvp*spb) * inv_ax
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_01_x * S_00_y * S_00_z + S_01_x * D_00_y * S_00_z + S_01_x * S_00_y * D_00_z )

case (02) ! | s     py    ( 3 )

  ! G1 (i=0, j=0, k=0)
  ! G2 (l=0, m=1, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_01_y = sy_int*ypb + yvp
  S_00_z = sz_int
  S_00_x = sx_int
  D_01_y = 4*beta2*sy_int*ypb3 + 12*yvp*beta2*ypb2 + 12*yvp2*beta2*ypb + 4*yvp3*beta2 - 6*beta*sy_int*ypb - 6*yvp*beta
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_00_x * S_01_y * S_00_z + S_00_x * D_01_y * S_00_z + S_00_x * S_01_y * D_00_z )

case (03) ! | s     pz    ( 4 )

  ! G1 (i=0, j=0, k=0)
  ! G2 (l=0, m=0, n=1)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_00_y = sy_int
  S_01_z = sz_int*zpb + zvp
  S_00_x = sx_int
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_01_z = 4*beta2*sz_int*zpb3 + 12*zvp*beta2*zpb2 + 12*zvp2*beta2*zpb + 4*zvp3*beta2 - 6*beta*sz_int*zpb - 6*zvp*beta

  K = -0.5d0 * ( D_00_x * S_00_y * S_01_z + S_00_x * D_00_y * S_01_z + S_00_x * S_00_y * D_01_z )

case (10) ! | px    s     ( 5 )

  ! G1 (i=1, j=0, k=0)
  ! G2 (l=0, m=0, n=0)

  D_10_x = (svp3*B2*ax2*cpa*cpb2 + 2*scvp21*B2*ax2*cpa*cpb*spb + scvp12*B2*ax2*cpa*spb2 + scvp21*B2*ax2*cpb2*spa + 2*scvp12*B2*ax2*cpb*spa*spb + cvp3*B2*ax2*spa*spb2 - scvp11*B*ax2*cpa*cpb + svp2*B*ax2*cpa*spb - cvp2*B*ax2*cpb*spa + scvp11*B*ax2*spa*spb) * inv_ax
  S_00_y = sy_int
  S_00_z = sz_int
  S_10_x = (svp*cpa + cvp*spa) * inv_ax
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_10_x * S_00_y * S_00_z + S_10_x * D_00_y * S_00_z + S_10_x * S_00_y * D_00_z )

case (11) ! | px    px    ( 6 )

  ! G1 (i=1, j=0, k=0)
  ! G2 (l=1, m=0, n=0)

  D_11_x = (svp4*B2*ax*cpa*cpb3 + 3*scvp31*B2*ax*cpa*cpb2*spb + 3*scvp22*B2*ax*cpa*cpb*spb2 + scvp13*B2*ax*cpa*spb3 + scvp31*B2*ax*cpb3*spa + 3*scvp22*B2*ax*cpb2*spa*spb + 3*scvp13*B2*ax*cpb*spa*spb2 + cvp4*B2*ax*spa*spb3 - 3*scvp21*B*ax*cpa*cpb2 - 3*scvp12*B*ax*cpa*cpb*spb + 3*svp3*B*ax*cpa*cpb*spb + 3*scvp21*B*ax*cpa*spb2 - 3*scvp12*B*ax*cpb2*spa - 3*cvp3*B*ax*cpb*spa*spb + 3*scvp21*B*ax*cpb*spa*spb + 3*scvp12*B*ax*spa*spb2 - svp2*ax*cpa*cpb - scvp11*ax*cpa*spb - scvp11*ax*cpb*spa - cvp2*ax*spa*spb) * inv_ax
  S_00_y = sy_int
  S_00_z = sz_int
  S_11_x = (svp2*cpa*cpb + scvp11*cpa*spb + scvp11*cpb*spa + cvp2*spa*spb) * inv_ax2
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_11_x * S_00_y * S_00_z + S_11_x * D_00_y * S_00_z + S_11_x * S_00_y * D_00_z )

case (12) ! | px    py    ( 7 )

  ! G1 (i=1, j=0, k=0)
  ! G2 (l=0, m=1, n=0)

  D_10_x = (svp3*B2*ax2*cpa*cpb2 + 2*scvp21*B2*ax2*cpa*cpb*spb + scvp12*B2*ax2*cpa*spb2 + scvp21*B2*ax2*cpb2*spa + 2*scvp12*B2*ax2*cpb*spa*spb + cvp3*B2*ax2*spa*spb2 - scvp11*B*ax2*cpa*cpb + svp2*B*ax2*cpa*spb - cvp2*B*ax2*cpb*spa + scvp11*B*ax2*spa*spb) * inv_ax
  S_01_y = sy_int*ypb + yvp
  S_00_z = sz_int
  S_10_x = (svp*cpa + cvp*spa) * inv_ax
  D_01_y = 4*beta2*sy_int*ypb3 + 12*yvp*beta2*ypb2 + 12*yvp2*beta2*ypb + 4*yvp3*beta2 - 6*beta*sy_int*ypb - 6*yvp*beta
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_10_x * S_01_y * S_00_z + S_10_x * D_01_y * S_00_z + S_10_x * S_01_y * D_00_z )

case (13) ! | px    pz    ( 8 )

  ! G1 (i=1, j=0, k=0)
  ! G2 (l=0, m=0, n=1)

  D_10_x = (svp3*B2*ax2*cpa*cpb2 + 2*scvp21*B2*ax2*cpa*cpb*spb + scvp12*B2*ax2*cpa*spb2 + scvp21*B2*ax2*cpb2*spa + 2*scvp12*B2*ax2*cpb*spa*spb + cvp3*B2*ax2*spa*spb2 - scvp11*B*ax2*cpa*cpb + svp2*B*ax2*cpa*spb - cvp2*B*ax2*cpb*spa + scvp11*B*ax2*spa*spb) * inv_ax
  S_00_y = sy_int
  S_01_z = sz_int*zpb + zvp
  S_10_x = (svp*cpa + cvp*spa) * inv_ax
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_01_z = 4*beta2*sz_int*zpb3 + 12*zvp*beta2*zpb2 + 12*zvp2*beta2*zpb + 4*zvp3*beta2 - 6*beta*sz_int*zpb - 6*zvp*beta

  K = -0.5d0 * ( D_10_x * S_00_y * S_01_z + S_10_x * D_00_y * S_01_z + S_10_x * S_00_y * D_01_z )

case (20) ! | py    s     ( 9 )

  ! G1 (i=0, j=1, k=0)
  ! G2 (l=0, m=0, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_10_y = sy_int*ypa + yvp
  S_00_z = sz_int
  S_00_x = sx_int
  D_10_y = 4*beta2*sy_int*ypa*ypb2 + 8*yvp*beta2*ypa*ypb + 4*yvp2*beta2*ypa + 4*yvp*beta2*ypb2 + 8*yvp2*beta2*ypb + 4*yvp3*beta2 - 2*beta*sy_int*ypa - 2*yvp*beta
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_00_x * S_10_y * S_00_z + S_00_x * D_10_y * S_00_z + S_00_x * S_10_y * D_00_z )

case (21) ! | py    px    ( 10 )

  ! G1 (i=0, j=1, k=0)
  ! G2 (l=1, m=0, n=0)

  D_01_x = svp3*B2*ax*cpb3 + 3*scvp21*B2*ax*cpb2*spb + 3*scvp12*B2*ax*cpb*spb2 + cvp3*B2*ax*spb3 - 3*scvp11*B*ax*cpb2 - 3*cvp2*B*ax*cpb*spb + 3*svp2*B*ax*cpb*spb + 3*scvp11*B*ax*spb2 - svp*ax*cpb - cvp*ax*spb
  S_10_y = sy_int*ypa + yvp
  S_00_z = sz_int
  S_01_x = (svp*cpb + cvp*spb) * inv_ax
  D_10_y = 4*beta2*sy_int*ypa*ypb2 + 8*yvp*beta2*ypa*ypb + 4*yvp2*beta2*ypa + 4*yvp*beta2*ypb2 + 8*yvp2*beta2*ypb + 4*yvp3*beta2 - 2*beta*sy_int*ypa - 2*yvp*beta
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_01_x * S_10_y * S_00_z + S_01_x * D_10_y * S_00_z + S_01_x * S_10_y * D_00_z )

case (22) ! | py    py    ( 11 )

  ! G1 (i=0, j=1, k=0)
  ! G2 (l=0, m=1, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_11_y = sy_int*ypa*ypb + yvp*ypa + yvp*ypb + yvp2
  S_00_z = sz_int
  S_00_x = sx_int
  D_11_y = 4*beta2*sy_int*ypa*ypb3 + 12*yvp*beta2*ypa*ypb2 + 12*yvp2*beta2*ypa*ypb + 4*yvp3*beta2*ypa + 4*yvp*beta2*ypb3 + 12*yvp2*beta2*ypb2 + 12*yvp3*beta2*ypb + 4*yvp4*beta2 - 6*beta*sy_int*ypa*ypb - 6*yvp*beta*ypa - 6*yvp*beta*ypb - 6*yvp2*beta
  D_00_z = 4*beta2*sz_int*zpb2 + 8*zvp*beta2*zpb + 4*zvp2*beta2 - 2*beta*sz_int

  K = -0.5d0 * ( D_00_x * S_11_y * S_00_z + S_00_x * D_11_y * S_00_z + S_00_x * S_11_y * D_00_z )

case (23) ! | py    pz    ( 12 )

  ! G1 (i=0, j=1, k=0)
  ! G2 (l=0, m=0, n=1)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_10_y = sy_int*ypa + yvp
  S_01_z = sz_int*zpb + zvp
  S_00_x = sx_int
  D_10_y = 4*beta2*sy_int*ypa*ypb2 + 8*yvp*beta2*ypa*ypb + 4*yvp2*beta2*ypa + 4*yvp*beta2*ypb2 + 8*yvp2*beta2*ypb + 4*yvp3*beta2 - 2*beta*sy_int*ypa - 2*yvp*beta
  D_01_z = 4*beta2*sz_int*zpb3 + 12*zvp*beta2*zpb2 + 12*zvp2*beta2*zpb + 4*zvp3*beta2 - 6*beta*sz_int*zpb - 6*zvp*beta

  K = -0.5d0 * ( D_00_x * S_10_y * S_01_z + S_00_x * D_10_y * S_01_z + S_00_x * S_10_y * D_01_z )

case (30) ! | pz    s     ( 13 )

  ! G1 (i=0, j=0, k=1)
  ! G2 (l=0, m=0, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_00_y = sy_int
  S_10_z = sz_int*zpa + zvp
  S_00_x = sx_int
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_10_z = 4*beta2*sz_int*zpa*zpb2 + 8*zvp*beta2*zpa*zpb + 4*zvp2*beta2*zpa + 4*zvp*beta2*zpb2 + 8*zvp2*beta2*zpb + 4*zvp3*beta2 - 2*beta*sz_int*zpa - 2*zvp*beta

  K = -0.5d0 * ( D_00_x * S_00_y * S_10_z + S_00_x * D_00_y * S_10_z + S_00_x * S_00_y * D_10_z )

case (31) ! | pz    px    ( 14 )

  ! G1 (i=0, j=0, k=1)
  ! G2 (l=1, m=0, n=0)

  D_01_x = svp3*B2*ax*cpb3 + 3*scvp21*B2*ax*cpb2*spb + 3*scvp12*B2*ax*cpb*spb2 + cvp3*B2*ax*spb3 - 3*scvp11*B*ax*cpb2 - 3*cvp2*B*ax*cpb*spb + 3*svp2*B*ax*cpb*spb + 3*scvp11*B*ax*spb2 - svp*ax*cpb - cvp*ax*spb
  S_00_y = sy_int
  S_10_z = sz_int*zpa + zvp
  S_01_x = (svp*cpb + cvp*spb) * inv_ax
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_10_z = 4*beta2*sz_int*zpa*zpb2 + 8*zvp*beta2*zpa*zpb + 4*zvp2*beta2*zpa + 4*zvp*beta2*zpb2 + 8*zvp2*beta2*zpb + 4*zvp3*beta2 - 2*beta*sz_int*zpa - 2*zvp*beta

  K = -0.5d0 * ( D_01_x * S_00_y * S_10_z + S_01_x * D_00_y * S_10_z + S_01_x * S_00_y * D_10_z )

case (32) ! | pz    py    ( 15 )

  ! G1 (i=0, j=0, k=1)
  ! G2 (l=0, m=1, n=0)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_01_y = sy_int*ypb + yvp
  S_10_z = sz_int*zpa + zvp
  S_00_x = sx_int
  D_01_y = 4*beta2*sy_int*ypb3 + 12*yvp*beta2*ypb2 + 12*yvp2*beta2*ypb + 4*yvp3*beta2 - 6*beta*sy_int*ypb - 6*yvp*beta
  D_10_z = 4*beta2*sz_int*zpa*zpb2 + 8*zvp*beta2*zpa*zpb + 4*zvp2*beta2*zpa + 4*zvp*beta2*zpb2 + 8*zvp2*beta2*zpb + 4*zvp3*beta2 - 2*beta*sz_int*zpa - 2*zvp*beta

  K = -0.5d0 * ( D_00_x * S_01_y * S_10_z + S_00_x * D_01_y * S_10_z + S_00_x * S_01_y * D_10_z )

case (33) ! | pz    pz    ( 16 )

  ! G1 (i=0, j=0, k=1)
  ! G2 (l=0, m=0, n=1)

  D_00_x = svp2*B2*ax2*cpb2 + 2*scvp11*B2*ax2*cpb*spb + cvp2*B2*ax2*spb2 - cvp*B*ax2*cpb + svp*B*ax2*spb
  S_00_y = sy_int
  S_11_z = sz_int*zpa*zpb + zvp*zpa + zvp*zpb + zvp2
  S_00_x = sx_int
  D_00_y = 4*beta2*sy_int*ypb2 + 8*yvp*beta2*ypb + 4*yvp2*beta2 - 2*beta*sy_int
  D_11_z = 4*beta2*sz_int*zpa*zpb3 + 12*zvp*beta2*zpa*zpb2 + 12*zvp2*beta2*zpa*zpb + 4*zvp3*beta2*zpa + 4*zvp*beta2*zpb3 + 12*zvp2*beta2*zpb2 + 12*zvp3*beta2*zpb + 4*zvp4*beta2 - 6*beta*sz_int*zpa*zpb - 6*zvp*beta*zpa - 6*zvp*beta*zpb - 6*zvp2*beta

  K = -0.5d0 * ( D_00_x * S_00_y * S_11_z + S_00_x * D_00_y * S_11_z + S_00_x * S_00_y * D_11_z )





       


      case default
        K = 0.0d0

      end select

        K_t  =  K_t + const * Kx * Ky * Kz * K 

        end do 
      end do

      !-----------------------------------------------------------------!

end subroutine kinetic_integral_toroidal_1D