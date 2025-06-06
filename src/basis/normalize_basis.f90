subroutine normalize_basis()

      implicit none 

      ! input ! 



      ! local ! 

      integer                       :: n_gaussian , n_contraction
      integer                       :: i , j , n_type
      character(len=100)            :: lines

      double precision,allocatable  :: exponent(:) , contraction(:,:) , contractionN(:,:) 

      open(1,file="./tmp/Basis_scratch")
      open(2,file="./tmp/Basis_normalized")

      do
      
        read(1,'(A)',end=3) lines 
      
        write(2,'(A)') trim(lines) 
      
        if (lines == "$ S-TYPE FUNCTIONS") then  
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 1 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
        ! P function ! 
      
        if (lines == "$ P-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 2 
          call norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
      end do 

3     close(1)

      close(2)

end subroutine normalize_basis


subroutine norm_orb(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN)

      implicit none 

      integer,intent(in)            :: n_gaussian , n_contraction

      double precision ,intent(in)  :: exponent(n_gaussian) , contraction(n_gaussian,n_contraction)

      double precision ,intent(out) :: contractionN(n_gaussian,n_contraction)

      ! local ! 

      integer          :: n , i , j , n_type 

      double precision :: D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0   ,&
                          & DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0      ,&
                          & THRMIN = 1.D-17

      double precision :: PIPPI , sum , T

      double precision , parameter :: pi = 3.14159265358979323846D00

      ! code !

      PIPPI = (DP5/PI)**DP75

      do n = 1 , n_contraction
        sum = D0 
        do i = 1 , n_gaussian
          do j = 1 , n_gaussian
            T = D2*SQRT(exponent(i)*exponent(j))/(exponent(i)+exponent(j))
            sum  = sum + contraction(i,n)*contraction(j,n)*(T**(n_type + DP5))
          end do 
        end do 
        IF (SQRT(sum) .LT. THRMIN) GOTO 10
        sum=D1/SQRT(sum)
        do j = 1, n_gaussian
          contractionN(j,n)= contraction(j,n)*sum*(D4*exponent(j))**(DP5*n_type+DP25)*PIPPI
        end do 
        
      end do 

      RETURN
10    continue 
      
      write(*,"(a)") " You Have Zero Norm" 

end subroutine


subroutine normalize_basis_tor()

      implicit none 

      ! input ! 

      ! local ! 

      integer                       :: n_gaussian , n_contraction
      integer                       :: i , j , n_type
      character(len=100)            :: lines

      double precision,allocatable  :: exponent(:) , contraction(:,:) , contractionN(:,:) 
      double precision              :: Lx , Ly , Lz 


      open(4,file="torus_parameters.inp")
        read(4,*) Lx 
      close(4)

      Ly = Lx ; Lz = Lx 

      open(1,file="./tmp/Basis_scratch")
      open(2,file="./tmp/Basis_normalized")
      open(3,file="./tmp/Basis_normalized_p")



      do
      
        read(1,'(A)',end=3) lines 
      
        write(2,'(A)') trim(lines) 
        write(3,'(A)') trim(lines) 
      
        if (lines == "$ S-TYPE FUNCTIONS") then  
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 1 
          call norm_orb_tor(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          write(3,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          do i = 1 , n_gaussian
            write(3,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
        ! P function ! 
      
        if (lines == "$ P-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 2 
          call norm_orb_tor(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          call norm_orb_tor_p(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(3,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(3,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
      end do 

3     close(1)

      close(2)

      close(3)

end subroutine normalize_basis_tor

subroutine normalize_basis_tor_2D()

      implicit none 

      ! input ! 

      ! local ! 

      integer                       :: n_gaussian , n_contraction
      integer                       :: i , j , n_type
      character(len=100)            :: lines

      double precision,allocatable  :: exponent(:) , contraction(:,:) , contractionN(:,:) 
      double precision              :: Lx , Ly , Lz 


      open(4,file="torus_parameters.inp")
        read(4,*) Lx 
      close(4)

      Ly = Lx ; Lz = Lx 

      open(1,file="./tmp/Basis_scratch")
      open(2,file="./tmp/Basis_normalized")

      do
      
        read(1,'(A)',end=3) lines 
      
        write(2,'(A)') trim(lines) 
      
        if (lines == "$ S-TYPE FUNCTIONS") then  
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 1 
          call norm_orb_tor_2D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,'(1000f16.8)') exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
        
      
      end do 

3     close(1)

      close(2)

end subroutine normalize_basis_tor_2D

!subroutine norm_orb_tor(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )
!
!      implicit none 
!
!      integer,intent(in)             :: n_contraction
!      integer,intent(in)             :: n_gaussian
!      integer,intent(in)             :: n_type 
!
!      double precision,intent(in)    :: exponent(n_gaussian)
!      double precision,intent(in)    :: Lx, Ly, Lz
!
!      double precision,intent(in)    :: contraction(n_gaussian,n_contraction)
!      double precision ,intent(out)  :: contractionN(n_gaussian,n_contraction)
!
!      integer                        :: n, i, j 
!      double precision               :: sum  , alpha , beta , gamma , c1 , c2 , const 
!      double precision               :: ax , ay , az 
!      double precision               :: THRMIN = 1.D-17
!      double precision , parameter   :: pi = 3.14159265358979323846D00
!      double precision               :: I_0_gamma_x, I_0_gamma_y, I_0_gamma_z
!
!
!      INTERFACE
!        
!        FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I0
!        END FUNCTION gsl_sf_bessel_I0
!
!        FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I1
!        END FUNCTION gsl_sf_bessel_I1
!
!        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
!        END FUNCTION gsl_sf_bessel_I0_scaled
!
!        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
!        END FUNCTION gsl_sf_bessel_I1_scaled
!
!      END INTERFACE
!
!
!      ax = 2.d0*pi/Lx
!      ay = 2.d0*pi/Ly
!      az = 2.d0*pi/Lz
!      
!      ! norm s orbital !
!
!      if (n_type == 1) then 
!
!        do n = 1 , n_contraction
!          sum = 0.d0
!          do i = 1 , n_gaussian
!            do j = 1 , n_gaussian
!
!              alpha = exponent(i)
!              beta  = exponent(j)
!              gamma = alpha + beta
!
!              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)
!              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
!              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)
!
!              c1 = contraction(i,n)
!              c2 = contraction(j,n)
!
!              const = c1*c2*Lx*Ly*Lz
!
!              sum  = sum + const*I_0_gamma_x*I_0_gamma_y*I_0_gamma_z
!            end do 
!          end do 
!
!          IF (SQRT(sum) .LT. THRMIN) GOTO 20
!          sum=1.d0/SQRT(sum)
!          do j = 1, n_gaussian
!            contractionN(j,n)= contraction(j,n)*sum
!          end do 
!
!        end do 
!
!      end if 
!      
!20    continue 
!
!      ! norm p orbital ! 
!
!      if (n_type == 2) then
!        do n = 1 , n_contraction
!          sum = 0.d0
!          do i = 1 , n_gaussian
!            do j = 1 , n_gaussian
!              alpha = exponent(i)
!              beta  = exponent(j)
!              gamma = alpha + beta
!              if ( gamma == 0.d0 ) cycle
!  
!              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma/ax**2)
!              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
!              I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)
!  
!              c1 = contraction(i,n)
!              c2 = contraction(j,n)
!  
!              const = c1*c2*Lx*Ly*Lz*(1/(2.d0*gamma))
!              sum  = sum + const*I_0_gamma_x*I_0_gamma_y*I_0_gamma_z
!            end do 
!          end do 
!
!
!        IF (SQRT(sum) .LT. THRMIN) GOTO 30
!          sum=1.d0/SQRT(sum)
!
!          do j = 1, n_gaussian
!            contractionN(j,n)= contraction(j,n)*sum
!          end do 
!
!        end do
!
!      end if 
!      
!30    continue 
!
!end subroutine  norm_orb_tor

!subroutine norm_orb_tor(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )
!
!      implicit none 
!
!      integer,intent(in)             :: n_contraction
!      integer,intent(in)             :: n_gaussian
!      integer,intent(in)             :: n_type 
!
!      double precision,intent(in)    :: exponent(n_gaussian)
!      double precision,intent(in)    :: Lx, Ly, Lz
!
!      double precision,intent(in)    :: contraction(n_gaussian,n_contraction)
!      double precision ,intent(out)  :: contractionN(n_gaussian,n_contraction)
!
!      integer                        :: n, i, j 
!      double precision               :: sum  , alpha , beta , gamma , c1 , c2 , const 
!      double precision               :: ax , ay , az 
!      double precision               :: THRMIN = 1.D-17
!      double precision , parameter   :: pi = 3.14159265358979323846D00
!      double precision               :: I_0_gamma_x, I_0_gamma_y, I_0_gamma_z
!
!
!      INTERFACE
!        
!        FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I0
!        END FUNCTION gsl_sf_bessel_I0
!
!        FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I1
!        END FUNCTION gsl_sf_bessel_I1
!
!        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
!        END FUNCTION gsl_sf_bessel_I0_scaled
!
!        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
!          USE iso_c_binding
!          REAL(C_DOUBLE), VALUE :: x_val
!          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
!        END FUNCTION gsl_sf_bessel_I1_scaled
!
!      END INTERFACE
!
!
!      ax = 2.d0*pi/Lx
!      !ay = 2.d0*pi/Ly
!      !az = 2.d0*pi/Lz
!      
!      ! norm s orbital !
!
!      if (n_type == 1) then 
!
!        do n = 1 , n_contraction
!          sum = 0.d0
!          do i = 1 , n_gaussian
!            do j = 1 , n_gaussian
!
!              alpha = exponent(i)
!              beta  = exponent(j)
!              gamma = alpha + beta
!
!              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)
!              !I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
!              !I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)
!
!              c1 = contraction(i,n)
!              c2 = contraction(j,n)
!
!              const = c1*c2*Lx!*Ly*Lz
!
!              sum  = sum + const*I_0_gamma_x*(pi/gamma)!*I_0_gamma_y*I_0_gamma_z
!              print*, const*I_0_gamma_x*(pi/gamma)
!            end do 
!          end do 
!
!          IF (SQRT(sum) .LT. THRMIN) GOTO 20
!          sum=1.d0/SQRT(sum)
!          do j = 1, n_gaussian
!            contractionN(j,n)= contraction(j,n)*sum
!          end do 
!
!        end do 
!
!      end if 
!      
!20    continue 
!
!      ! norm p orbital ! 
!
!      if (n_type == 2) then
!        do n = 1 , n_contraction
!          sum = 0.d0
!          do i = 1 , n_gaussian
!            do j = 1 , n_gaussian
!              alpha = exponent(i)
!              beta  = exponent(j)
!              gamma = alpha + beta
!              if ( gamma == 0.d0 ) cycle
!  
!              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma/ax**2)
!              !I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
!              !I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)
!  
!              c1 = contraction(i,n)
!              c2 = contraction(j,n)
!  
!              !const = c1*c2*Lx*Ly*Lz*(1/(2.d0*gamma))
!              const = c1*c2*Lx*(1/(2.d0*gamma)) !Ly*Lz*(1/(2.d0*gamma))
!              !sum  = sum + const*I_0_gamma_x*I_0_gamma_y*I_0_gamma_z
!              sum  = sum + const*I_0_gamma_x * (pi/gamma)!*I_0_gamma_y*I_0_gamma_z
!            end do 
!          end do 
!
!
!        IF (SQRT(sum) .LT. THRMIN) GOTO 30
!          sum=1.d0/SQRT(sum)
!
!          do j = 1, n_gaussian
!            contractionN(j,n)= contraction(j,n)*sum
!          end do 
!
!        end do
!
!      end if 
!      
!30    continue 
!
!end subroutine  norm_orb_tor

subroutine norm_orb_tor(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      implicit none 

      integer,intent(in)             :: n_contraction
      integer,intent(in)             :: n_gaussian
      integer,intent(in)             :: n_type 

      double precision,intent(in)    :: exponent(n_gaussian)
      double precision,intent(in)    :: Lx, Ly, Lz

      double precision,intent(in)    :: contraction(n_gaussian,n_contraction)
      double precision ,intent(out)  :: contractionN(n_gaussian,n_contraction)

      integer                        :: n, i, j 
      double precision               :: sum  , alpha , beta , gamma , c1 , c2 , const 
      double precision               :: ax 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x


      INTERFACE
        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
        END FUNCTION gsl_sf_bessel_I0_scaled

        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
        END FUNCTION gsl_sf_bessel_I1_scaled

      END INTERFACE


      ax = 2.d0*pi/Lx

      ! norm s orbital !

      if (n_type == 1) then 

      do n = 1 , n_contraction
        do i = 1 , n_gaussian

        Norm(i) = 1.d0 / dsqrt( Lx * (pi/(2.d0*exponent(i))) * gsl_sf_bessel_I0_scaled(4.d0*exponent(i)/(ax*ax)) ) 

        contractionN(i,n) =  contraction(i,n) * Norm(i) 

        end do 
      end do 
        
        do n = 1 , n_contraction
          
          sum = 0.d0

          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)

              c1 = contractionN(i,n)
              c2 = contractionN(j,n)

              const = c1*c2*Lx

              sum  = sum + const*I_0_gamma_x*(pi/gamma)
            end do 
          end do 

          IF (SQRT(sum) .LT. THRMIN) GOTO 20
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contractionN(j,n)*sum
          end do 

        end do 

      end if 
      
20    continue 

      ! norm p orbital ! 

      if (n_type == 2) then

        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian
              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              if ( gamma == 0.d0 ) cycle
  
              I_0_gamma_x = gsl_sf_bessel_I1_scaled(2.d0*gamma/ax**2)
  
              c1 = contraction(i,n)
              c2 = contraction(j,n)
                
              const = c1*c2*Lx*(1/(2.d0*gamma))
              sum  = sum + const*I_0_gamma_x * (pi/gamma)

            end do 
          end do 


        IF (SQRT(sum) .LT. THRMIN) GOTO 30
          sum=1.d0/SQRT(sum)

          do j = 1, n_gaussian
            contractionN(j,n)= contraction(j,n)*sum
          end do 

        end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor

subroutine norm_orb_tor_2D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      implicit none 

      integer,intent(in)             :: n_contraction
      integer,intent(in)             :: n_gaussian
      integer,intent(in)             :: n_type 

      double precision,intent(in)    :: exponent(n_gaussian)
      double precision,intent(in)    :: Lx, Ly, Lz

      double precision,intent(in)    :: contraction(n_gaussian,n_contraction)
      double precision ,intent(out)  :: contractionN(n_gaussian,n_contraction)

      integer                        :: n, i, j 
      double precision               :: sum  , alpha , beta , gamma , c1 , c2 , const 
      double precision               :: ax , ay 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x , I_0_gamma_y


      INTERFACE
        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
        END FUNCTION gsl_sf_bessel_I0_scaled

        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
        END FUNCTION gsl_sf_bessel_I1_scaled

      END INTERFACE


      ax = 2.d0*pi/Lx
      ay = 2.d0*pi/Ly

      ! norm s orbital !

      if (n_type == 1) then 

      do n = 1 , n_contraction
        do i = 1 , n_gaussian

        Norm(i) = 1.d0 / dsqrt( Lx * gsl_sf_bessel_I0_scaled(4.d0*exponent(i)/(ax*ax)) * Ly * gsl_sf_bessel_I0_scaled(4.d0*exponent(i)/(ay*ay)) ) 

        contractionN(i,n) =  contraction(i,n) * Norm(i) 

        end do 
      end do 
        
        do n = 1 , n_contraction
          
          sum = 0.d0

          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)
              I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)


              c1 = contractionN(i,n)
              c2 = contractionN(j,n)

              const = c1*c2*Lx*Ly 

              sum  = sum + const*I_0_gamma_x*I_0_gamma_y
              
            end do 
          end do 

          IF (SQRT(sum) .LT. THRMIN) GOTO 20
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contractionN(j,n)*sum
          end do 

        end do 

      end if 
      
20    continue 

end subroutine  norm_orb_tor_2D


subroutine norm_orb_tor_p(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      implicit none 

      integer,intent(in)             :: n_contraction
      integer,intent(in)             :: n_gaussian
      integer,intent(in)             :: n_type 

      double precision,intent(in)    :: exponent(n_gaussian)
      double precision,intent(in)    :: Lx, Ly, Lz

      double precision,intent(in)    :: contraction(n_gaussian,n_contraction)
      double precision ,intent(out)  :: contractionN(n_gaussian,n_contraction)

      integer                        :: n, i, j 
      double precision               :: sum  , alpha , beta , gamma , c1 , c2 , const 
      double precision               :: ax , ay , az 
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x, I_0_gamma_y, I_0_gamma_z


      INTERFACE
        
        FUNCTION gsl_sf_bessel_I0(x_val) BIND(C, NAME="gsl_sf_bessel_I0")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0
        END FUNCTION gsl_sf_bessel_I0

        FUNCTION gsl_sf_bessel_I1(x_val) BIND(C, NAME="gsl_sf_bessel_I1")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1
        END FUNCTION gsl_sf_bessel_I1

        FUNCTION gsl_sf_bessel_I0_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I0_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I0_scaled
        END FUNCTION gsl_sf_bessel_I0_scaled

        FUNCTION gsl_sf_bessel_I1_scaled(x_val) BIND(C, NAME="gsl_sf_bessel_I1_scaled")
          USE iso_c_binding
          REAL(C_DOUBLE), VALUE :: x_val
          REAL(C_DOUBLE) :: gsl_sf_bessel_I1_scaled
        END FUNCTION gsl_sf_bessel_I1_scaled

      END INTERFACE


      ax = 2.d0*pi/Lx
      !ay = 2.d0*pi/Ly
      !az = 2.d0*pi/Lz
      
      ! norm s orbital !

      if (n_type == 1) then 

        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)
              !I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
              !I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)

              c1 = contraction(i,n)
              c2 = contraction(j,n)

              const = c1*c2*Lx!*Ly*Lz

              sum  = sum + const*I_0_gamma_x*(pi/gamma)!*I_0_gamma_y*I_0_gamma_z
            end do 
          end do 

          IF (SQRT(sum) .LT. THRMIN) GOTO 20
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contraction(j,n)*sum
          end do 

        end do 

      end if 
      
20    continue 

      ! norm p orbital ! 

      if (n_type == 2) then
        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian
              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta
              if ( gamma == 0.d0 ) cycle
  
              I_0_gamma_x = gsl_sf_bessel_I0_scaled(2.d0*gamma/ax**2)
              !I_0_gamma_y = gsl_sf_bessel_I0_scaled(2.d0*gamma/ay**2)
              !I_0_gamma_z = gsl_sf_bessel_I0_scaled(2.d0*gamma/az**2)
  
              c1 = contraction(i,n)
              c2 = contraction(j,n)
  
              !const = c1*c2*Lx*Ly*Lz*(1/(2.d0*gamma))
              const = c1*c2*Lx*(1/(2.d0*gamma)) !Ly*Lz*(1/(2.d0*gamma))
              !sum  = sum + const*I_0_gamma_x*I_0_gamma_y*I_0_gamma_z
              sum  = sum + const*I_0_gamma_x * (pi/gamma)!*I_0_gamma_y*I_0_gamma_z
            end do 
          end do 


        IF (SQRT(sum) .LT. THRMIN) GOTO 30
          sum=1.d0/SQRT(sum)

          do j = 1, n_gaussian
            contractionN(j,n)= contraction(j,n)*sum
          end do 

        end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor_p