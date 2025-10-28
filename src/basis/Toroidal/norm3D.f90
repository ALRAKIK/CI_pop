subroutine normalize_basis_tor_3D()

      use files

      implicit none 

      ! input ! 

      ! local ! 

      integer                       :: n_gaussian , n_contraction
      integer                       :: i , j , n_type
      character(len=100)            :: lines

      double precision,allocatable  :: exponent(:) , contraction(:,:) , contractionN(:,:) 
      double precision              :: Lx , Ly , Lz 


      open(4,file="torus_parameters.inp")
        read(4,*) Lx , Ly , Lz 
      close(4)

      open(1,file=trim(tmp_file_name)//"/Basis_scratch")
      open(2,file=trim(tmp_file_name)//"/Basis_normalized")

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
          call norm_orb_tor_3D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        end if
      
      end do 

3     close(1)
      close(2)

end subroutine normalize_basis_tor_3D


subroutine norm_orb_tor_3D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      use gsl_bessel_mod

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
      double precision               :: ax  , ay  , az 
      double precision               :: ax2 , ay2 , az2 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x , I_0_gamma_y , I_0_gamma_z  

      ax = 2.d0*pi/Lx
      ay = 2.d0*pi/Ly
      az = 2.d0*pi/Lz 

      ax2 = ax*ax 
      ay2 = ay*ay
      az2 = az*az 

      ! norm s orbital !

      if (n_type == 1) then 

      do n = 1 , n_contraction
        do i = 1 , n_gaussian

        Norm(i) = 1.d0 / dsqrt( Lx * bessi_scaled(0,4.d0*exponent(i)/(ax2)) &
                              * Ly * bessi_scaled(0,4.d0*exponent(i)/(ay2)) & 
                              * Lz * bessi_scaled(0,4.d0*exponent(i)/(az2)) ) 

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

              I_0_gamma_x = bessi_scaled(0,2.d0*gamma/ax2)
              I_0_gamma_y = bessi_scaled(0,2.d0*gamma/ay2)
              I_0_gamma_z = bessi_scaled(0,2.d0*gamma/az2)


              c1 = contractionN(i,n)
              c2 = contractionN(j,n)

              const = c1*c2*Lx*Ly*Lz 

              sum  = sum + const*I_0_gamma_x*I_0_gamma_y*I_0_gamma_z 
              
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

end subroutine  norm_orb_tor_3D