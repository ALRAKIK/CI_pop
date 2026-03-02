subroutine normalize_basis_tor_1D()

      use files 

      implicit none 

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

        !------------------------------------------------------
        ! S-TYPE
        !------------------------------------------------------
      
        if (lines == "$ S-TYPE FUNCTIONS") then
          write(2,'(A)') trim(lines)
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 1 
          call norm_orb_tor_1D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian    ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
      
        !------------------------------------------------------
        ! P-TYPE
        !------------------------------------------------------
      
      else if (lines == "$ P-TYPE FUNCTIONS") then 
          read(1,*) n_gaussian , n_contraction
          allocate(exponent(n_gaussian),contraction(n_gaussian,n_contraction),contractionN(n_gaussian,n_contraction))
          do i = 1 , n_gaussian
            read(1,*) exponent(i) , (contraction(i,j),j=1,n_contraction)
          end do 
          n_type = 2 
          ! --- Px --- ! 
          write(2,'(A)') "$ PX-TYPE FUNCTIONS"
          call norm_orb_tor_px_1D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          ! --- Py --- ! 
          write(2,'(A)') "$ PY-TYPE FUNCTIONS"
          call norm_orb_tor_py_1D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          ! --- Pz --- !
          write(2,'(A)') "$ PZ-TYPE FUNCTIONS"
          call norm_orb_tor_pz_1D(n_gaussian,n_contraction,exponent,contraction,n_type,contractionN,Lx,Ly,Lz)
          write(2,'(4I4)') n_gaussian ,  n_contraction
          do i = 1 , n_gaussian
            write(2,*) exponent(i) , (contractionN(i,j),j=1,n_contraction)
          end do
          deallocate(exponent,contraction,contractionN)
        else
          write(2,'(A)') trim(lines)
        end if
      
      end do 

3     close(1)

      close(2)

end subroutine normalize_basis_tor_1D

subroutine norm_orb_tor_1D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )


      use gsl_bessel_mod
      use bessel_functions

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
      double precision               :: ax   , ax2 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846264338327950288419716939937510d0
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x

      ax  = 2.d0*pi/Lx
      ax2 = ax * ax 

      ! norm s orbital !

      if (n_type == 1) then 

      do n = 1 , n_contraction
        do i = 1 , n_gaussian
         
          Norm(i) = 1.d0 / dsqrt( Lx * (pi/(2.d0*exponent(i))) *  iv_scaled(0,4.d0*exponent(i)/(ax2)) )  

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

               I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)
             
               c1 = contractionN(i,n)
               c2 = contractionN(j,n)

               const = c1 * c2 * Lx

               sum  = sum + const * I_0_gamma_x * (pi/gamma)

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
        do i = 1 , n_gaussian

           Norm(i) = 1.d0 / dsqrt( Lx  * (pi/(2.d0*exponent(i))) *  iv_scaled(1,4.d0*exponent(i)/(ax2)) * (1.d0/(4.d0*exponent(i))) )

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
                
              I_1_gamma_x = iv_scaled(1,2.d0*gamma/ax2)
              
              c1 = contractionN(i,n)
              c2 = contractionN(j,n)
                
              const = c1 * c2 * Lx * (1.d0/(2.d0*gamma))
              sum   = sum + const * I_1_gamma_x * (pi/gamma)
              
            end do 
          end do 
        
        IF (SQRT(sum) .LT. THRMIN) GOTO 30
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contractionN(j,n)*sum
          end do 
        
      end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor_1D


subroutine norm_orb_tor_px_1D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      use gsl_bessel_mod
      use bessel_functions
      
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
      double precision               :: ax   , ax2 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x
      double precision               :: I_1_gamma_x

      ax  = 2.d0*pi/Lx
      ax2 = ax * ax 
      
      ! norm s orbital !

      if (n_type == 1) then 

        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)

              c1 = contraction(i,n)
              c2 = contraction(j,n)

              const = c1 * c2 * Lx

              sum  = sum + const * I_0_gamma_x * (pi/gamma)
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
        do i = 1 , n_gaussian

          Norm(i) = 1.d0 / dsqrt( Lx  * (pi/(2.d0*exponent(i))) *  iv_scaled(0,4.d0*exponent(i)/(ax2)) * (1.d0/(4.d0*exponent(i))) )
          !Norm(i) = 1.d0 / dsqrt( 0.5d0 * Lx / ax2 * (pi/(2.d0*exponent(i))) * (iv_scaled(0,4.d0*exponent(i)/(ax2)) - iv_scaled(2,4.d0*exponent(i)/(ax2)) ) )
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
               if ( gamma == 0.d0 ) cycle

               I_1_gamma_x = iv_scaled(1,2.d0*gamma/ax2)

               c1 = contractionN(i,n)
               c2 = contractionN(j,n)

              const = c1 * c2 * Lx * (1/(2.d0*gamma))
              sum   = sum + const  *  I_1_gamma_x * (pi/gamma)
              
             end do 
           end do 

         IF (SQRT(sum) .LT. THRMIN) GOTO 30
           sum=1.d0/SQRT(sum)
           do j = 1, n_gaussian
             contractionN(j,n)= contractionN(j,n)*sum
           end do 
         end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor_px_1D

subroutine norm_orb_tor_py_1D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      use gsl_bessel_mod
      use bessel_functions

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
      double precision               :: ax   , ax2 
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x

      ax  = 2.d0*pi/Lx
      ax2 = ax * ax 
      
      ! norm s orbital !

      if (n_type == 1) then 

        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)

              c1 = contraction(i,n)
              c2 = contraction(j,n)

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
        do i = 1 , n_gaussian

           Norm(i) = 1.d0 / dsqrt( Lx  * (pi/(2.d0*exponent(i))) *  iv_scaled(0,4.d0*exponent(i)/(ax2)) * (1.d0/(4.d0*exponent(i))) )
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
              if ( gamma == 0.d0 ) cycle

              I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)

              c1 = contractionN(i,n)
              c2 = contractionN(j,n)

              const = c1*c2*Lx*(1/(2.d0*gamma))
              sum  = sum + const*I_0_gamma_x * (pi/gamma)
              
            end do 
          end do 

        IF (SQRT(sum) .LT. THRMIN) GOTO 30
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contractionN(j,n)*sum
          end do 
        end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor_py_1D

subroutine norm_orb_tor_pz_1D(n_gaussian , n_contraction, exponent, contraction  , n_type , contractionN, Lx , Ly , Lz )

      use gsl_bessel_mod
      use bessel_functions

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
      double precision               :: ax   , ax2
      double precision               :: norm(n_gaussian)
      double precision               :: THRMIN = 1.D-17
      double precision , parameter   :: pi = 3.14159265358979323846D00
      double precision               :: I_0_gamma_x

      ax  = 2.d0*pi/Lx
      ax2 = ax * ax 
      
      ! norm s orbital !

      if (n_type == 1) then 

        do n = 1 , n_contraction
          sum = 0.d0
          do i = 1 , n_gaussian
            do j = 1 , n_gaussian

              alpha = exponent(i)
              beta  = exponent(j)
              gamma = alpha + beta

              I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)

              c1 = contraction(i,n)
              c2 = contraction(j,n)

              const = c1 * c2 * Lx

              sum  = sum + const * I_0_gamma_x * (pi/gamma)
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
        do i = 1 , n_gaussian

           Norm(i) = 1.d0 / dsqrt( Lx  * (pi/(2.d0*exponent(i))) *  iv_scaled(0,4.d0*exponent(i)/(ax2)) * (1.d0/(4.d0*exponent(i))) )
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
              if ( gamma == 0.d0 ) cycle

              I_0_gamma_x = iv_scaled(0,2.d0*gamma/ax2)

              c1 = contractionN(i,n)
              c2 = contractionN(j,n)

              const = c1 * c2 * Lx * (1/(2.d0*gamma))
              sum  = sum + const * I_0_gamma_x * (pi/gamma)
              
            end do 
          end do 

        IF (SQRT(sum) .LT. THRMIN) GOTO 30
          sum=1.d0/SQRT(sum)
          do j = 1, n_gaussian
            contractionN(j,n)= contractionN(j,n)*sum
          end do 
        end do

      end if 
      
30    continue 

end subroutine  norm_orb_tor_pz_1D