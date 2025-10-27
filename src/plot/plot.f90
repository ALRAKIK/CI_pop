subroutine plot(n_atoms, geometry, atoms)

      use atom_basis
      use classification_ERI
      use torus_init
      use files 
      implicit none 

      integer, intent(in)          :: n_atoms 
      double precision, intent(in) :: geometry(n_atoms, 3)
      type(atom), intent(in)       :: atoms(n_atoms)

      ! local variables
      integer                       :: i , j , k , gr_column , gt_column
      double precision              :: x, x_min, x_max, step
      double precision,allocatable  ::  gr(:),  gt(:)
      double precision,allocatable  :: grp(:), gtp(:)
      double precision              :: center_x
      integer                       :: center_index
      double precision, parameter   :: pi = 3.14159265358979323846D00
      integer                       :: n_points

      double precision            :: exponent , coefficient

      ! --------------------------------------------------------------- !

      center_index = n_atoms/2


      print*, atoms(center_index)%num_exponent_s
      print*, atoms(center_index)%exponent_s
      print*, atoms(center_index)%num_s_function
      print*, atoms(center_index)%coefficient_s(:,1)

      allocate( gr(atoms(center_index)%num_s_function))
      allocate( gt(atoms(center_index)%num_s_function))
      allocate(grp(atoms(center_index)%num_p_function))
      allocate(gtp(atoms(center_index)%num_p_function))

      step     = 0.01d0
      center_x = geometry(center_index+1, 1)
      x_min    = center_x - Lx/16.0d0
      x_max    = center_x + Lx/16.0d0
      n_points = nint((x_max - x_min) / step)

      open(1, file=trim(tmp_file_name)//"/plot.dat")

      do k = 0, n_points
        x = x_min + k * step
      
        do i = 1, atoms(center_index)%num_s_function
          gr(i) = 0.d0  
          gt(i) = 0.d0  
          do j = 1, atoms(center_index)%num_exponent_s
            exponent    = atoms(center_index)%exponent_s(j)
            coefficient = atoms(center_index)%coefficient_s(j,i)
          
            ! Regular Gaussian
            gr(i) = gr(i) + coefficient * dexp(-exponent * (x - center_x)**2)
            ! Toroidal Gaussian 
            gt(i) = gt(i) + coefficient * dexp(-exponent * Lx**2 / pi**2 * dsin( pi/Lx * (x - center_x))**2)

          end do
        end do
      
        write(1, '(100F16.8)') x, (gr(i), i=1,atoms(center_index)%num_s_function), & 
                                  (gt(i), i=1,atoms(center_index)%num_s_function)
      end do

      close(1)

      do i = 1, atoms(center_index)%num_s_function

      gr_column = 1 + i                    
      gt_column = 1 + atoms(center_index)%num_s_function + i  ! gt(i) is after all gr columns
  

      open(2, file=trim(tmp_file_name)//"/plot_s_function_"//trim(adjustl(str(i)))//".gnu")
        write(2, '(A)') 'set terminal pngcairo enhanced font "Arial,14" size 2560,1440'
        write(2, '(A)') 'set output "s_function_'//trim(adjustl(str(i)))//'.png"'
        write(2, '(A)') ''
        write(2, '(A)') '# Professional styling'
        write(2, '(A)') 'set border linewidth 3.0'
        write(2, '(A)') 'set style line 1 linecolor rgb "#0060ad" linetype 1 linewidth 4'
        write(2, '(A)') 'set style line 2 linecolor rgb "#dd181f" linetype 1 linewidth 4'
        write(2, '(A)') 'set style line 3 linecolor rgb "#000000" linetype 2 linewidth 2'
        write(2, '(A)') ''
        write(2, '(A)') '# Axis labels with better formatting'
        write(2, '(A)') 'set xlabel "x coordinate" font "Arial,32"'
        write(2, '(A)') 'set ylabel "Gaussian value" font "Arial,32"'
        write(2, '(A)') ''
        write(2, '(A)') '# Title and key'
        write(2, '(A,I0,A)') 'set title "Contracted s-Function ', i, '" font "Arial,36"'
        write(2, '(A)') 'set key top right box linestyle 3 spacing 1.5 font "Arial,24"'
        write(2, '(A)') ''
        write(2, '(A)') '# Grid for better readability'
        write(2, '(A)') 'set grid linecolor rgb "#cccccc" linetype 2 linewidth 1.0'
        write(2, '(A)') ''
        write(2, '(A)') '# Plot this specific contracted Gaussian'
        write(2, '(A,I0,A)') 'plot \'
        write(2, '(A,I0,A,I0,A)') '  "'//trim(tmp_file_name)//'/plot.dat" using 1:', gr_column, &
                                  ' with lines linestyle 1 title "Regular Gaussian (Function ', i, ')", \'
        write(2, '(A,I0,A,I0,A)') '  "'//trim(tmp_file_name)//'/plot.dat" using 1:', gt_column, &
                                  ' with lines linestyle 2 title "Toroidal Gaussian (Function ', i, ')"'
        write(2, '(A)') ''
      close(2)
        call system("gnuplot "//trim(tmp_file_name)//"/plot_s_function_"//trim(adjustl(str(i)))//".gnu")
      end do

      print *, "High-quality plotting completed!"
      print *, "Files created:"
      print *, "  - s_functions.png (s-type Gaussians)"
      print *, "  - p_functions.png (p-type Gaussians)" 
      print *, "  - "//trim(tmp_file_name)//"/plot.dat (data file)"

      deallocate(gr)
      deallocate(gt)

      if (atoms(center_index)%num_p_function > 0) then 

      open(1, file=trim(tmp_file_name)//"/plot.dat")

      do k = 0, n_points
        x = x_min + k * step
      
        do i = 1, atoms(center_index)%num_p_function
          grp(i) = 0.d0  
          gtp(i) = 0.d0  
          do j = 1, atoms(center_index)%num_exponent_p
            exponent    = atoms(center_index)%exponent_p(j)
            coefficient = atoms(center_index)%coefficient_p(j,i)
          
            ! Regular Gaussian
            grp(i) = grp(i) + coefficient * (x-center_x)                                      * dexp(-exponent * (x - center_x)**2)
            ! Toroidal Gaussian 
            gtp(i) = gtp(i) + coefficient * ((Lx/(2.d0*pi)) * dsin(2.d0*pi/Lx *(x-center_x))) * dexp(-exponent * Lx**2 / pi**2 * dsin( pi/Lx * (x - center_x))**2)

          end do
        end do
      
        write(1, '(100F16.8)') x, (grp(i), i=1,atoms(center_index)%num_p_function), & 
                                  (gtp(i), i=1,atoms(center_index)%num_p_function)
      end do

      close(1)


      do i = 1, atoms(center_index)%num_p_function

        gr_column = 1 + i                    
        gt_column = 1 + atoms(center_index)%num_p_function + i  ! gt(i) is after all gr columns
  

      open(2, file=trim(tmp_file_name)//"/plot_p_function_"//trim(adjustl(str(i)))//".gnu")
        write(2, '(A)') 'set terminal pngcairo enhanced font "Arial,14" size 2560,1440'
        write(2, '(A)') 'set output "p_function_'//trim(adjustl(str(i)))//'.png"'
        write(2, '(A)') ''
        write(2, '(A)') '# Professional styling'
        write(2, '(A)') 'set border linewidth 3.0'
        write(2, '(A)') 'set style line 1 linecolor rgb "#0060ad" linetype 1 linewidth 4'
        write(2, '(A)') 'set style line 2 linecolor rgb "#dd181f" linetype 1 linewidth 4'
        write(2, '(A)') 'set style line 3 linecolor rgb "#000000" linetype 2 linewidth 2'
        write(2, '(A)') ''
        write(2, '(A)') '# Axis labels with better formatting'
        write(2, '(A)') 'set xlabel "x coordinate" font "Arial,32"'
        write(2, '(A)') 'set ylabel "Gaussian value" font "Arial,32"'
        write(2, '(A)') ''
        write(2, '(A)') '# Title and key'
        write(2, '(A,I0,A)') 'set title "Contracted p-Function ', i, '" font "Arial,36"'
        write(2, '(A)') 'set key top right box linestyle 3 spacing 1.5 font "Arial,24"'
        write(2, '(A)') ''
        write(2, '(A)') '# Grid for better readability'
        write(2, '(A)') 'set grid linecolor rgb "#cccccc" linetype 2 linewidth 1.0'
        write(2, '(A)') ''
        write(2, '(A)') '# Plot this specific contracted Gaussian'
        write(2, '(A,I0,A)') 'plot \'
        write(2, '(A,I0,A,I0,A)') '  "'//trim(tmp_file_name)//'/plot.dat" using 1:', gr_column, &
                                  ' with lines linestyle 1 title "Regular p Gaussian (Function ', i, ')", \'
        write(2, '(A,I0,A,I0,A)') '  "'//trim(tmp_file_name)//'/plot.dat" using 1:', gt_column, &
                                  ' with lines linestyle 2 title "Toroidal p Gaussian (Function ', i, ')"'
        write(2, '(A)') ''
      close(2)
        call system("gnuplot "//trim(tmp_file_name)//"/plot_p_function_"//trim(adjustl(str(i)))//".gnu")
      end do

      end if 
        
      stop


      contains

      character(len=20) function str(o)
        integer, intent(in) :: o
        write(str, *) o
        str = adjustl(str)
      end function str

end subroutine plot