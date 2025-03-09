module filter_module

  implicit none

contains

      subroutine filter(array_coef,array_expo,f_array_coef,f_array_expo)
        
        implicit none 
        
        double precision, dimension(:), intent(in) :: array_coef , array_expo

        double precision, dimension(:), allocatable, intent(out) :: f_array_coef, f_array_expo

        integer, allocatable        :: zero_indices(:)

        logical                     :: mask(size(array_coef))

        integer                     :: i, n , count_nonzeros , num_zeros

        double precision, parameter :: THRESHOLD = 1.0e-15

        
        n = size(array_coef)

        if (size(array_expo) /= n) then
          write(*,*) "ERROR: Arrays must have the same size!"
          return
        end if

        mask = (abs(array_coef) < THRESHOLD)

        num_zeros = count(mask)

        count_nonzeros = n - num_zeros

        allocate(f_array_coef(count_nonzeros))
        allocate(f_array_expo(count_nonzeros))

        allocate(zero_indices(num_zeros))

        zero_indices = pack([(i, i=1,n)], mask)

        f_array_coef = pack(array_coef, .not. mask)
        f_array_expo = pack(array_expo, .not. mask)
        
      end subroutine filter
end module filter_module