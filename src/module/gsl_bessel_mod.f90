module gsl_bessel_mod
  use, intrinsic :: iso_c_binding
  implicit none
  public :: bessi_scaled

      interface
      ! Scaled version: returns exp(-|x|) * I_n(x)
      function gsl_sf_bessel_In_scaled(n, x) bind(C, name="gsl_sf_bessel_In_scaled")
       use, intrinsic :: iso_c_binding
       integer(c_int),  value :: n
       real(c_double),  value :: x
       real(c_double)         :: gsl_sf_bessel_In_scaled
      end function gsl_sf_bessel_In_scaled
     
      ! Normal (unscaled) version: returns I_n(x)
      function gsl_sf_bessel_In(n, x) bind(C, name="gsl_sf_bessel_In")
       use, intrinsic :: iso_c_binding
       integer(c_int),  value :: n
       real(c_double),  value :: x
       real(c_double)         :: gsl_sf_bessel_In
      end function gsl_sf_bessel_In

      function gsl_sf_bessel_In_scaled_array(nmin, nmax, x, result_array) &
               bind(C, name="gsl_sf_bessel_In_scaled_array")
       use, intrinsic :: iso_c_binding
       integer(c_int), value :: nmin
       integer(c_int), value :: nmax
       real(c_double), value :: x
       type(c_ptr)   , value :: result_array
       integer(c_int)        :: gsl_sf_bessel_In_scaled_array
      end function gsl_sf_bessel_In_scaled_array


      end interface

contains

      ! Scaled Bessel function wrapper
      function bessi_scaled(n, x) result(res)

      use, intrinsic :: iso_c_binding
      integer, intent(in)        :: n
      real(c_double), intent(in) :: x
      real(c_double)             :: res
      real(c_double)             :: pi     = 3.14159265358979323846D00

      res = gsl_sf_bessel_In_scaled(int(n, c_int), x)

      end function bessi_scaled

      function bessi(n, x) result(res)

      use, intrinsic :: iso_c_binding
      integer, intent(in)        :: n
      real(c_double), intent(in) :: x
      real(c_double)             :: res
      real(c_double)             :: pi     = 3.14159265358979323846D00

      res = gsl_sf_bessel_In(int(n, c_int), x)

      end function bessi

      function bessi_scaled_array(nmin, nmax, x, result_array) result(ierr)
      use, intrinsic :: iso_c_binding
      integer, intent(in) :: nmin
      integer, intent(in) :: nmax
      real(c_double), intent(in) :: x
      real(c_double), intent(out), target :: result_array(nmin:nmax)
      integer :: ierr
      
      ierr = gsl_sf_bessel_In_scaled_array( &
             int(nmin, c_int), &
             int(nmax, c_int), &
             x, &
             c_loc(result_array) )
             
      end function bessi_scaled_array
  
end module gsl_bessel_mod
