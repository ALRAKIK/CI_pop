module gsl_bessel_mod
  use, intrinsic :: iso_c_binding
  implicit none
  public :: bessi_scaled

  interface
     function gsl_sf_bessel_In_scaled(n, x) bind(C, name="gsl_sf_bessel_In_scaled")
       use, intrinsic :: iso_c_binding
       integer(c_int),  value :: n
       real(c_double),  value :: x
       real(c_double)         :: gsl_sf_bessel_In_scaled
     end function gsl_sf_bessel_In_scaled
  end interface

contains
  function bessi_scaled(n, x) result(res)
    use, intrinsic :: iso_c_binding
    integer, intent(in)      :: n
    real(c_double), intent(in) :: x
    real(c_double)           :: res

    res = gsl_sf_bessel_In_scaled(int(n, c_int), x)
  end function bessi_scaled
end module gsl_bessel_mod
