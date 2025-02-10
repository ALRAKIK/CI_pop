function boys_function(n, T) result(Fn)

      implicit none
      integer, intent(in)          :: n
      double precision, intent(in) :: T
      double precision             :: Fn, exp_T, sqrt_T, erf_T
      integer                      :: k
      double precision,parameter   :: pi = dacos(-1.0d0)

      ! Base case for F_0(T)
      sqrt_T = sqrt(T)
      erf_T = erf(sqrt_T)  ! Built-in error function
      Fn = 0.5d0 * sqrt(pi) * erf_T / sqrt_T

      ! Recurrence for higher n
      exp_T = exp(-T)
      do k = 1, n
          Fn = (2.0d0*k - 1.0d0)*Fn + exp_T
          Fn = Fn / (2.0d0 * T)
      end do

      if (T == 0.0d0) then
          Fn = 1.0d0 / (2.0d0 * n + 1.0d0)
          return
      end if

end function boys_function
