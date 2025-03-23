function Kronecker_delta(i) result(delta)

      ! Kronecker Delta

      implicit none

      ! Input variables

      integer,intent(in)            :: i

      ! Output variables

      double precision              :: delta

      if(i == 0) then
        delta = 1d0
      else
        delta = 0d0
      endif

end function Kronecker_delta