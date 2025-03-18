function Kronecker_delta(i,j) result(delta)

  ! Kronecker Delta

  implicit none

  ! Input variables

  integer,intent(in)            :: i,j

  ! Output variables

  double precision              :: delta

  if(i == j) then
    delta = 1d0
  else
    delta = 0d0
  endif

end function Kronecker_delta