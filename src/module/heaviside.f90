MODULE HeavisideModule
  IMPLICIT NONE
  CONTAINS
    FUNCTION Heaviside(x) RESULT(H)
      IMPLICIT NONE
      double precision, INTENT(IN) :: x
      double precision :: H
      H = MERGE(1.0, 0.0, x >= 0.0)
    END FUNCTION Heaviside
END MODULE HeavisideModule