subroutine check_symmetric_matrix(nBas,S,T,V,HC)

      implicit none 

      integer, intent(in) :: nBas
      integer             :: i , j

      double precision, intent(in) :: S (nBas,nBas)
      double precision, intent(in) :: T (nBas,nBas)
      double precision, intent(in) :: V (nBas,nBas)
      double precision, intent(in) :: HC(nBas,nBas) 


      do i = 1 , nBas-1
        do j = i , nBas 
          if (abs(S(i,j)-S(j,i)) > 1.0d-10) then 
            print * , "S is not symmetric"
            stop 
          end if 
          if (abs(T(i,j)-T(j,i)) > 1.0d-10) then 
            print * , "T is not symmetric"
            stop 
          end if 
          if (abs(V(i,j)-V(j,i)) > 1.0d-10) then 
            print * , "V is not symmetric"
            stop 
          end if 
          if (abs(HC(i,j)-HC(j,i)) > 1.0d-10) then 
            print * , "HC is not symmetric"
            stop 
          end if
        end do 
      end do 

end subroutine check_symmetric_matrix