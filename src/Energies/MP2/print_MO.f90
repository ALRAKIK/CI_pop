subroutine print_MO_file(nBas,ERI)

      use files 
      implicit none 

      integer , intent(in)         ::                nBas
      double precision             :: ERI(nBas,nBas,nBas,nBas)

      integer                      :: i , j , k , l 


      open(1,file=trim(tmp_file_name)//"/ERI_MO.dat")
        do i = 1, nBas
          do j = 1 , nBas
            do k = 1 , nBas
              do l = 1 , nBas
                if (abs(ERI(i,j,k,l)) > 1e-14 ) write(1,'(4I4,f24.16)') i , j , k , l , ERI(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine


subroutine print_MO_file_SAO(nBas,ERI)

      use files 
      use constants_module
      implicit none 

      integer , intent(in)         ::                    nBas
      complex(dpc)                 :: ERI(nBas,nBas,nBas,nBas)

      integer                      :: i , j , k , l 


      open(1,file=trim(tmp_file_name)//"/ERI_MO.dat")
        do i = 1, nBas
          do j = 1 , nBas
            do k = 1 , nBas
              do l = 1 , nBas
                if (abs(ERI(i,j,k,l)) > 1e-14 ) write(1,'(4I4,f24.16,f24.16)') i , j , k , l , ERI(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine