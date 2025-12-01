subroutine print_MO_file(n_f,HC,ERI)

      use files 
      implicit none 

      integer , intent(in)         ::                 n_f
      double precision             ::          HC(n_f,n_f)
      double precision             :: ERI(n_f,n_f,n_f,n_f)

      integer                      :: i , j , k , l 

      open(1,file=trim(tmp_file_name)//"/HC_MO_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(HC,1))
      do i = 1 , size(HC,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (HC(i,j),j=1,size(HC,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/HC_MO.dat ")
        do i = 1 , size(HC,1)
          do j = i , size(HC,1)
            if (abs(HC(i,j)) > 1e-15 ) write(1,*) i , j ,  HC(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/ERI_MO.dat")
        do i = 1, n_f
          do j = 1 , n_f
            do k = 1 , n_f
              do l = 1 , n_f
                if (abs(ERI(i,j,k,l)) > 1e-24 ) write(1,'(4I4,f24.16)') i , j , k , l , ERI(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine