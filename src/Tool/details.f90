subroutine details_integrals(n_f,S,T,V,ERI)

      use files 
      implicit none 

      integer , intent(in)         ::                 n_f 
      double precision             ::           S(n_f,n_f)
      double precision             ::           T(n_f,n_f)
      double precision             ::           V(n_f,n_f)
      double precision             :: ERI(n_f,n_f,n_f,n_f)


      ! local ! 

      integer                      :: i, j, k, l 


      open(1,file=trim(tmp_file_name)//"/OV.dat ")
        do i = 1 , size(S,1)
          do j = i , size(S,1)
            if (abs(S(i,j)) > 1e-15 ) write(1,*) i , j , S(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/OV_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(S,1))
      do i = 1 , size(S,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (S(i,j),j=1,size(S,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI.dat ")
      do i = 1 , size(T,1)
        do j = i , size(T,1)
          if (abs(T(i,j)) > 1e-15 ) write(1,*) i, j , T(i,j)
        end do 
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/KI_matrix.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(T,1))
        do i = 1 , size(T,1)
          write(1,'(i3,6x,1000(f16.12,2x))') i ,   (T(i,j),j=1,size(T,1))
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA.dat ")
        do i = 1 , size(V,1)
          do j = i , size(V,1)
            if (abs(V(i,j)) > 1e-15 ) write(1,*) i , j ,  V(i,j)
          end do 
        end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/NA_matrix.dat ")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(V,1))
      do i = 1 , size(V,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i , (V(i,j),j=1,size(V,1))
      end do 
      close(1)

      open(1,file=trim(tmp_file_name)//"/ERI.dat")
        !do i = 1, number_of_functions_per_unitcell
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