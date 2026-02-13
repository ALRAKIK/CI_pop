subroutine Print_the_input_file()

      use files 

      implicit none

      character(len=100)   :: lines

      write(outfile,'(2a)') "      ",repeat("=",67)
      
      Call HEADER ('The Input File',-1)

      open(1,file="unitcell.mol")
      do
        read(1,'(A)',end=3) lines
        write(outfile,'(6x,A)') trim(lines)
      end do 

3      close(1)

end subroutine Print_the_input_file