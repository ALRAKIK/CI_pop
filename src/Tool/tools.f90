subroutine Title()

      use files 

      implicit none 

      character(len=8)  :: date
      character(len=10) :: time

      ! local ! 

      write(outfile,*) 
      write(outfile,*) 

      write(outfile,"(8x,a)") "__| |___________________________________________________| |__"
      write(outfile,"(8x,a)") "__   ___________________________________________________   __"
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "  | |  ______     __        ______   ______     ______  | |  "
      write(outfile,"(8x,a)") "  | | /\  ___\   /\ \      /\  == \ /\  __ \   /\  == \ | |  "
      write(outfile,"(8x,a)") "  | | \ \ \____  \ \ \     \ \  _-/ \ \ \/\ \  \ \  _-/ | |  "
      write(outfile,"(8x,a)") "  | |  \ \_____\  \ \_\     \ \_\    \ \_____\  \ \_\   | |  "
      write(outfile,"(8x,a)") "  | |   \/_____/   \/_/      \/_/     \/_____/   \/_/   | |  "
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "__| |___________________________________________________| |__"
      write(outfile,"(8x,a)") "__   ___________________________________________________   __"
      write(outfile,"(8x,a)") "  | |                                                   | |  "

      write(outfile,'(a)')
      write(outfile,'(a)')
      write(outfile,'(a)')  "      CI_POP: Calculate Integrals, is a program to calculate   the" 
      write(outfile,'(a)')  "              molecular integrals using Clifford Toroidal Gaussian"
      write(outfile,'(a)')  "              functions for advanced quantum chemistry computations"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(2a)') "      ", "Amer Alrakik        Email  : alrakikamer@gmail.com"
      write(outfile,'(2a)') "      ", "                    Github : https://github.com/ALRAKIK/CI_pop"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ", "Contributors        Arjan Berger"
      write(outfile,'(2a)') "      ", "                    Stefano Evangelisti"
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  "                          Version: 0.0.2 | Quantum Chemistry | 2025"
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  ""


      call date_and_time(date, time)

      write(outfile,'(6x,10a)') "Executing Date : ", date(7:8), "/", date(5:6), "/", date(1:4)
      write(outfile,'(6x,10a)') "Executing time : ", time(1:2), ":", time(3:4), ":", time(5:6)


end subroutine Title 

subroutine header(HEAD,IN)

      use files 
      implicit none 

      CHARACTER :: HEAD*(*)
      integer   :: in , indent , length , I 
  
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
        INDENT = IN + 1
      ELSE
        INDENT = (72 - LENGTH)/2 + 1
      END IF
      
      WRITE (outfile, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (outfile, '(80A)') (' ',I=1,INDENT), ('=',I=1,LENGTH)
      WRITE (outfile, '()')
      FLUSH(outfile)
  
end subroutine Header

subroutine header_under(HEAD,IN)

      use files
      implicit none 

      CHARACTER :: HEAD*(*)
      integer   :: in , indent , length , I 
  
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
        INDENT = IN + 1
      ELSE
        INDENT = (72 - LENGTH)/2 + 1
      END IF
      
      WRITE (outfile, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (outfile, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (outfile, '()')
    
end subroutine Header_under

subroutine header_HF(HEAD,IN)

      use files
      implicit none 

      CHARACTER :: HEAD*(*)
      integer   :: in , indent , length , I 
  
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
        INDENT = IN + 1
      ELSE
        INDENT = (72 - LENGTH)/2 + 1
      END IF
      
      WRITE (HFfile, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (HFfile, '(80A)') (' ',I=1,INDENT), ('_',I=1,LENGTH)
      WRITE (HFfile, '()')
    
end subroutine Header_HF

subroutine print_orbital_table(ERI,number_of_functions)

      use files 
      use classification_ERI

      implicit none

      type(ERI_function)  :: ERI(number_of_functions)
      integer, intent(in) :: number_of_functions
      integer             :: i, j
        
      call header_under("Atomic Orbitals", -1)
        
      ! Print header
      write(outfile,'(A)') "|-------------------------------------------------------------------------------|"
      write(outfile,'(A)') "|  Idx         X         Y        Z      Type        Exponent      Coefficient  |"
      write(outfile,'(A)') "|-------------------------------------------------------------------------------|"

      do i = 1 , number_of_functions
        write(outfile,'(a1,I5,3x,3f10.6,3x,a2,6x,f12.6,2x,f12.6,4x,a1)')"|", i  , ERI(i)%x , ERI(i)%y , ERI(i)%z , ERI(i)%orbital , ERI(i)%exponent(1) , ERI(i)%coefficient(1) , "|"
        do j = 2 , size(ERI(i)%exponent)
          write(outfile,'(a1,49x,f12.6,2x,f12.6,4x,a1)')"|" , ERI(i)%exponent(j) , ERI(i)%coefficient(j) , "|"
        end do
        write(outfile,'(A)') "|-------------------------------------------------------------------------------|"
      end do

        write(outfile,'(A)') ""

        FLUSH(outfile)

end subroutine print_orbital_table