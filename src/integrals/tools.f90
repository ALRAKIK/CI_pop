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

end subroutine print_orbital_table