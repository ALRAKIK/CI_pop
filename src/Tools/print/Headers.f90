subroutine header(HEAD,IN)

      use files 

      implicit none 

      character :: HEAD*(*)
      integer   :: in , indent , length , I 

      ! --------------------------------------------------------------- !

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

      ! --------------------------------------------------------------- !
  
end subroutine Header

subroutine header_under(HEAD,IN)

      use files

      implicit none 

      character :: HEAD*(*)
      integer   :: in , indent , length , I 

      ! --------------------------------------------------------------- !
  
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
        INDENT = IN + 1
      ELSE
        INDENT = (72 - LENGTH)/2 + 1
      END IF
      
      WRITE (outfile, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (outfile, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (outfile, '()')
      FLUSH(outfile)

      ! --------------------------------------------------------------- !
    
end subroutine Header_under

subroutine header_HF(HEAD,IN)

      use files

      implicit none 

      character :: HEAD*(*)
      integer   :: in , indent , length , I 

      ! --------------------------------------------------------------- !
  
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
        INDENT = IN + 1
      ELSE
        INDENT = (72 - LENGTH)/2 + 1
      END IF
      
      WRITE (HFfile, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (HFfile, '(80A)') (' ',I=1,INDENT), ('_',I=1,LENGTH)
      WRITE (HFfile, '()')
      FLUSH(outfile)

      ! --------------------------------------------------------------- !
    
end subroutine Header_HF