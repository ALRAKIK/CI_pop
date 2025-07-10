subroutine initialize_ff(calculation_type,n_atoms)

      use files

      implicit none

      ! --------------------------------------------------------------- !

      ! input ! 

      character(len=10) ,intent(in)     :: calculation_type 
      integer           ,intent(in)     :: n_atoms
      
      ! local ! 

      character(len=100)                :: command

      ! --------------------------------------------------------------- !


      write(tmp_file_name,'(A,A,A,I0,A)') "tmp_",trim(calculation_type),"_",n_atoms,"/"
      write(command,'(A,A)') 'rm -r  ' , trim(tmp_file_name)
      call system(command)
      write(command,'(A,A)') 'mkdir  ' , trim(tmp_file_name)
      call system(command)
      write(output_file_name,'(A,A,A,I0,A)') "results_",trim(calculation_type),"_",n_atoms,".out"
      open (outfile,file=trim(output_file_name))


end subroutine initialize_ff