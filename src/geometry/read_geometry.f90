subroutine read_geometry(number_of_atoms,charge,geometry,calculation_type,type)

      use constants_module

      implicit none 


      integer                       :: i 
      integer                       :: number_of_atoms
      double precision              :: geometry(max_atom,3)
      integer                       :: charge(max_atom)
      character(len=10),intent(out) :: calculation_type
      character(len=2)              :: type(max_atom)
      
      
      open(1,file="supermolecule.mol")

      number_of_atoms = 0

      do
        read(1,"(A)",end=3)
        number_of_atoms = number_of_atoms + 1
      end do 

3     close(1)


      open(1,file="supermolecule.mol")
      do i = 1 , number_of_atoms
        read(1,*) type(i) , geometry(i,1) , geometry(i,2) , geometry(i,3)
      end do 
      close(1)


      call system("rm supermolecule.mol")

      open(4, file = "general_parameters.dat")
        read (4,*) calculation_type
      close(4)

      call system("rm general_parameters.dat")

      call get_charge(number_of_atoms,type,charge)

end subroutine