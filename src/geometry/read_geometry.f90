subroutine read_geometry(number_of_atoms,charge,geometry)

      implicit none 


      integer                      :: i 
      integer                      :: number_of_atoms
      double precision             :: geometry(100,3)
      integer                      :: charge(100)
      character(len=2),allocatable :: type(:)
      
      
!      open(1,file="mol.mol")
      open(1,file="supermolecule.mol")

      number_of_atoms = 0

      do
        read(1,"(A)",end=3)
        number_of_atoms = number_of_atoms + 1
      end do 


3     close(1)


      allocate(type(number_of_atoms))      

!      open(1,file="mol.mol")
      open(1,file="supermolecule.mol")
      do i = 1 , number_of_atoms
        read(1,*) type(i) , geometry(i,1) , geometry(i,2) , geometry(i,3)
      end do 
      close(1)

      call get_charge(number_of_atoms,type,charge)



end subroutine