subroutine basis(number_of_atoms,charge,atoms)

      use atom_basis

      implicit none 

      integer,intent(in)      :: number_of_atoms
      integer,intent(in)      :: charge(number_of_atoms)
      type(atom),intent(out)  :: atoms(number_of_atoms)

      integer                 :: i 

      call system('rm -r tmp')
      call system('mkdir tmp')

      call extract_basis   (number_of_atoms,charge)

      call normalize_basis()

      call system("python src/basis/clean_lines.py")

      do i = 1, number_of_atoms
        atoms(i)%charge = charge(i)
        call read_basis_class(atoms(i))
      end do

end subroutine basis 

subroutine basis_tor(number_of_atoms,charge,atoms,norm_helper)

      use atom_basis

      implicit none 

      integer,intent(in)        :: number_of_atoms
      integer,intent(in)        :: charge(number_of_atoms)
      type(atom),intent(out)    :: atoms(number_of_atoms)
      type(atom),intent(out)    :: norm_helper(number_of_atoms)

      integer                   :: i 

      call system('rm -r tmp')
      call system('mkdir tmp')

      call extract_basis_tor   (number_of_atoms,charge)

      call normalize_basis_tor()

      call system("python src/basis/clean_lines_tor.py")
      call system("python src/basis/clean_lines_tor_help.py")

      do i = 1, number_of_atoms
        atoms(i)%charge = charge(i)
        call read_basis_class_tor(atoms(i),1)
        norm_helper(i)%charge = charge(i)
        call read_basis_class_tor(norm_helper(i),2)
      end do

end subroutine basis_tor