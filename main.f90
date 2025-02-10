program CI 

      use atom_basis

      implicit none 

      integer                      :: i , j 
      integer                      :: n_atoms

      double precision             :: geometry_tmp(100,3)
      integer                      :: charge_tmp(100)


      double precision,allocatable :: geometry(:,:)
      integer         ,allocatable :: charge(:)
      type(atom)      ,allocatable :: atoms(:)

      call read_geometry(n_atoms,charge_tmp,geometry_tmp)

      allocate(geometry(n_atoms,3))
      allocate(charge(n_atoms))

      do i = 1 , n_atoms
        charge    (i) =     charge_tmp(i)
        geometry(i,1) = geometry_tmp(i,1)
        geometry(i,2) = geometry_tmp(i,2)
        geometry(i,3) = geometry_tmp(i,3)
      end do 

      allocate(atoms(n_atoms))

      call basis(n_atoms,charge,atoms)

      do i = 1 , n_atoms
        write(*,"(I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      
!     -------------------------------------------------------------------     !
!                    Here we start calculate the integrals 
!     -------------------------------------------------------------------     !

      call overlap_matrix(n_atoms,geometry,atoms)

      call kinetic_matrix(n_atoms,geometry,atoms)

      call nuclear_attraction_matrix(n_atoms,geometry,atoms)













end program 