subroutine build_super_molecule(keyword,num_atom_per_unitcell)

      use constants_module
      use unitcell_module

      implicit none 
      
      !-----------------------------------------------------------------!

      ! local ! 

      integer                        :: number_of_unitcell
      integer                        :: num_atoms
      integer                        :: io_stat, i, j, k, l, atom_index

      double precision               :: distance_between_unitcells
      double precision               :: Lx , Ly , Lz
      double precision               :: rx, theta

      double precision,allocatable   ::  geometry_unitcell(:,:)
      double precision,allocatable   ::           unitcell(:,:)
      double precision,allocatable   ::     super_geometry(:,:)

      character(len=10)              ::    type_of_calculation
      character(len=100)             ::                   line

      character(len=2),allocatable   ::          atom_names(:)
      character(len=2),allocatable   ::             a_names(:)
      character(len=2),allocatable   ::         super_atoms(:)

      logical                        ::           reading_cell


      ! output !

      character(len=10),intent(out)  ::            keyword(20)
      integer                        :: num_atom_per_unitcell

      !-----------------------------------------------------------------!

      open(1,file="unitcell.mol")

      read(1,*) type_of_calculation

      call read_2nd_line(type_of_calculation, number_of_unitcell,  &
                        & distance_between_unitcells, Lx , Ly , Lz)

      call read_keywords(keyword)

      allocate(atom_names(max_atom))
      allocate(geometry_unitcell(max_atom,3))

      num_atoms              = 0
      geometry_unitcell(:,3) = 0.d0 

      ! --------------------------------------------------------------- !
      !               Read unitcell geometry from file                  !
      ! --------------------------------------------------------------- !

      do

        read(1, '(A)', iostat=io_stat) line
        if (trim(line) == '$$') then
          reading_cell = .false.
          exit
        end if

        num_atoms = num_atoms + 1

        read(line, *, iostat=io_stat) atom_names(num_atoms),            &
        &                             geometry_unitcell(num_atoms,1),   &
        &                             geometry_unitcell(num_atoms,2),   &
        &                             geometry_unitcell(num_atoms,3)

        if (io_stat /= 0) then
          print *, 'Error parsing atom data:', trim(line)
          num_atoms = num_atoms - 1
        end if

      end do

      close(1)

      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !


      num_atom_per_unitcell = num_atoms

      allocate(unitcell(num_atoms,3))
      allocate(a_names(num_atoms))

      do i = 1 , num_atoms
        unitcell(i,1) = geometry_unitcell(i,1)
        unitcell(i,2) = geometry_unitcell(i,2)
        unitcell(i,3) = geometry_unitcell(i,3)
        a_names(i)    = atom_names(i)
      end do

      deallocate(atom_names)
      deallocate(geometry_unitcell)
      ! --------------------------------------------------------------- !

      if (type_of_calculation == 'Tori1D') then 

      allocate(super_atoms(num_atoms * number_of_unitcell))
      allocate(super_geometry(num_atoms * number_of_unitcell, 3))

      atom_index = 0

      do i = 0, number_of_unitcell-1
        do j = 1, num_atoms
            atom_index = atom_index + 1
            super_atoms(atom_index) = a_names(j)
            ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
            super_geometry(atom_index, 1) = unitcell(j, 1)              &
          &                           + i * distance_between_unitcells
            ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
            super_geometry(atom_index, 2) = unitcell(j, 2)
            super_geometry(atom_index, 3) = unitcell(j, 3) 
        end do
      end do

      open(2,file="supermolecule.mol")
        do i = 1, num_atoms * number_of_unitcell
          write(2, *)       super_atoms(i), super_geometry(i, 1),       &
          &           super_geometry(i, 2), super_geometry(i, 3)
        end do
      close(2)

      open(3,file="torus_parameters.inp")
        if (keyword(4) == 'Angstrom') then 
          write(3,*) Lx * Ang_par , Ly * Ang_par , Lz * Ang_par 
        else
          write(3,*) Lx , Ly , Lz
        end if
        write(3,*) num_atoms
      close(3)

      end if 



      if (type_of_calculation == "Ring") then 

        allocate(super_atoms(num_atoms * number_of_unitcell))
        allocate(super_geometry(num_atoms * number_of_unitcell, 3))

        atom_index = 0
        rx         = Lx / (2.d0*pi) 

        do i = 0, number_of_unitcell-1
          do j = 1, num_atoms
              atom_index = atom_index + 1
              super_atoms(atom_index) = a_names(j)
              theta = (i * distance_between_unitcells + unitcell(j, 1)) &
            & * 2.d0*pi/(number_of_unitcell * distance_between_unitcells)
              super_geometry(atom_index, 1) = rx * cos(theta)
              super_geometry(atom_index, 2) = rx * sin(theta)
              super_geometry(atom_index, 3) = 0.d0
          end do
        end do

        open(2,file="supermolecule.mol")
          do i = 1, num_atoms * number_of_unitcell
            write(2,*) super_atoms(i), super_geometry(i, 1), &
                       super_geometry(i, 2), super_geometry(i, 3)
          end do
        close(2)

        open(3,file="torus_parameters.inp")
          if (keyword(4) == 'Angstrom') then 
            write(3,*) Lx * Ang_par , Ly * Ang_par , Lz * Ang_par
          else
            write(3,*) Lx , Ly , Lz
          end if
          write(3,*) num_atoms
        close(3)

      end if 



      if (type_of_calculation == 'OBC') then 

      allocate(super_atoms(num_atoms))
      allocate(super_geometry(num_atoms, 3))

      atom_index = 0

      do j = 1, num_atoms
          atom_index = atom_index + 1
          super_atoms(atom_index) = a_names(j)
          super_geometry(atom_index, 1) = unitcell(j, 1)
          super_geometry(atom_index, 2) = unitcell(j, 2)
          super_geometry(atom_index, 3) = unitcell(j, 3) 
      end do

      open(2,file="supermolecule.mol")
        do i = 1, num_atoms * number_of_unitcell
          write(2, *)       super_atoms(i), super_geometry(i, 1),       &
          &           super_geometry(i, 2), super_geometry(i, 3)
        end do
      close(2)

      open(3,file="torus_parameters.inp")
        if (keyword(4) == 'Angstrom') then 
          write(3,*) Lx * Ang_par , Ly * Ang_par , Lz * Ang_par
        else
          write(3,*) Lx , Ly , Lz
        end if
        write(3,*) num_atoms
      close(3)

      end if 


      if (type_of_calculation == 'Tori3D') then 

      allocate(super_atoms(num_atoms * nx * ny * nz))
      allocate(super_geometry(num_atoms * nx * ny * nz, 3))

      atom_index = 0

      do i = 0, nx - 1 
        do j = 0, ny - 1 
          do k = 0 , nz -1 
            do l = 1 , num_atoms
            atom_index = atom_index + 1
            super_atoms(atom_index) = a_names(l)
            ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
            super_geometry(atom_index, 1) = unitcell(l, 1)              &
          &                           + i * distance_between_unitcells
            ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
            super_geometry(atom_index, 2) = unitcell(l, 2)              &
          &                           + j * distance_between_unitcells
          ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - ! - !
            super_geometry(atom_index, 3) = unitcell(l, 3)            & 
          &                           + k * distance_between_unitcells
            end do 
          end do 
        end do
      end do

      open(2,file="supermolecule.mol")
        do i = 1, num_atoms * nx * ny * nz
          write(2, *)       super_atoms(i), super_geometry(i, 1),       &
          &           super_geometry(i, 2), super_geometry(i, 3)
        end do
      close(2)

      open(3,file="torus_parameters.inp")
        if (keyword(4) == 'Angstrom') then 
          write(3,*) Lx * Ang_par , Ly * Ang_par , Lz * Ang_par
        else
          write(3,*) Lx , Ly , Lz
        end if
        write(3,*) num_atoms
      close(3)

      end if 



      ! --------------------------------------------------------------- !

      deallocate(super_atoms)
      deallocate(super_geometry)

      open(4 , file =  "general_parameters.dat" )
        write (4,'(a10)') type_of_calculation
      close(4)

      ! --------------------------------------------------------------- !

end subroutine


subroutine read_2nd_line(type_of_calculation, number_of_unitcell,     &
                         & distance_between_unitcells, Lx , Ly , Lz)
      
      use unitcell_module
      implicit none

      character(len=10)              :: type_of_calculation
      integer                        :: number_of_unitcell
      double precision               :: distance_between_unitcells
      double precision               :: Lx , Ly , Lz

        if (type_of_calculation == "Ring") then

        write(*,'(a)') "Type of calculation: Ring"
        read(1,*) number_of_unitcell , distance_between_unitcells, Lx

      else if (type_of_calculation == "Torus") then

        write(*,'(a)') "Type of calculation: Torus"
        read(1,*) number_of_unitcell , distance_between_unitcells, Lx , Ly , Lz 

      else if (type_of_calculation == "OBC") then 

        write(*,'(a)') "Type of calculation: OBC"

      else if (type_of_calculation == "Tori1D") then 

        write(*,'(a)') "Type of calculation: Toroidal"
        read(1,*) number_of_unitcell , distance_between_unitcells, Lx

      else if (type_of_calculation == "Tori2D") then 

        write(*,'(a)') "Type of calculation: Toroidal real 2D"
        read(1,*) number_of_unitcell , distance_between_unitcells, Lx , Ly

      else if (type_of_calculation == "Tori3D") then 

        write(*,'(a)') "Type of calculation: Toroidal real 3D"
        read(1,*) nx , ny , nz , distance_between_unitcells, Lx , Ly , Lz

      else if (type_of_calculation == "OBC2D") then 

        write(*,'(a)') "Type of calculation: OBC with FCI"
        read(1,*) number_of_unitcell , distance_between_unitcells

      else 
        write(*,'(a)') ""
        write(*,'(a)') "Error: Unknown type of calculation"
        write(*,'(a)') ""
        write(*,'(a)') "Please use either    'OBC',   'Ring',  'Torus'"
        write(*,'(a)') "                  'Tori1D', 'Tori2D', 'Tori3D'."
        write(*,'(a)') ""
        stop
      end if

end subroutine


subroutine read_keywords(keyword)

      implicit none

      character(len=100)             :: line
      integer                        :: io_stat

      character(len=10),intent(out)  :: keyword(20)
      logical                        :: reading_cell

      do
        read(1, '(A)',iostat=io_stat) line
        if (io_stat /= 0) then
          write(*,'(a)') ""
          write(*,'(a)')'Error: Reached end of file without finding'
          write(*,'(a)')'       opening $$ marker in your unitcell'
          write(*,'(a)') ""
          close(1)
          stop
        end if

        if (trim(line) == '$$') then
          reading_cell = .true.
          exit
        end if

        if (trim(line) == 'Read'     )   keyword(1) = 'Read'

        if (trim(line) == 'Integrals')   keyword(2) = 'Integrals'

        if (trim(line) == 'Trexio'   )   keyword(3) = 'Trexio'

        if (trim(line) == 'Angstrom' )   keyword(4) = 'Angstrom'

        if (trim(line) == 'Plot'     )   keyword(5) = 'Plot'

        if (trim(line) == 'Details'  )   keyword(6) = 'Details'

        if (trim(line) == 'ERI_a'  )     keyword(7) = 'ERI_a'

      end do

end subroutine