subroutine build_super_molecule(keyword)

      implicit none 
      
      !-----------------------------------------------------------------!

      ! input ! 


      ! local ! 

      integer                        :: number_of_unitcell , max_atoms , num_atoms
      integer                        :: io_stat , i , j , atom_index


      double precision,parameter     :: pi  = 3.14159265358979323846D00
      double precision               :: distance_between_unitcells , L
      double precision               :: rx , theta 
      

      double precision,allocatable   ::  geometry_unitcell(:,:)
      double precision,allocatable   ::           unitcell(:,:)
      double precision,allocatable   ::     super_geometry(:,:)


      logical                        ::           reading_cell
      character(len=10)              ::    type_of_calculation
      character(len=100)             ::                   line

      character(len=2),allocatable   ::          atom_names(:)
      character(len=2),allocatable   ::             a_names(:)
      character(len=2),allocatable   ::         super_atoms(:)


      ! output !

      character(len=10),intent(out)  ::         keyword(10)

      !-----------------------------------------------------------------!

      open(1,file="unitcell.mol")

      read(1,*) type_of_calculation , number_of_unitcell , distance_between_unitcells , L 

      if (type_of_calculation == "Ring") then
        write(*,'(a)') "Type of calculation: Ring"
      else if (type_of_calculation == "Torus") then
        write(*,'(a)') "Type of calculation: Torus"
      else if (type_of_calculation == "OBC") then 
        write(*,'(a)') "Type of calculation: OBC"
      else if (type_of_calculation == "Tori") then 
        write(*,'(a)') "Type of calculation: Toroidal"
      else if (type_of_calculation == "Tori2D") then 
        write(*,'(a)') "Type of calculation: Toroidal real 2D"
      else if (type_of_calculation == "Tori3D") then 
        write(*,'(a)') "Type of calculation: Toroidal real 3D"
      else if (type_of_calculation == "OBC2D") then 
        write(*,'(a)') "Type of calculation: OBC with FCI"
      else 
        write(*,'(a)') "Error: Unknown type of calculation. Please use either 'Torus','Ring','Tori' or 'OBC' ."
        stop
      end if

      do
        read(1, '(A)',iostat=io_stat) line
        if (io_stat /= 0) then
          write(*,'(a)') ""
          write(*,'(a)')'Error: Reached end of file without finding opening $$ marker in your unitcell'
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

      end do

      max_atoms = 100

      allocate(atom_names(max_atoms))
      allocate(geometry_unitcell(max_atoms,3))

      num_atoms = 0
      geometry_unitcell(:,3) = 0.d0 

      do
        read(1, '(A)', iostat=io_stat) line

        if (trim(line) == '$$') then
          reading_cell = .false.
          exit
        end if

        num_atoms = num_atoms + 1

        read(line, *, iostat=io_stat) atom_names(num_atoms), geometry_unitcell(num_atoms,1), &
                            geometry_unitcell(num_atoms,2), geometry_unitcell(num_atoms,3)
        if (io_stat /= 0) then
        print *, 'Error parsing atom data:', trim(line)
        num_atoms = num_atoms - 1
        end if

      end do

      close(1)

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

      allocate(super_atoms(num_atoms * number_of_unitcell))
      allocate(super_geometry(num_atoms * number_of_unitcell, 3))

      atom_index = 0

      do i = 0, number_of_unitcell-1
        do j = 1, num_atoms
            atom_index = atom_index + 1
            super_atoms(atom_index) = a_names(j)
            super_geometry(atom_index, 1) = unitcell(j, 1) + i * distance_between_unitcells
            super_geometry(atom_index, 2) = unitcell(j, 2)
            super_geometry(atom_index, 3) = unitcell(j, 3) 
        end do
      end do

      open(2,file="supermolecule.mol")
        do i = 1, num_atoms * number_of_unitcell
          write(2, *) super_atoms(i), super_geometry(i, 1), &
                                   super_geometry(i, 2), super_geometry(i, 3)
        end do
      close(2)

      open(3,file="torus_parameters.inp")
        if (keyword(4) == 'Angstrom') then 
          write(3,*) L * 1.8897261249935897D00
        else
          write(3,*) L
        end if
        !write(3,*) L
        write(3,*) num_atoms
      close(3)


      if (type_of_calculation == "Ring") then 

        atom_index = 0
        rx         = L / (2.d0*pi) 

        do i = 0, number_of_unitcell-1
          do j = 1, num_atoms
              atom_index = atom_index + 1
              super_atoms(atom_index) = a_names(j)
              theta = (i * distance_between_unitcells + unitcell(j, 1) )* 2.d0*pi / (number_of_unitcell * distance_between_unitcells)
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
            write(3,*) L * 1.8897261249935897D00
          else
            write(3,*) L
          end if
!          write(3,*) L
          write(3,*) num_atoms
        close(3)

      end if 

      open(4 , file =  "general_parameters.dat" )
        write (4,'(a10)') type_of_calculation
      close(4)

end subroutine