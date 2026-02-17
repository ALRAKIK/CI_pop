program CI 

      use files
      use constants_module
      use torus_init
      use atom_basis
      use classification_ERI
      use trexio
      use keywords
      use table_lookup_module

      implicit none

      !-----------------------------------------------------------------!
      !                         The  variables                          !
      !-----------------------------------------------------------------!


      integer                         ::                  i ,    j
      integer                         ::            n_atoms , nBAS
      integer                         ::           nO , n_electron 

      double precision                ::  geometry_tmp(max_atom,3)
      integer                         ::      charge_tmp(max_atom)
      character*(2)                   ::       label_tmp(max_atom)
      integer                         ::           n_atom_unitcell 
      integer                         ::       number_of_functions
      integer                         ::      number_of_primitives
      integer                         ::          number_of_shells
      character*(10)                  ::               keyword(20)

      double precision  ,allocatable  ::             geometry(:,:)
      integer           ,allocatable  ::                 charge(:)
      type(atom)        ,allocatable  ::                  atoms(:)
      type(atom)        ,allocatable  ::         norm_helper_px(:)
      type(atom)        ,allocatable  ::         norm_helper_py(:)
      type(atom)        ,allocatable  ::         norm_helper_pz(:)
      type(ERI_function),allocatable  ::                    AO (:)

      double precision,allocatable    ::                    S(:,:)
      double precision,allocatable    ::                    T(:,:)
      double precision,allocatable    ::                    V(:,:)
      double precision,allocatable    ::                   Hc(:,:)
      double precision,allocatable    ::                Hc_MO(:,:)
      double precision,allocatable    ::                    X(:,:)
      double precision,allocatable    ::              ERI(:,:,:,:)
      double precision,allocatable    ::           ERI_MO(:,:,:,:)
      double precision,allocatable    ::                      e(:)
      double precision,allocatable    ::                    c(:,:)

      double precision                ::               E_nuc , EHF
      double precision                ::            start,end,time

      character(len=10)               ::          calculation_type 
      character(len=2),allocatable    ::                  label(:)

      integer                         :: io_stat
      character(len=100)              :: line
      integer                         :: n_alpha , n_beta

      !-----------------------------------------------------------------!
      !                        END variables                            !
      !-----------------------------------------------------------------!

      ! --------------------------------------------------------------- !
      !                   build the super molecule                      !
      ! --------------------------------------------------------------- !
      
      call build_super_molecule(keyword,n_atom_unitcell) ! build the super molecule from the unitcell file 

      ! --------------------------------------------------------------- !
      !                        Read key words                           !
      ! --------------------------------------------------------------- !

      c_integral = any(keyword == 'Integrals')
      c_read     = any(keyword == 'Read'     )
      c_trexio   = any(keyword == 'Trexio'   )
      c_Angstrom = any(keyword == 'Angstrom' )
      c_plot     = any(keyword == 'Plot'     )
      c_details  = any(keyword == 'Details'  )
      c_MO       = any(keyword == 'MO'       )
      c_UHF      = any(keyword == 'UHF'      )
      c_Huckel   = any(keyword == 'Huckel'   )
      c_MP2      = any(keyword == 'MP2'      )
      c_OV       = any(keyword == 'OV'       )
      c_K        = any(keyword == 'K'        )
      c_NA       = any(keyword == 'NA'       )
      c_ERI      = any(keyword == 'ERI'      )
      c_One      = any(keyword == 'One'      )


      ! --------------------------------------------------------------- !
      ! *************************************************************** !

      call load_table('table.bin')
      !call print_table_info()


      if (c_UHF) then 
        n_alpha = 0 
        n_beta  = 0 
        open(1,file="unitcell.mol")
        do
        read(1, '(A)',iostat=io_stat) line
        if (io_stat /= 0) then
          write(*,'(a)') ""
          write(*,'(a)')'Error: Reached end of file without finding'
          write(*,'(a)')'              UHF keyword in your unitcell'
          write(*,'(a)') ""
          close(1)
          stop
        end if

        if (trim(line) == 'UHF') then
          read(1, *,iostat=io_stat) n_alpha , n_beta
          exit
        end if

        end do
        close(1)
      end if 

      ! --------------------------------------------------------------- !
      ! *************************************************************** !

      ! --------------------------------------------------------------- !
      ! Read the Geometry, Label and the Charge from the super molecule !
      ! --------------------------------------------------------------- !

      call read_geometry(n_atoms,charge_tmp,geometry_tmp,               &
      &                  calculation_type,label_tmp)                      

      ! --------------------------------------------------------------- !
      !           prepare the folder to write the output files          !
      ! --------------------------------------------------------------- !

      call initialize_ff(calculation_type,n_atom_unitcell,label_tmp,n_atoms)

      ! --------------------------------------------------------------- !
      !       If the calculation is on a torus read the parameters      !  
      ! --------------------------------------------------------------- !
       
      if (calculation_type ==  "Torus" .or.                             & 
      &   calculation_type == "Tori1D" .or.                             &                            
      &   calculation_type == "Tori2D" .or.                             & 
      &   calculation_type == "Tori3D")                                 &
      &   call Torus_def()

      ! --------------------------------------------------------------- !

      allocate(geometry    (n_atoms,3))
      allocate(charge        (n_atoms))
      allocate(label         (n_atoms))
      allocate(atoms         (n_atoms))
      allocate(norm_helper_px(n_atoms))
      allocate(norm_helper_py(n_atoms))
      allocate(norm_helper_pz(n_atoms))

      ! --------------------------------------------------------------- !
      !        Convert the geometry to Bohr if in Angstrom              !
      ! --------------------------------------------------------------- !

      do i = 1 , n_atoms
        label(i)      =      label_tmp(i)
        charge    (i) =     charge_tmp(i)
        if (c_Angstrom) then 
          geometry(i,1) = geometry_tmp(i,1) * Ang_par
          geometry(i,2) = geometry_tmp(i,2) * Ang_par
          geometry(i,3) = geometry_tmp(i,3) * Ang_par
        else
          geometry(i,1) = geometry_tmp(i,1)
          geometry(i,2) = geometry_tmp(i,2)
          geometry(i,3) = geometry_tmp(i,3)
        end if
      end do 



      ! --------------------------------------------------------------- !
      !                        Start Printing                           !
      ! --------------------------------------------------------------- !

      Call Title()

      Call Print_the_input_file() 

      call print_basis(n_atoms,charge)

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !             Print the geometry of the supermolecule             !
      ! --------------------------------------------------------------- !

      Call HEADER ('The Geometry',-1)

      do i = 1 , n_atoms
        write(outfile,"(6x,I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      write(outfile,*)
      write(outfile,*)

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !                  build the atomic basis set                     !
      ! --------------------------------------------------------------- !

      if (calculation_type == "Tori1D" .or. &
      &   calculation_type == "Tori2D" .or. &
          calculation_type == "Tori3D" ) then 

        call basis_tor(n_atoms,charge,atoms,norm_helper_px,norm_helper_py,norm_helper_pz,calculation_type)

      else

        call basis(n_atoms,charge,atoms)

      end if 

      ! --------------------------------------------------------------- !
      !                 Information from the Basis set                  !
      ! --------------------------------------------------------------- !

      n_electron = 0 
      do i = 1 , n_atoms
        n_electron = n_electron + atoms(i)%charge
      end do 

      nO = n_electron/2

      number_of_functions = 0 
      do i = 1 , n_atoms
        number_of_functions = number_of_functions                       &
        &  + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do

      nBAS = 0 
      do i = 1 , n_atoms
        nBAS = nBAS                                                     & 
        &  + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do

      number_of_primitives = 0
      do i = 1 , n_atoms
        number_of_primitives = number_of_primitives                     &
        &  + atoms(i)%num_exponent_s + atoms(i)%num_exponent_p
      end do

      number_of_shells = 0
      do i = 1 , n_atoms
        number_of_shells = number_of_shells                             &
        &  + atoms(i)%num_s_function + atoms(i)%num_p_function
      end do

      ! --------------------------------------------------------------- !
      !          Allocate the memory for the atomic orbitals            !
      ! --------------------------------------------------------------- !

      allocate(AO(number_of_functions))

      ! --------------------------------------------------------------- !
      !          Classify the atomic orbitals and print them            !
      ! --------------------------------------------------------------- !

      if (calculation_type == "Tori1D" .or.                             &
      &   calculation_type == "Tori2D" .or.                             & 
      &   calculation_type == "Tori3D") then 
        call classification_orbital_tor(n_atoms,number_of_functions,    &
        &                               geometry,atoms,norm_helper_px,norm_helper_py,norm_helper_pz,AO)
      else
        call classification_orbital(n_atoms,number_of_functions,        &
        &                               geometry,atoms,AO)
      end if

      call print_orbital_table(AO,number_of_functions)


      ! --------------------------------------------------------------- !
      !           print  Informations from the Basis set                !
      ! --------------------------------------------------------------- !

      call basis_info(number_of_functions,number_of_primitives,number_of_shells,n_atom_unitcell)

      ! --------------------------------------------------------------- !
      !                    Nuclear repulsion energy                     !
      ! --------------------------------------------------------------- !

      call NRE(calculation_type,n_atoms,geometry,atoms,E_nuc)

      ! --------------------------------------------------------------- !
      !             Write the geometry to a TREXIO file                 !
      ! --------------------------------------------------------------- !


      if (c_trexio) then
        call trexio_conv_init(calculation_type,n_atom_unitcell,label_tmp,n_atoms)
        call trexio_conv_global(n_atoms,label,geometry,charge,          &
        &                       E_nuc,n_electron,number_of_functions,   &
        &                       number_of_primitives,number_of_shells)
      end if


      ! --------------------------------------------------------------- !

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !               ---------------------------------                 !
      !                       Allocate the memory                       !
      !               ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas))
      
      allocate(ERI  (nBas,nBas,nBas,nBas))
      
      allocate(X(nBas,nBas),e(nBas),c(nBas,nBas))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !               ---------------------------------                 !
      !                       Plot the gussians                         !
      !               ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (c_plot) then 
        call plot(n_atoms,geometry,atoms)
      end if 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !               ---------------------------------                 !
      !                   calculate the integrals                       !
      !               ---------------------------------                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (c_read) then

        write(outfile,'(A)') 'The integrals will be read from file'

      else 

        write(outfile,'(A)') 'The integrals will be calculated'

        select case (trim(calculation_type))
          case ("OBC", "Ring", "OBC2D")
            call molecule(n_atoms,number_of_functions,atoms,geometry,S,T,V,ERI)            ! Molecule 
          case ("Torus")
            call Torus_PBC(n_atoms,number_of_functions,atoms,AO,geometry)                  ! Torus with PBC 
          case ("Tori1D")
            call Tori1D(n_atoms,number_of_functions,atoms,AO,geometry,S,T,V,ERI)           ! Toroidal 1D Gaussian TRR
          case ("Tori2D")
            call Tori2D(n_atoms,number_of_functions,atoms,AO,geometry)                     ! Real Toroidal 2D Gaussian
          case ("Tori3D")
            call Tori3D(n_atoms,number_of_functions,atoms,AO,geometry,S,T,V,ERI)           ! Real Toroidal 3D Gaussian
          case default
            write(outfile,'(A)') 'Unknown calculation type: ',          &
            &                     trim(calculation_type)
            stop
        end select         

      end if


      ! --------------------------------------------------------------- !
      !            Read the one and the two electron integrals     
      ! --------------------------------------------------------------- !      

      if (c_Integral) then
        write(outfile,'(A)') 'All of The integrals calculated'
        stop 
      end if

      if (c_read) then
        call read_integrals_from_file(nBas,S,T,V,Hc,ERI,calculation_type)
      else 
        HC (:,:) = T(:,:) + V(:,:)
      end if


      ! --------------------------------------------------------------- !
      !                  write the integrals to outputfile              !
      ! --------------------------------------------------------------- !

      if (c_details) then 
        call details_integrals(nBas,S,T,V,ERI)
      end if 

      ! --------------------------------------------------------------- !

      if (c_trexio) then
        call trexio_conv_integrals(nBas,S,T,V,Hc,ERI)
      end if

      ! --------------------------------------------------------------- !

      !------------------------------------------------------!
      !                                  (-1/2)         t    !
      !  orthogonalization and get  X = S       =  U s U     !
      !                                                      !
      !------------------------------------------------------!

      !CALL HEADER ('The Overlap Matrix',-1)

      !call matout(nBas,nBas,S)

      ! --------------------------------------------------------------- !


      ! --------------------------------------------------------------- ! 
      !        pure 2D toroidal case (not 2D in 3D normal case)         !
      ! --------------------------------------------------------------- !

      if (calculation_type == "Tori2D" ) then 
        call get_X_from_overlap_2D(nBAS,S,X)
      else 
        call get_X_from_overlap(nBAS,S,X)
      end if 
       
      ! ---------------------------------------------------------------- !
      !                                                                  !
      !                           HF code start                          !
      !                                                                  !
      ! ---------------------------------------------------------------- ! 
 
      !if (calculation_type == "Tori2D" .or. calculation_type == "Tori3D" ) E_nuc = 0.d0 
      
      if (.not. c_UHF) then 
         call cpu_time(start)
           call RHF(nBas,c_details,c_Huckel,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
         call cpu_time(end)
      else 
        call cpu_time(start)
          call UHF(nBas,c_details,n_alpha,n_beta,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
        call cpu_time(end)
      end if 

        time = end - start

        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',time,' seconds'

      write(outfile,*)
      
      if (c_details) then 
        call system("tar -czf " // trim(output_file_name) // ".tar.gz "  // trim(tmp_file_name) )
      end if 

      
      !-----------------------------------------------------------------!
      !-----------------------------------------------------------------!


      !if (calculation_type == "Tori2D" .or. calculation_type == "OBC2D" .or. calculation_type == "Tori3D" ) then

      !  Hc(:,:) = T(:,:)

      !-----------------------------------------------------------------!
      ! AO to MO transformation
      !-----------------------------------------------------------------!

      if (c_MO .or. c_MP2) then 

      allocate(ERI_MO(nBas,nBas,nBas,nBas))
      allocate(Hc_MO(nBas,nBas))
      call cpu_time(start)
        call AO_to_MO_HC (nBas,c,T,HC_MO)
        call AO_to_MO_ERI(nBas,c,ERI,ERI_MO)
      call cpu_time(end)

      call print_MO_file(nBas,HC_MO,ERI_MO)
    
      time = end - start
      

      end if 

      !-----------------------------------------------------------------!
      ! FCI Energy calculation
      !-----------------------------------------------------------------!
      if (c_MP2) then 
        call cpu_time(start)
          call MP2(nBas,nO,e,ERI_MO,E_nuc,EHF)
        call cpu_time(end)

        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',time,' seconds'
        write(outfile,*)

        time = end - start
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 calculation = ',time,' seconds'
        write(outfile,*)

      end if

      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !


    
      !-----------------------------------------------------------------!
      ! FCI Energy calculation
      !-----------------------------------------------------------------!

      !call FCI(Hc_MO,ERI_MO,nBAS,E_nuc)
        
      !end if

      ! --------------------------------------------------------------- !
      ! --------------------------------------------------------------- !





      close(outfile)

      call system("rm torus_parameters.inp")
      if ( .not. c_details) then 
        call system("rm -r " // trim(tmp_file_name))
      end if 


      ! --------------------------------------------------------------- !
      !            close the TREXIO file and exit                       !
      ! --------------------------------------------------------------- !

      if (c_trexio) then
        call trexio_conv_close()
      end if

      ! --------------------------------------------------------------- !

end program CI