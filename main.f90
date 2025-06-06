program CI 

      use files
      use torus_init
      use atom_basis
      use classification_ERI

      implicit none 

      !-----------------------------------------------------------------!
      !                         START variables                         !
      !-----------------------------------------------------------------!


      integer                         :: i , j
      integer                         :: n_atoms , nBAS
      integer                         :: nO      , n_electron 

      double precision                :: geometry_tmp(100,3)
      integer                         :: charge_tmp(100)
      integer                         :: number_of_functions

      double precision,allocatable    :: geometry(:,:)
      integer         ,allocatable    :: charge(:)
      type(atom)      ,allocatable    :: atoms(:)
      type(atom)      ,allocatable    :: norm_helper(:)
      type(ERI_function),allocatable  :: AO (:)

      double precision,allocatable    ::          S(:,:)
      double precision,allocatable    ::          T(:,:)
      double precision,allocatable    ::          V(:,:)
      double precision,allocatable    ::         Hc(:,:)
      double precision,allocatable    ::      Hc_MO(:,:)
      double precision,allocatable    ::          X(:,:)
      double precision,allocatable    ::    ERI(:,:,:,:)
      double precision,allocatable    :: ERI_MO(:,:,:,:)
      double precision,allocatable    ::            e(:)
      double precision,allocatable    ::          c(:,:)
      double precision,allocatable    :: S_SAO_diag(:,:)


      double precision                :: E_nuc , EHF
      double precision                :: start_HF,end_HF,t_HF

      character(len=10)               :: calculation_type 


      !-----------------------------------------------------------------!
      !                        END variables                            !
      !-----------------------------------------------------------------!    

      call build_super_molecule()

      call read_geometry(n_atoms,charge_tmp,geometry_tmp,calculation_type)

      if (calculation_type == "Torus" .or. calculation_type == "Tori" .or. calculation_type == "Tori2D" ) call Torus_def()

      allocate(geometry(n_atoms,3))
      allocate(charge(n_atoms))

      do i = 1 , n_atoms
        charge    (i) =     charge_tmp(i)
        geometry(i,1) = geometry_tmp(i,1)
        geometry(i,2) = geometry_tmp(i,2)
        geometry(i,3) = geometry_tmp(i,3)
      end do 

      allocate(atoms    (n_atoms))
      allocate(norm_helper(n_atoms))

      if (calculation_type == "Tori" .or. calculation_type == "Tori2D" ) then 
        call basis_tor(n_atoms,charge,atoms,norm_helper,calculation_type)
      else
        call basis(n_atoms,charge,atoms)
      end if 

      open (outfile,file="results.out")

      CALL HEADER ('The Geometry',-1)

      do i = 1 , n_atoms
        write(outfile,"(I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      number_of_functions = 0 
      do i = 1 , n_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      nBAS = 0 
      do i = 1 , n_atoms
        nBAS = nBAS + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(AO(number_of_functions))

      if (calculation_type == "Tori" .or. calculation_type == "Tori2D") then 
        call classification_orbital_tor(n_atoms,number_of_functions,geometry,atoms,norm_helper,AO)
      else
        call classification_orbital(n_atoms,number_of_functions,geometry,atoms,AO)
      end if

!      call classification_orbital(n_atoms,number_of_functions,geometry,atoms,AO)
      call print_orbital_table(AO,number_of_functions)

!      call check_openmp_enabled()

      !               ---------------------------------                 !
      !                       Allocate the memory                       !
      !               ---------------------------------                 !

      allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas),ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas))

!     -------------------------------------------------------------------     !
!                            Plot the gussians                                !
!     -------------------------------------------------------------------     !

!      call plot(n_atoms,geometry,calculation_type)

!     -------------------------------------------------------------------     !
!                         calculate the integrals 
!     -------------------------------------------------------------------     !

      !-----------------------------------------------------------------!
      !                    calculate molecule                           !
      !-----------------------------------------------------------------!

      if (calculation_type == "OBC" .or. calculation_type == "Ring" ) then 
        call header("calculation with OBC",-1)
        call header_under("Calculate the integerals",-1)
        call cpu_time(start_HF)
        call overlap_matrix(n_atoms,geometry,atoms)
        call kinetic_matrix(n_atoms,geometry,atoms)
        call nuclear_attraction_matrix(n_atoms,geometry,atoms)  
        call ERI_integral(n_atoms,geometry,atoms)
        call cpu_time(end_HF)
        t_HF = end_HF - start_HF
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',t_HF,' seconds'
        write(outfile,*)
      end if 

      !-----------------------------------------------------------------!
      !                        calculate SAO                            !
      !-----------------------------------------------------------------!

      if (calculation_type == "SAO") then 
        call header("Torus with SAO",-1)
        call header_under("Calculate the integerals",-1)
        call Torus_def()
        allocate(S_SAO_diag(number_of_functions,number_of_functions))
        call cpu_time(start_HF)
        call overlap_matrix_SAO(n_atoms,number_of_functions,atoms,AO,S_SAO_diag)
        call kinetic_matrix_SAO(n_atoms,number_of_functions,atoms,AO,S_SAO_diag)
        call nuclear_attraction_matrix_SAO(n_atoms,number_of_functions,geometry,atoms,AO,S_SAO_diag)
        call ERI_integral_SAO(n_atoms,number_of_functions,geometry,atoms,S_SAO_diag)
        call cpu_time(end_HF)
        t_HF = end_HF - start_HF
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',t_HF,' seconds'
        write(outfile,*)
        deallocate(S_SAO_diag)
      end if 

      !-----------------------------------------------------------------!
      !           calculate the with Torus translational symmetry       !
      !-----------------------------------------------------------------!


      if (calculation_type == "Torus") then 

        call header("Torus with PBC",-1)
        call header_under("Calculate the integerals",-1)
        call Torus_def()
        write(outfile,*) ""
        write(outfile,'(a,f12.8)') "The length of the box:  Lx = ", Lx 
        write(outfile,*) ""
        call overlap_matrix_torus(n_atoms,number_of_functions,atoms,AO)
        call kinetic_matrix_torus(n_atoms,number_of_functions,atoms,AO)
        call nuclear_attraction_matrix_torus(n_atoms,number_of_functions,geometry,atoms,AO)
        call ERI_integral_torus(n_atoms,geometry,atoms)
        write(outfile,*) ""
        call system("rm torus_parameters.inp")
      end if 

      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!
      !        calculate the with Toroidal translational symmetry       !
      !-----------------------------------------------------------------!


      if (calculation_type == "Tori") then 

        call header("Toroidal Gaussian",-1)
        call header_under("Calculate the integerals",-1)
        call Torus_def()
        write(outfile,*) ""
        write(outfile,'(a,f12.8)') "The length of the box:  Lx = ", Lx 
        write(outfile,*) ""
        call cpu_time(start_HF)
        call overlap_matrix_toroidal(n_atoms,number_of_functions,atoms,AO)
        call kinetic_matrix_toroidal(n_atoms,number_of_functions,atoms,AO)
        call nuclear_attraction_matrix_toroidal(n_atoms,number_of_functions,geometry,atoms,AO)
        call cpu_time(end_HF)
        t_HF = end_HF - start_HF
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for NA  integrals = ',t_HF,' seconds'
        call ERI_integral_toroidal(n_atoms,geometry,atoms)
!        call overlap_matrix_toroidal_num(n_atoms,number_of_functions,atoms,AO)
!        write(outfile,*) ""
        call system("rm torus_parameters.inp")
      end if 
      
      !-----------------------------------------------------------------!

      !-----------------------------------------------------------------!
      !        calculate the with Toroidal translational symmetry       !
      !-----------------------------------------------------------------!

      if (calculation_type == "Tori2D") then 

        call header("Toroidal  2D Real Gaussian",-1)
        call header_under("Calculate the integerals",-1)
        call Torus_def()
        write(outfile,*) ""
        write(outfile,'(a,f12.8,a,f12.8,a)') "The length of the box:  Lx , Ly = (", Lx ," , ", Ly , " )"  
        write(outfile,*) ""
        call cpu_time(start_HF)
        call overlap_matrix_toroidal_2D(n_atoms,number_of_functions,atoms,AO)
        call kinetic_matrix_toroidal_2D(n_atoms,number_of_functions,atoms,AO)
        !call nuclear_attraction_matrix_toroidal(n_atoms,number_of_functions,geometry,atoms,AO)
        call cpu_time(end_HF)
        t_HF = end_HF - start_HF
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for NA  integrals = ',t_HF,' seconds'
        call ERI_integral_toroidal_2D(n_atoms,geometry,atoms)
        call system("rm torus_parameters.inp")
      end if 
      
      !-----------------------------------------------------------------!

      
!     -------------------------------------------------------------------     !
!                        Nuclear repulsion energy  
!     -------------------------------------------------------------------     !

      call NRE(n_atoms,geometry,atoms,E_nuc)

!     -------------------------------------------------------------------     !
!                Read the one and the two electron integrals     
!     -------------------------------------------------------------------     !

      n_electron = 0 
      do i = 1 , n_atoms
        n_electron = n_electron + atoms(i)%charge
      end do 

      nO = n_electron/2 

      call read_integrals(nBas,S,T,V,Hc,ERI)

!      call split_matrix(nBas,S)
!      call split_matrix(nBas,T)
!      call split_matrix(nBas,V)
!      call split_matrix(nBas,HC)
!      call split_matrix_ERI(nBas,ERI)

!      call read_overlap_T(nBas,S_T)

!      call check_symmetric_matrix(nBas,S,T,V,HC)

      !------------------------------------------------------!
      !                                  (-1/2)         t    !
      !  orthogonalization and get  X = S       =  U s U     !
      !                                                      !
      !------------------------------------------------------!

      CALL HEADER ('The Overlap Matrix',-1)

      call matout(nBas,nBas,S)

      call get_X_from_overlap(nBAS,S,X)
       
      ! ---------------------------------------------------------------- !
      !                                                                  !
      !                           HF code start                          !
      !                                                                  !
      ! ---------------------------------------------------------------- ! 

      if (calculation_type == "SAO") then 
      call cpu_time(start_HF)
        call RHF_SAO(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
      call cpu_time(end_HF)
      else 
        if (calculation_type == "Tori2D") E_nuc = 0.d0 
        call cpu_time(start_HF)
        call RHF(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
      call cpu_time(end_HF)
      end if 
      t_HF = end_HF - start_HF
      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
      write(outfile,*)

      if (calculation_type == "Tori2D") then
        
        
      !-----------------------------------------------------------------!
      ! AO to MO transformation
      !-----------------------------------------------------------------!

        allocate(ERI_MO(nBas,nBas,nBas,nBas))
        allocate(Hc_MO(nBas,nBas))

        call cpu_time(start_HF)
        call AO_to_MO_HC (nBas,c,HC,HC_MO)
        call AO_to_MO_ERI(nBas,c,ERI,ERI_MO)
        call cpu_time(end_HF)
    
        t_HF = end_HF - start_HF
        write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_HF,' seconds'
        write(outfile,*)
    

      !-----------------------------------------------------------------!
      ! FCI Energy calculation
      !-----------------------------------------------------------------!


      call FCI(Hc_MO,ERI_MO,nBAS)
        
      end if

      ! ---------------------------------------------------------------- !
      ! ---------------------------------------------------------------- !

      close(outfile)


end program CI 