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
      type(atom)      ,allocatable    :: atoms_tor(:)
      type(ERI_function),allocatable  :: AO (:)

      double precision,allocatable    ::       S(:,:)
      double precision,allocatable    ::     S_T(:,:)
      double precision,allocatable    ::       T(:,:)
      double precision,allocatable    ::       V(:,:)
      double precision,allocatable    ::      Hc(:,:)
      double precision,allocatable    ::       X(:,:)
      double precision,allocatable    :: ERI(:,:,:,:)
      double precision,allocatable    ::         e(:)
      double precision,allocatable    ::       c(:,:)


      double precision                :: E_nuc , EHF
      double precision                :: start_HF,end_HF,t_HF

      character(len=5)                :: calculation_type 


      !-----------------------------------------------------------------!
      !                        END variables                            !
      !-----------------------------------------------------------------!    

      call read_geometry(n_atoms,charge_tmp,geometry_tmp)

      allocate(geometry(n_atoms,3))
      allocate(charge(n_atoms))

      do i = 1 , n_atoms
        charge    (i) =     charge_tmp(i)
        geometry(i,1) = geometry_tmp(i,1)
        geometry(i,2) = geometry_tmp(i,2)
        geometry(i,3) = geometry_tmp(i,3)
      end do 

      allocate(atoms    (n_atoms))
      allocate(atoms_tor(n_atoms))

      call basis(n_atoms,charge,atoms)

      call basis_tor(n_atoms,charge,atoms_tor)

      open (outfile,file="results.out")

      CALL HEADER ('The Geometry',-1)

      do i = 1 , n_atoms
        write(outfile,"(I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      number_of_functions = 0 
      do i = 1 , n_atoms
        number_of_functions = number_of_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      allocate(AO(number_of_functions))

      call classification_orbital(n_atoms,number_of_functions,geometry,atoms,AO)
      call print_orbital_table(AO,number_of_functions)


!     -------------------------------------------------------------------     !
!                         Parameters for the calculation 
!     -------------------------------------------------------------------     !

      calculation_type = "OBC"

!     -------------------------------------------------------------------     !
!                            Plot the gussians                                !
!     -------------------------------------------------------------------     !


!      call plot(n_atoms,geometry)


!     -------------------------------------------------------------------     !
!                         calculate the integrals 
!     -------------------------------------------------------------------     !

      call header_under("Calculate the integerals",-1)

      call cpu_time(start_HF)

!      call overlap_matrix_tor(n_atoms,geometry,atoms_tor)

      !-----------------------------------------------------------------!
      !                    calculate molecule                           !
      !-----------------------------------------------------------------!

      if (calculation_type == "OBC") then 
        
        call header("calculation with OBC",-1)

        call overlap_matrix(n_atoms,geometry,atoms)
        call kinetic_matrix(n_atoms,geometry,atoms)
        call nuclear_attraction_matrix(n_atoms,geometry,atoms)  
        call ERI_integral(n_atoms,geometry,atoms)

      end if 

      !-----------------------------------------------------------------!
      !                        calculate SAO                            !
      !-----------------------------------------------------------------!

      if (calculation_type == "SAO") then 
        
        call overlap_matrix(n_atoms,geometry,atoms)
        call SAO_S_matrix(number_of_functions)
        call kinetic_matrix(n_atoms,geometry,atoms)
        call SAO_k_matrix(number_of_functions)
        call nuclear_attraction_matrix(n_atoms,geometry,atoms)  
        call SAO_NA_matrix(number_of_functions)
!        call ERI_integral(n_atoms,geometry,atoms)

      end if 

      !-----------------------------------------------------------------!
      !           calculate the with Torus translational symmetry       !
      !-----------------------------------------------------------------!


      if (calculation_type == "Torus") then 

        call header("Torus with PBC",-1)
        call Torus_def()

        call overlap_matrix_torus(n_atoms,number_of_functions,atoms,AO)
        call kinetic_matrix_torus(n_atoms,number_of_functions,atoms,AO)
        call nuclear_attraction_matrix_torus(n_atoms,number_of_functions,geometry,atoms,AO)
        call ERI_integral_torus(n_atoms,geometry,atoms)

        write(outfile,*) "Translation symmetry applied to integrals"
        write(outfile,*) ""

      end if 

      !-----------------------------------------------------------------!

      call cpu_time(end_HF)

      t_HF = end_HF - start_HF

      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',t_HF,' seconds'
      write(outfile,*)
 
!     -------------------------------------------------------------------     !
!                        Nuclear repulsion energy  
!     -------------------------------------------------------------------     !

      call NRE(n_atoms,geometry,atoms,E_nuc)

!     -------------------------------------------------------------------     !
!                Read the one and the two electron integrals     
!     -------------------------------------------------------------------     !

      nBAS = 0 
      do i = 1 , n_atoms
        nBAS = nBAS + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 

      n_electron = 0 
      do i = 1 , n_atoms
        n_electron = n_electron + atoms(i)%charge
      end do 

      nO = n_electron/2 

      !               ---------------------------------                 !
      !                       Allocate the memory                       !
      !               ---------------------------------                 !

      allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas),ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas))
      allocate(S_T(nBas,nBas))

!     -------------------------------------------------------------------     !

      call read_integrals(nBas,S,T,V,Hc,ERI)

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

      call cpu_time(start_HF)
        call RHF(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
      call cpu_time(end_HF)

      t_HF = end_HF - start_HF
      write(outfile,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
      write(outfile,*)

      ! ---------------------------------------------------------------- !
      ! ---------------------------------------------------------------- !

      close(outfile)

end program CI 