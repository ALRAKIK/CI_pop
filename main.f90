program CI 

      use atom_basis

      implicit none 

      integer                      :: i , j 
      integer                      :: n_atoms , nBAS
      integer                      :: nO      , n_electron 

      double precision             :: geometry_tmp(100,3)
      integer                      :: charge_tmp(100)


      double precision,allocatable :: geometry(:,:)
      integer         ,allocatable :: charge(:)
      type(atom)      ,allocatable :: atoms(:)

      double precision,allocatable :: S(:,:)
      double precision,allocatable :: T(:,:)
      double precision,allocatable :: V(:,:)
      double precision,allocatable :: Hc(:,:)
      double precision,allocatable :: X(:,:)
      double precision,allocatable :: ERI(:,:,:,:)
      double precision,allocatable :: e(:)
      double precision,allocatable :: c(:,:)


      double precision             :: E_nuc , EHF
      double precision             :: start_HF,end_HF,t_HF

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

      call system("clear")

      CALL HEADER ('The Geometry',-1)

      do i = 1 , n_atoms
        write(*,"(I2,3f16.8)") charge(i), (geometry(i,j),j=1,3)
      end do 

      

      
!     -------------------------------------------------------------------     !
!                         calculate the integrals 
!     -------------------------------------------------------------------     !

      call cpu_time(start_HF)

      call overlap_matrix(n_atoms,geometry,atoms)

      call kinetic_matrix(n_atoms,geometry,atoms)

      call nuclear_attraction_matrix(n_atoms,geometry,atoms)

      call ERI_integral(n_atoms,geometry,atoms)

      call cpu_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for integrals = ',t_HF,' seconds'
      write(*,*)

      

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

      allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas),ERI(nBas,nBas,nBas,nBas),e(nBas),c(nBas,nBas))

      call read_integrals(nBas,S,T,V,Hc,ERI)

      !------------------------------------------------------!
      !                                  (-1/2)         t    !
      !  orthogonalization and get  X = S       =  U s U     !
      !                                                      !
      !------------------------------------------------------!

      CALL HEADER ('The Overlap Matrix',-1)

      call matout(nBas,nBas,S)

      call get_X_from_overlap(nBAS,S,X)

      nO = n_electron/2 
      
      ! ---------------------------------------------------------------- !
      !                                                                  !
      !                           HF code start                          !
      !                                                                  !
      ! ---------------------------------------------------------------- ! 

      call cpu_time(start_HF)
        call RHF(nBas,nO,S,T,V,Hc,ERI,X,E_nuc,EHF,e,c)
      call cpu_time(end_HF)

      t_HF = end_HF - start_HF
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for HF = ',t_HF,' seconds'
      write(*,*)

      ! ---------------------------------------------------------------- !
      ! ---------------------------------------------------------------- !

end program CI 