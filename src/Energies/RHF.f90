subroutine RHF(nBas,nO,S,T,V,Hc,ERI,X,ENuc,EHF,e,c)

      ! Perform a restricted Hartree-Fock calculation

      use files 
      implicit none

      ! Input variables

      integer,intent(in)            :: nBas
  
      integer,intent(in)            :: nO
      double precision,intent(in)   :: S(nBas,nBas)
      double precision,intent(in)   :: T(nBas,nBas)
      double precision,intent(in)   :: V(nBas,nBas)
      double precision,intent(in)   :: Hc(nBas,nBas) 
      double precision,intent(in)   :: X(nBas,nBas) 
      double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
      double precision,intent(in)   :: ENuc
  
      ! Local variables
  
      integer,parameter             :: maxSCF = 100
      double precision,parameter    :: thresh = 1d-5
      integer                       :: nSCF
      double precision              :: Conv
      double precision              :: Gap
      double precision              :: ET,EV,EJ , EHF_old
      double precision              :: EK
      double precision,allocatable  :: cp(:,:)
      double precision,allocatable  :: P(:,:)
      double precision,allocatable  :: J(:,:)
      double precision,allocatable  :: K(:,:)
      double precision,allocatable  :: F(:,:),Fp(:,:)
      double precision,allocatable  :: error(:,:)
      double precision,external     :: trace_matrix

      integer                       :: max_diis = 5
      integer                       :: n_diis
      integer                       :: i  , o
      integer                       :: mu , nu 
      double precision              :: sum 
      double precision              :: sum_o(nbas)
      double precision              :: rcond
      double precision,allocatable  :: err_diis(:,:)
      double precision,allocatable  :: F_diis(:,:)
  
      ! Output variables
  
      double precision,intent(out)  :: EHF 
      double precision,intent(out)  :: e(nBas)
      double precision,intent(out)  :: c(nBas,nBas)
    
      write(outfile,*)
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)'|                          Restricted Hartree-Fock calculation                           |'
      write(outfile,*)'******************************************************************************************'
      write(outfile,*)
  
      ! Memory allocation
  
      allocate(cp(nBas,nBas),P(nBas,nBas),                         &
             J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
             error(nBas,nBas))

      allocate(err_diis(nBas*nBas,max_diis))
      allocate(F_diis(nBas*nBas,max_diis))
  
      ! Guess coefficients and eigenvalues
  
      F(:,:) = Hc(:,:)

      open(1,file="./tmp/FOCK_matrix.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(F,1))
      do i = 1 , size(F,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (F(i,o),o=1,size(F,1))
      end do 
      close(1)

      open(1,file="./tmp/X.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
      do i = 1 , size(X,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
      end do 
      close(1)

      cp(:,:) = matmul(transpose(X(:,:)), matmul(F(:,:), X(:,:)))

      open(1,file="./tmp/C_prime.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(cp,1))
      do i = 1 , size(cp,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (cp(i,o),o=1,size(cp,1))
      end do 
      close(1)

      call diagonalize_matrix(nBAS, cp, e)
      
      c(:,:) = matmul(X(:,:), cp(:,:))

      open(1,file="./tmp/C.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(c,1))
      do i = 1 , size(c,1)
        write(1,'(i3,6x,1000(f16.10,2x))') i ,  (c(i,o),o=1,size(c,1))
      end do 
      close(1)
      
      P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))

      open(1,file="./tmp/density_matrix.dat")
        write(1,'(15x,1000(i3,15x))') (i,i=1,size(P,1))
        do i = 1 , size(P,1)
          write(1,'(i3,6x,1000(f16.10,2x))') i ,  (P(i,o),o=1,size(P,1))
        end do 
      close(1)



      ! --------------------------------------------------------------- !
      !              check that P_{mu nu} S_{mu nu} = N 
      ! --------------------------------------------------------------- !

      sum = 0.d0 
      do mu = 1 , nbas 
        do nu = 1 , nbas 
          sum = sum + P(mu,nu)*S(mu,nu)
        end do 
      end do 

      call header_under("P_{mu nu} S_{mu nu} = N",-1)

      write(outfile,"(a,f16.8)") "P_{mu nu} S_{mu nu} =  " , sum 

      ! --------------------------------------------------------------- !

      ! --------------------------------------------------------------- !
      !              apply the density on the functions 
      ! --------------------------------------------------------------- !

      sum_o(:) = 0.d0 

      do o = 1 , nbas
        do mu = 1 , nbas 
          do nu = 1 , nbas 
            sum_o(o) = sum_o(o) + P(mu,nu)*S(mu,o)*S(nu,o)
          end do 
        end do 
      end do 

      call header_under("P_{mu nu} S_{mu s} S_{nu s}",-1)

      write(outfile,"(a)") "P_{mu nu} S_{mu s} S_{nu s} :  "
      do i = 1 , nbas 
        write(outfile,"(40x,i3,1x,f16.8)") i , sum_o(i)
      end do 

      ! --------------------------------------------------------------- !

      ! Initialization
  
      nSCF = 0
      Conv = 1d0

      n_diis        = 0
      F_diis(:,:)   = 0d0
      err_diis(:,:) = 0d0
      rcond         = 0d0

      EHF_old       = 0d0 
      EHF           = 0d0 
  
      !------------------------------------------------------------------------
      ! Main SCF loop
      !------------------------------------------------------------------------
  
      write(outfile,*)
      write(outfile,*) repeat('-', 110)
      write(outfile,*) "|",repeat(' ', 47),"RHF calculation",repeat(' ', 46),"|"
      write(outfile,*) repeat('-', 110)

      write(outfile,'(1x,a1,a2,1x,a1,a,a1,a,a1,a,a1,a,a1,a,a1,a,a1)') &
      "|","#","|","      HF energy  ","|", "    Conv   ","|","   HL Gap  ","|", "    T contribution  ",&
      "|","   V contribution   ","|","   Two contribution ","|"
      write(outfile,*) repeat('-', 110)
      
      do while(nSCF < maxSCF)
  
      !   Increment 
  
      nSCF = nSCF + 1
  
      !   Compute the Hartree potential J
      call hartree_potential(nBas,P,ERI,J)
      ! ****************** !
  
      !   Compute the exchange potential K
      call exchange_potential(nBas,P,ERI,K)
      ! ****************** !

      !   Build Fock operator
    
      F(:,:) = Hc(:,:) + J(:,:) + K(:,:)
  
      !   Compute the error vector and extract the convergence criterion
      error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
      if(nSCF > 1) Conv = maxval(abs(error)) 
      if (Conv < thresh) exit 
      ! ****************** !
  
      !------------------------------------------------------------------------
      !   Compute HF energy
      !------------------------------------------------------------------------
  
      !   Compute the kinetic energy
      ET = trace_matrix(nBas,matmul(P,T))
      ! ****************** !

      !   Compute the potential energy
      EV = trace_matrix(nBas,matmul(P,V))
      ! ****************** !

      !   Compute the Hartree energy
      EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))
      ! ****************** !

      !   Compute the exchange energy
      EK = 0.5d0*trace_matrix(nBas,matmul(P,K))
      ! ****************** !

      !   Total HF energy
  
      EHF = ET + EV + EJ + EK

      if (abs(EHF - EHF_old) < thresh ) exit  !.and. nSCF > maxSCF/4 ) exit

      if (nSCF > 2) then 
        EHF_old = EHF 
      end if 

      !   Compute HOMO-LUMO gap
          if(nBas > nO) then
            Gap = e(nO+1) - e(nO)
          else
            Gap = 0d0
          endif
      ! ****************** !

      ! DIIS extrapolation ! 

      if(max_diis > 1) then
        n_diis = min(n_diis+1,max_diis)
        call DIIS_extrapolation(rcond,nBas*nBas,nBas*nBas,n_diis,err_diis,F_diis,error,F)
      end if

      !   Transform for the Fock matrix F in the orthogonal basis 
      
      Fp(:,:) = matmul(transpose(X(:,:)), matmul(F(:,:),X(:,:)))
            
      !   Diagonalize F' to get MO coefficients (eigenvectors in the orthogonal basis) c' and MO energies (eigenvalues) e

      cp(:,:) = Fp(:,:)

      call diagonalize_matrix(nBas,cp,e)

      !   Back-transform the MO coefficients c in the original non-orthogonal basis

      c = matmul(X,cp)

      !   Compute the density matrix P 

      P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

      !   Dump results
   
      write(outfile,"(1x,a1,i2,1x,a1,f16.8,1x,a1,f10.6,1x,a1,f10.6,1x,a1,f16.8,4x,a1,f16.8,4x,a1,f16.8,4x,a1)")      & 
      "|",nSCF,"|",EHF+ENuc,"|",Conv,"|",Gap,"|",ET,"|",EV,"|",Ej+EK,"|"
      enddo

      write(outfile,*) repeat('-', 110)

      !------------------------------------------------------------------------
      ! End of SCF loop
      !------------------------------------------------------------------------
  
      if(nSCF == maxSCF) then
  
        write(outfile,*)
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)'                 Convergence failed                 '
        write(outfile,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(outfile,*)

        call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)
  
        stop
  
      endif
  
      ! Compute final HF energy

        call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)
  
end subroutine RHF
  