subroutine guess_RHF(nBas,nO,HC,X,ENuc,T,V,P)

      use files

      implicit none 

      ! --------------------------------------------------------------- !
      integer         , intent(in)  :: nBas , nO 
      double precision, intent(in)  :: HC(nBas,nBas), X(nBas,nBas)
      double precision, intent(in)  ::  T(nBas,nBas), V(nBas,nBas)
      double precision, intent(in)  :: ENuc
      double precision,external     :: trace_matrix
      ! --------------------------------------------------------------- !
      double precision,allocatable  :: cp(:,:), c(:,:) , e(:) , ct(:,:)
      integer                       :: i      , o 
      double precision              :: ET     , EV
      ! --------------------------------------------------------------- !
      double precision, intent(out) :: P(nbas,nbas)
      ! --------------------------------------------------------------- !


      allocate(cp(nbas,nbas),c(nbas,nbas),ct(nbas,nbas), e(nbas))

      cp(:,:) = matmul(transpose(X(:,:)), matmul(HC(:,:), X(:,:)))

      call diagonalize_matrix(nbas, cp, e)

      c(:,:)  = matmul(X(:,:), cp(:,:))

      ct(:,:) = transpose(c)

      P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))


      ! --------------------------------------------------------------- !
      !                          Print                                  !
      ! --------------------------------------------------------------- !


      write(HFfile,'(a)') "!--------------------------------------------!"
      write(HFfile,'(a)') "                 Iter  =   0"
      write(HFfile,'(a)') "!--------------------------------------------!"
      write(HFfile,'(a)') ""
      
      call header_HF("Guess F Matrix, F = HC ", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(HC,1))
      do i = 1 , size(HC,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (HC(i,o),o=1,size(HC,1))
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess X Matrix, X = S^(-1/2)", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(X,1))
      do i = 1 , size(X,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (X(i,o),o=1,size(X,1))
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess CP Matrix, CP = X^t F X", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(cp,1))
      do i = 1 , size(cp,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (cp(i,o),o=1,size(cp,1))
      end do 
      write(HFfile,'(a)') ""

      call header_HF("Guess C Matrix, C = X CP", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(c,1))
      do i = 1 , size(c,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (c(i,o),o=1,size(c,1))
      end do 
      write(HFfile,'(a)') ""

      call header_HF(" C^t Matrix, C^t = transpose(C)", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(c,1))
      do i = 1 , size(c,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  ((ct(i,o)),o=1,size(transpose(c),1))
      end do 
      write(HFfile,'(a)') ""


      call header_HF("Guess Density Matrix, P = 2 C C^t", -1)
      write(HFfile,'(15x,1000(i3,15x))') (i,i=1,size(P,1))
      do i = 1 , size(P,1)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  (P(i,o),o=1,size(P,1))
      end do 
      write(HFfile,'(a)') ""
      
      call header_HF("Guess Orbital Energies", -1)
      do i = 1 , size(e)
        write(HFfile,'(i3,6x,1000(f16.10,2x))') i ,  e(i)
      end do 
      write(HFfile,'(a)') ""



      write(HFfile,'(a)') ""
      ET = trace_matrix(nBas,matmul(P,T))
      EV = trace_matrix(nBas,matmul(P,V))
      write(HFfile,'(a,f16.10)')   " The Kinetic   Energy    = ", ET
      write(HFfile,'(a,f16.10)')   " The Potential Energy    = ", EV
      write(HFfile,'(a,f16.10,a)') " The Nuclear   Energy    = ", ENuc , "  +"
      write(HFfile,'(a,f16.10)')   "-----------------------------------"
      write(HFfile,'(a,f16.10)') " The Energy of the guess  = ", ET+EV+ENuc
      write(HFfile,'(a)') ""

      write(HFfile,'(2a)')  repeat('*_',36) , "*"
      write(HFfile,'(a)')   repeat('_',73)


      ! --------------------------------------------------------------- !
      deallocate(cp, e)
      ! --------------------------------------------------------------- !


end subroutine guess_RHF 