subroutine read_integrals(nBas,S,T,V,Hc,G)

      ! Read one- and two-electron integrals from files

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas

      ! Local variables

      integer                       :: mu,nu,la,si
      double precision              :: Ov,Kin,Nuc,ERI

      double precision,intent(out)  :: S(nBas,nBas)
      double precision,intent(out)  :: T(nBas,nBas)
      double precision,intent(out)  :: V(nBas,nBas)
      double precision,intent(out)  :: Hc(nBas,nBas)
      double precision,intent(out)  :: G(nBas,nBas,nBas,nBas)

      open( 8,file='./tmp/OV.dat' )
      open( 9,file='./tmp/KI.dat' )
      open(10,file='./tmp/NA.dat' )
      open(11,file='./tmp/ERI.dat')

      S(:,:) = 0d0
      do 
        read(8,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(8)

      ! Read kinetic integrals

      T(:,:) = 0d0
      do 
        read(9,*,end=9) mu,nu,Kin
        T(mu,nu) = Kin
        T(nu,mu) = Kin
      enddo
      9 close(9)

      ! Read nuclear integrals

      V(:,:) = 0d0
      do 
        read(10,*,end=10) mu,nu,Nuc
        V(mu,nu) = Nuc
        V(nu,mu) = Nuc
      enddo
      10 close(unit=10)

      ! Define core Hamiltonian

      Hc(:,:) = T(:,:) + V(:,:)

      ! Read ERI integrals

      G(:,:,:,:) = 0.d0

      do 
      read(11,*,end=11) mu,nu,la,si,ERI
        G(mu, nu, la, si) = ERI
      enddo
      11 close(unit=11)


end subroutine


subroutine read_overlap_T(nBas,S)

      ! Read one- and two-electron integrals from files

      implicit none

      ! Input variables

      integer,intent(in)            :: nBas

      ! Local variables

      integer                       :: mu,nu
      double precision              :: Ov

      double precision,intent(out)  :: S(nBas,nBas)
      
      open( 8,file='./tmp/OV_tor.dat' )

      S(:,:) = 0d0
      do 
        read(8,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(8)

end subroutine

!------------------------------------------------------------------------!
subroutine get_X_from_overlap(N,over,X)

      use files 
      implicit none 

      ! input !

      integer,intent(in)            :: N
      double precision, intent(in)  :: over(N,N)

      ! local !

      integer                       :: i
      double precision              :: Evec(N,N)
      double precision              :: Eval(N)
      double precision,parameter    :: thresh = 1d-6

      ! output !

      double precision, intent(out) :: X(N,N)

      Evec(:,:) = over(:,:)

      call diagonalize_matrix(N,Evec,Eval)

      ! //////////////////////////////////////////////////////////////// !
      call header_under("The Overlap Eigenvalues",-1)

      write(outfile, "(5f16.8)") (Eval(i), i=1,n)

      write(outfile,*) ""
      write(outfile,*) ""
      write(outfile,'(a,g16.8)') " The smallest eigenvalue : " , MINVAL(Eval)
      write(outfile,*) ""

      if (Minval(Eval) < 0.d0) then 
        call header("Error",-1)
          write(outfile,'(a80)')'*******************************************************************************************'
          write(outfile,'(a80)')           "* The smallest eigenvalue of the overlap is negative, exiting the program      *"
          write(outfile,'(a80)')'*******************************************************************************************'
          stop            
      end if
      ! //////////////////////////////////////////////////////////////// !

      do i = 1 , N
          if (Eval(i) > thresh) then
            Eval(i) = dsqrt(1.d0 / Eval(i))
          else
            write(outfile,"(a,I3,a)") 'Eigenvalue ',i,' is too small for Lowdin orthogonalization'
          end if
      end do
    
      call ADAt(N,Evec,Eval,X)
  
end subroutine get_X_from_overlap

!------------------------------------------------------------------------
subroutine AtDA(N,A,D,B)

      ! Perform B = At.D.A where A is a NxN matrix and D is a diagonal matrix given
      ! as a vector of length N

      implicit none

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(in)   :: A(N,N),D(N)

      ! Local viaruabkes

      integer                       :: i,j,k

      ! Output variables
  
      double precision,intent(out)  :: B(N,N)
  
      B = 0d0
  
      do i=1,N
        do j=1,N
          do k=1,N
            B(i,k) = B(i,k) + A(j,i)*D(j)*A(j,k)
          enddo
        enddo
      enddo
  
end subroutine AtDA
  
!------------------------------------------------------------------------
subroutine ADAt(N,A,D,B)
  
      ! Perform B = A.D.At where A is a NxN matrix and D is a diagonal matrix given 
      ! as a vector of length N

      implicit none

      ! Input variables
  
      integer,intent(in)            :: N
      double precision,intent(in)   :: A(N,N),D(N)
    
      ! Local viaruabkes
  
      integer                       :: i,j,k
  
      ! Output variables
  
      double precision,intent(out)  :: B(N,N)
  
      B = 0d0

      do i=1,N
        do j=1,N
          do k=1,N
            B(i,k) = B(i,k) + A(i,j)*D(j)*A(k,j)
          enddo
        enddo
      enddo
  
end subroutine ADAt
!------------------------------------------------------------------------
subroutine diagonalize_matrix(N,A,e)

      ! Diagonalize a square matrix

      implicit none

      ! Input variables

      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)

      ! Local variables

      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)


      ! Memory allocation

      allocate(work(3*N))
      lwork = size(work)

      call dsyev('V','U',N,A,N,e,work,lwork,info)

      if(info /= 0) then
        write(20,'(a)') 'Problem in diagonalize_matrix (dsyev)!!'
        stop
      endif

      do i = 1 , N
        if (abs(e(i)) < 1e-10) e(i) = 0
      end do

end subroutine diagonalize_matrix
!------------------------------------------------------------------------
subroutine matout(m,n,A)

      ! Print the MxN array A
  
      use files
      implicit none

      integer,parameter             :: ncol = 5
      double precision,parameter    :: small = 1d-10
      integer,intent(in)            :: m,n
      double precision,intent(in)   :: A(m,n)
      double precision              :: B(ncol)
      integer                       :: ilower,iupper,num,i,j
    
      do ilower=1,n,ncol
        iupper = min(ilower + ncol - 1,n)
        num = iupper - ilower + 1
        write(outfile,'(3X,10(9X,I6))') (j,j=ilower,iupper)
        do i=1,m
          do j=ilower,iupper
            B(j-ilower+1) = A(i,j)
          enddo
          do j=1,num
            if(abs(B(j)) < small) B(j) = 0d0
          enddo
          write(outfile,'(I7,10F15.8)') i,(B(j),j=1,num)
        enddo
      enddo
  
end subroutine matout

!------------------------------------------------------------------------
function trace_matrix(n,A) result(Tr)

      ! Calculate the trace of the square matrix A
  
      implicit none
  
      ! Input variables
  
      integer,intent(in)            :: n
      double precision,intent(in)   :: A(n,n)
  
      ! Local variables
  
      integer                       :: i
  
      ! Output variables

      double precision              :: Tr
  
      Tr = 0d0
      do i=1,n
        Tr = Tr + A(i,i)
      enddo
  
end function trace_matrix
!------------------------------------------------------------------------