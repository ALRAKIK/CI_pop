subroutine DIIS_extrapolation(rcond,n_err,n_e,n_diis,error,e,error_in,e_inout)

      ! Perform DIIS extrapolation
  
      implicit none
    
      ! Input variables
  
      integer,intent(in)             :: n_err,n_e
      double precision,intent(in)    :: error_in(n_err)
      double precision,intent(inout) :: error(n_err,n_diis)
      double precision,intent(inout) :: e(n_e,n_diis)
  
      ! Local variables
  
      double precision,allocatable  :: A(:,:)
      double precision,allocatable  :: b(:)
      double precision,allocatable  :: w(:)
  
      ! Output variables
  
      double precision,intent(out)  :: rcond
      integer,intent(inout)         :: n_diis
      double precision,intent(inout):: e_inout(n_e)
  
      ! Memory allocaiton
  
      allocate(A(n_diis+1,n_diis+1),b(n_diis+1),w(n_diis+1))
  
      ! Update DIIS "history"
  
      call prepend(n_err,n_diis,error,error_in)
      call prepend(n_e,n_diis,e,e_inout)
  
      !  Build A matrix
    
      A(1:n_diis,1:n_diis) = matmul(transpose(error),error)
    
      A(1:n_diis,n_diis+1) = -1d0
      A(n_diis+1,1:n_diis) = -1d0
      A(n_diis+1,n_diis+1) = +0d0
  
      ! Build x matrix
  
      b(1:n_diis) = +0d0
      b(n_diis+1) = -1d0
  
      ! Solve linear system
  
      call linear_solve(n_diis+1,A,b,w,rcond)
  
      ! Extrapolate
  
      e_inout(:) = matmul(w(1:n_diis),transpose(e(:,1:n_diis)))
  
end subroutine 

subroutine linear_solve(N,A,b,x,rcond)

      ! Solve the linear system A.x = b where A is a NxN matrix
      ! and x and x are vectors of size N
    
      implicit none
    
      integer,intent(in)             :: N
      double precision,intent(out)   :: A(N,N),b(N),rcond
      double precision,intent(out)   :: x(N)
    
      integer                        :: info,lwork
      double precision               :: ferr,berr
      integer,allocatable            :: ipiv(:),iwork(:)
      double precision,allocatable   :: AF(:,:),work(:)
    
      ! Find optimal size for temporary arrays
    
      allocate(work(1))
      allocate(AF(N,N),ipiv(N),iwork(N))
    
      lwork = -1
      call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)
      lwork = int(work(1))
    
      deallocate(work)
    
      allocate(work(lwork))
    
      call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)
    
      ! if (info /= 0) then

      !   print *,  info
      !   stop 'error in linear_solve (dsysvx)!!'

      ! end if
    
end subroutine

subroutine prepend(N,M,A,b)

      ! Prepend the vector b of size N into the matrix A of size NxM
  
      implicit none
  
      ! Input variables
  
      integer,intent(in)            :: N,M
      double precision,intent(in)   :: b(N)
  
      ! Local viaruabkes
  
      integer                       :: i,j
  
      ! Output variables
  
      double precision,intent(out)  :: A(N,M)
  
  
      ! print*,'b in append'
      ! call matout(N,1,b)
  
      do i=1,N
        do j=M-1,1,-1
          A(i,j+1) = A(i,j)
        end do
        A(i,1) = b(i)
      end do
  
end subroutine