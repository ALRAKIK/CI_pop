subroutine AO_to_MO_ERI(nBas,c,ERI_AO,ERI_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  integer                       :: i,j,k,l
  double precision,allocatable  :: scr(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas,nBas,nBas,nBas)

! Memory allocation

  allocate(scr(nBas,nBas,nBas,nBas))
  
!--------------------------------------
! AO to MO transformation starts here !
!--------------------------------------

      ! Transform 4th index

      scr(:,:,:,:) = 0d0

      do mu=1,nBas
         do nu=1,nBas
           do la=1,nBas
             do si=1,nBas
            
               do s=1,nBas
              
                 scr(mu,nu,la,s) = scr(mu,nu,la,s) + c(si,s)*ERI_AO(mu,nu,la,si)
              
               enddo    
             
             enddo
          enddo
        enddo
      enddo

      ! Transform 3rd index

      ERI_MO(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
          do la=1,nBas
          
            do r=1,nBas
              do s=1,nBas
              
                ERI_MO(mu,nu,r,s) = ERI_MO(mu,nu,r,s) + c(la,r)*scr(mu,nu,la,s)
              
              enddo
            enddo    
          
          enddo
        enddo
      enddo

      ! Transform 2nd index

      scr(:,:,:,:) = 0d0

      do mu=1,nBas
        do nu=1,nBas
        
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
              
                scr(mu,q,r,s) = scr(mu,q,r,s) + c(nu,q)*ERI_MO(mu,nu,r,s)
              
              enddo
            enddo
          enddo    
        
        enddo
      enddo

      ! Transform 1st (and last) index

      ERI_MO(:,:,:,:) = 0d0

      do mu=1,nBas
      
        do p=1,nBas
          do q=1,nBas
            do r=1,nBas
              do s=1,nBas
              
                ERI_MO(p,q,r,s) = ERI_MO(p,q,r,s) + c(mu,p)*scr(mu,q,r,s)
              
              enddo
            enddo
          enddo
        enddo
      
      enddo

      open(1,file="./tmp/ERI_MO.dat")
        do i = 1, nbas
          do j = 1 , nbas
            do k = 1 , nbas
              do l = 1 , nbas
                if (abs(ERI_MO(i,j,k,l)) > 1e-8 ) write(1,"(I5,I5,I5,I5,f16.10)") i , j , k , l ,  ERI_MO(i,j,k,l)
              end do 
            end do 
          end do 
        end do 
      close(1)

end subroutine AO_to_MO_ERI

subroutine AO_to_MO_HC(nBas,c,HC_AO,HC_MO)

! Expression of bi-electronic integrals in the MO basis set

  implicit none
!  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: HC_AO(nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si
  integer                       :: p,q,r,s
  integer                       :: i,j,k,l

! Output variables

  double precision,intent(out)  :: HC_MO(nBas,nBas)

! Memory allocation

  
      !--------------------------------------
      ! AO to MO transformation starts here !
      !--------------------------------------

      HC_MO(:,:) = 0d0

      do i=1,2
         do j=1,2
            HC_MO(i,j) = 0.d0 
            do k=1,2
               do l=1,2
                  HC_MO(i,j) = HC_MO(i,j)+c(k,i)*c(l,j)*HC_AO(k,l)
               enddo
            enddo
         enddo
      enddo



      open(1,file="./tmp/HC_MO.dat")
      write(1,'(15x,1000(i3,15x))') (i,i=1,size(HC_MO,1))
      do i = 1 , size(HC_MO,1)
        write(1,'(i3,6x,1000(f16.12,2x))') i ,  (HC_MO(i,j),j=1,size(HC_MO,1))
      end do 
      close(1)

end subroutine AO_to_MO_HC
