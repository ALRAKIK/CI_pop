subroutine SAO_S_matrix(nbas)

      implicit none 

      integer ,intent(in) :: nbas 

      double precision             :: SAO_S(nbas,nbas) , S(nbas,nbas)

      double precision             :: OV , Kronecker_delta
      double precision,parameter   :: pi = 3.14159265358979323846D00

      integer             :: mu , nu , k , k_prime , i , j 

      open(1,file="./tmp/OV.dat")

      S(:,:) = 0d0

      do 
        read(1,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(1)

      SAO_S(:,:) = 0.d0

      write(*,*) "" 

      do k = 1 , nbas 
        do k_prime  = 1 , nbas 
          do mu = 1 , nbas 
            SAO_S(k,k_prime) = SAO_S(k,k_prime) + Kronecker_delta(k,k_prime) * dcos(((2*pi)/nbas)*k_prime*(mu-1)) * S(1,mu)
          end do 
        end do 
      end do 

      open(1,file="./tmp/OV.dat")
      do i = 1 , size(SAO_S,1)
        do j = i , size(SAO_S,1)
          write(1,'(I5,I5,f16.8)') i , j , SAO_S(i,j)
        end do 
      end do 
      close(1)

end subroutine 

subroutine SAO_K_matrix(nbas)

      implicit none 

      integer ,intent(in)                     :: nbas 

      double precision             :: SAO_S(nbas,nbas) , S(nbas,nbas)

      double precision             :: OV , Kronecker_delta
      double precision,parameter   :: pi = 3.14159265358979323846D00

      integer                      :: mu , nu , k , k_prime , i , j 

      open(1,file="./tmp/KI.dat")

      S(:,:) = 0d0

      do 
        read(1,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(1)

      SAO_S(:,:) = 0.d0

      do k = 1 , nbas 
        do k_prime  = 1 , nbas 
          do mu = 1 , nbas 
            SAO_S(k,k_prime) = SAO_S(k,k_prime) + Kronecker_delta(k,k_prime) * dcos(((2*pi)/nbas)*k_prime*(mu-1)) * S(1,mu)
          end do 
        end do 
      end do 

      open(1,file="./tmp/KI.dat")
      do i = 1 , size(SAO_S,1)
        do j = i , size(SAO_S,1)
          write(1,'(I5,I5,f16.8)') i , j , SAO_S(i,j)
        end do 
      end do 
      close(1)

end subroutine 


subroutine SAO_NA_matrix(nbas)

      implicit none 

      integer ,intent(in) :: nbas 
      double precision             :: SAO_S(nbas,nbas) , S(nbas,nbas)
      double precision             :: OV , Kronecker_delta
      double precision,parameter   :: pi = 3.14159265358979323846D00

      integer             :: mu , nu , k , k_prime , i , j 

      open(1,file="./tmp/NA.dat")

      S(:,:) = 0d0

      do 
        read(1,*,end=8) mu,nu,Ov
        S(mu,nu) = Ov
        S(nu,mu) = Ov
      enddo
      8 close(1)


      SAO_S(:,:) = 0.d0

      do k = 1 , nbas 
        do k_prime  = 1 , nbas 
          do mu = 1 , nbas 
            SAO_S(k,k_prime) = SAO_S(k,k_prime) + Kronecker_delta(k,k_prime) * dcos(((2*pi)/nbas)*k_prime*(mu-1)) * S(1,mu)
          end do 
        end do 
      end do 

      open(1,file="./tmp/NA.dat")
      do i = 1 , size(SAO_S,1)
        do j = i , size(SAO_S,1)
          write(1,'(I5,I5,f16.8)') i , j , SAO_S(i,j)
        end do 
      end do 
      close(1)

end subroutine 

function Kronecker_delta(i,j) result(delta)

      ! Kronecker Delta

      implicit none

      ! Input variables

      integer,intent(in)            :: i,j

      ! Output variables

      double precision              :: delta

      if(i == j) then
        delta = 1d0
      else
        delta = 0d0
      endif
  
end function Kronecker_delta