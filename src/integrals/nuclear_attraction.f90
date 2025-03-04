subroutine nuclear_attraction_matrix(number_of_atoms,geometry,atoms)

      use atom_basis

      implicit none 

      integer                      :: i , j , k , l 
      integer                      :: func_index(1:number_of_atoms)
      integer                      :: number_of_atoms
      integer                      :: idx_s  , idx_p 
      integer                      :: idx_p1 , idx_p2
      type(atom)                   :: atoms(number_of_atoms)
      type(atom)                   :: atom1 , atom2 

      double precision             :: geometry(number_of_atoms,3)

      double precision,allocatable :: NA(:,:)
      double precision             :: r1(3) , r2(3)

      double precision,parameter   :: pi = dacos(-1.d0)


      double precision             :: SS , SP(3) , PS(3) , PP(3,3)
      integer                      :: total_functions



      total_functions = 0 

      do i = 1 , number_of_atoms
        total_functions = total_functions + atoms(i)%num_s_function + 3 * atoms(i)%num_p_function
      end do 
    
      allocate(NA(total_functions,total_functions))
    
      NA(:,:) = 0.d0 
    
      func_index(1) = 1
      do i = 2, number_of_atoms
        func_index(i) = func_index(i-1) + atoms(i-1)%num_s_function + 3 * atoms(i-1)%num_p_function
      end do
    
      do i = 1 , number_of_atoms
        do j = i , number_of_atoms

          r1(1) = geometry(i,1) ; r2(1) = geometry(j,1)
          r1(2) = geometry(i,2) ; r2(2) = geometry(j,2)
          r1(3) = geometry(i,3) ; r2(3) = geometry(j,3)
        
          atom1 = atoms(i)
          atom2 = atoms(j)
        
          do k = 1 , atom1%num_s_function
            do l = 1 , atom2%num_s_function
              call nuclear_attraction_integral_ss(number_of_atoms,geometry,r1,r2,atom1,atom2,atoms,k,l,SS)
              NA(func_index(i) + k - 1, func_index(j) + l - 1) = SS    ! SS 
              NA(func_index(j) + l - 1, func_index(i) + k - 1) = SS    ! SS 
            end do 
          end do 

          do k = 1 , atom1%num_s_function
            do l = 1 , atom2%num_p_function
              call nuclear_attraction_integral_sp(number_of_atoms,geometry,r1,r2,atom1,atom2,atoms,k,l,SP)
              idx_s = func_index(i) + k - 1 
              idx_p = func_index(j) + atom2%num_s_function - 1 + (l-1) * 3 
              NA(  idx_s, idx_p+1) = SP(1)                             !  SPx
              NA(idx_p+1,   idx_s) = SP(1)
              NA(  idx_s, idx_p+2) = SP(2)                             !  SPy
              NA(idx_p+2,   idx_s) = SP(2)  
              NA(  idx_s, idx_p+3) = SP(3)                             !  SPz
              NA(idx_p+3,   idx_s) = SP(3)  
            end do 
          end do 

          do k = 1 , atom1%num_p_function
            do l = 1 , atom2%num_s_function
              call nuclear_attraction_integral_ps(number_of_atoms,geometry,r1,r2,atom1,atom2,atoms,k,l,PS)
              idx_p = func_index(i) + atom1%num_s_function - 1 + (k-1) * 3  
              idx_s = func_index(j) + l - 1                                 
              NA(idx_p+1,   idx_s) = PS(1)                               ! PxS
              NA(  idx_s, idx_p+1) = PS(1)                               
              NA(idx_p+2,   idx_s) = PS(2)                               ! PyS
              NA(  idx_s, idx_p+2) = PS(2)                               
              NA(idx_p+3,   idx_s) = PS(3)                               ! PzS
              NA(  idx_s, idx_p+3) = PS(3)                               
            end do 
          end do

          do k = 1 , atom1%num_p_function
            do l = 1 , atom2%num_p_function
              call nuclear_attraction_integral_pp(number_of_atoms,geometry,r1,r2,atom1,atom2,atoms,k,l,PP)
              idx_p1 = func_index(i) + atom1%num_s_function - 1 + (k-1) * 3  ! p-function index in atom1
              idx_p2 = func_index(j) + atom2%num_s_function - 1 + (l-1) * 3  ! p-function index in atom2
              
              NA(idx_p1+1, idx_p2+1) = PP(1,1)                            ! PxPx
              NA(idx_p2+1, idx_p1+1) = PP(1,1)                            ! PxPx
              NA(idx_p1+1, idx_p2+2) = PP(1,2)                            ! PxPy
              NA(idx_p2+2, idx_p1+1) = PP(1,2)                            ! PxPy
              NA(idx_p1+1, idx_p2+3) = PP(1,3)                            ! PxPz
              NA(idx_p2+3, idx_p1+1) = PP(1,3)                            ! PxPz

              NA(idx_p1+2, idx_p2+1) = PP(2,1)                            ! PyPx
              NA(idx_p2+1, idx_p1+2) = PP(2,1)                            ! PyPx
              NA(idx_p1+2, idx_p2+2) = PP(2,2)                            ! PyPy
              NA(idx_p2+2, idx_p1+2) = PP(2,2)                            ! PyPy
              NA(idx_p1+2, idx_p2+3) = PP(2,3)                            ! PyPz
              NA(idx_p2+3, idx_p1+2) = PP(2,3)                            ! PyPz

              NA(idx_p1+3, idx_p2+1) = PP(3,1)                            ! PzPx
              NA(idx_p2+1, idx_p1+3) = PP(3,1)                            ! PzPx
              NA(idx_p1+3, idx_p2+2) = PP(3,2)                            ! PzPy
              NA(idx_p2+2, idx_p1+3) = PP(3,2)                            ! PzPy
              NA(idx_p1+3, idx_p2+3) = PP(3,3)                            ! PzPz
              NA(idx_p2+3, idx_p1+3) = PP(3,3)                            ! PzPz

            end do 
          end do

        end do 
      end do 



      open(1,file="./tmp/NA.dat")
        do i = 1 , size(NA,1)
          do j = i , size(NA,1)
            write(1,'(I5,I5,f16.8)') i , j ,  NA(i,j)
          end do 
        end do 
      close(1)

      open(1,file="./tmp/NA_matrix.dat")
      do i = 1 , size(NA,1)
        write(1,'(1000(f16.12,2x))')  (NA(i,j),j=1,size(NA,1))
      end do 
      close(1)


end subroutine