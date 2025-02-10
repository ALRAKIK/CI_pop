subroutine nuclear_attraction_matrix(number_of_atoms,geometry,atoms)

      use atom_basis

      implicit none 

      integer                      :: i , j , k , l 
      integer                      :: func_index(1:number_of_atoms)
      integer                      :: number_of_atoms
      type(atom)                   :: atoms(number_of_atoms)
      type(atom)                   :: atom1 , atom2 

      double precision             :: geometry(number_of_atoms,3)

      double precision,allocatable :: NA(:,:)
      double precision             :: r1(3) , r2(3)

      double precision,parameter   :: pi = dacos(-1.d0)


      double precision             :: SS 
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










          

        end do 
      end do 



      open(1,file="NA_matrix")
        do i = 1 , size(NA,1)
          write(1,'(1000(f16.12,2x))')  (NA(i,j),j=1,size(NA,1))
        end do 
      close(1)


end subroutine