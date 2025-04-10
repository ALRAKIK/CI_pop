subroutine classification(number_of_atoms,number_of_functions,geometry,atoms,ERI)
    
      use atom_basis
      use classification_ERI
      use filter_module

      implicit none 

      integer                       :: i , j , k
      integer                       :: number_of_functions
      integer                       :: number_of_atoms 
      double precision              :: geometry(number_of_atoms,3)
      type(atom)                    :: atoms(number_of_atoms)
      type(ERI_function)            :: ERI  (number_of_functions)
      double precision, allocatable :: f_array_coef(:), f_array_expo(:)
      
      j = 0 
      
      do k = 1 , number_of_atoms

        do i = 1 , atoms(k)%num_s_function
              
          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "s"

          call filter(atoms(k)%coefficient_s(:,i),atoms(k)%exponent_s,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

        end do 

        do i = 1 , atoms(k)%num_p_function

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "px"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "py"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          j = j + 1

          ERI(j)%x = geometry(k,1) 
          ERI(j)%y = geometry(k,2) 
          ERI(j)%z = geometry(k,3)

          ERI(j)%orbital     = "pz"

          call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)

          ERI(j)%exponent    = f_array_expo
          ERI(j)%coefficient = f_array_coef

          if (allocated(f_array_coef)) deallocate(f_array_coef)
          if (allocated(f_array_expo)) deallocate(f_array_expo)

          end do 

      end do 


end subroutine


subroutine classification_orbital(number_of_atoms,number_of_functions,geometry,atoms,AO)
    
        use atom_basis
        use classification_ERI
        use filter_module
  
        implicit none 
  

        integer           ,intent(in)   :: number_of_functions
        integer           ,intent(in)   :: number_of_atoms 
        double precision  ,intent(in)   :: geometry(number_of_atoms,3)
        type(atom)        ,intent(in)   :: atoms(number_of_atoms)
        type(ERI_function),intent(out)  :: AO  (number_of_functions)

        integer                         :: i , j , k
        double precision, allocatable   :: f_array_coef(:), f_array_expo(:)
        
        j = 0 
        
        do k = 1 , number_of_atoms
  
          do i = 1 , atoms(k)%num_s_function
                
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "s"
  
            call filter(atoms(k)%coefficient_s(:,i),atoms(k)%exponent_s,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
          end do 
  
          do i = 1 , atoms(k)%num_p_function
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "px"
  
            call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "py"
  
            call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "pz"
  
            call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            end do 
  
        end do 
  
  
end subroutine classification_orbital

subroutine classification_orbital_tor(number_of_atoms,number_of_functions,geometry,atoms,norm_helper,AO)
    
        use atom_basis
        use classification_ERI
        use filter_module
  
        implicit none 
  

        integer           ,intent(in)   :: number_of_functions
        integer           ,intent(in)   :: number_of_atoms 
        double precision  ,intent(in)   :: geometry(number_of_atoms,3)
        type(atom)        ,intent(in)   :: atoms(number_of_atoms)
        type(atom)        ,intent(in)   :: norm_helper(number_of_atoms)
        type(ERI_function),intent(out)  :: AO  (number_of_functions)

        integer                         :: i , j , k
        double precision, allocatable   :: f_array_coef(:), f_array_expo(:)
        
        j = 0 
        
        do k = 1 , number_of_atoms
  
          do i = 1 , atoms(k)%num_s_function
                
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "s"
  
            call filter(atoms(k)%coefficient_s(:,i),atoms(k)%exponent_s,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
          end do 
  
          do i = 1 , atoms(k)%num_p_function
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "px"
  
            call filter(atoms(k)%coefficient_p(:,i),atoms(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "py"
  
            call filter(norm_helper(k)%coefficient_p(:,i),norm_helper(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            j = j + 1
  
            AO(j)%x = geometry(k,1) 
            AO(j)%y = geometry(k,2) 
            AO(j)%z = geometry(k,3)
  
            AO(j)%orbital     = "pz"
  
            call filter(norm_helper(k)%coefficient_p(:,i),norm_helper(k)%exponent_p,f_array_coef,f_array_expo)
  
            AO(j)%exponent    = f_array_expo
            AO(j)%coefficient = f_array_coef
  
            if (allocated(f_array_coef)) deallocate(f_array_coef)
            if (allocated(f_array_expo)) deallocate(f_array_expo)
  
            end do 
  
        end do 
  
  
end subroutine classification_orbital_tor 

