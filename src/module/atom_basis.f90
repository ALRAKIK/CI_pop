module atom_basis
      
      implicit none

      type :: atom
        integer                         :: charge 
        integer                         :: num_s_function
        integer                         :: num_p_function
        integer                         :: num_exponent_s
        integer                         :: num_exponent_p 
        double precision  , allocatable :: exponent_s(:)
        double precision  , allocatable :: exponent_p(:)
        double precision  , allocatable :: coefficient_s(:,:)
        double precision  , allocatable :: coefficient_p(:,:)
      end type atom

end module 