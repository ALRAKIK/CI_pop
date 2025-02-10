subroutine nuclear_attraction_integral_ss(n_atoms,geometry,r1,r2,atom1,atom2,atoms,index1,index2,S_ss_normal)

  use atom_basis
  implicit none 

  double precision,intent(in)  :: r1(3) , r2(3)
  type(atom),intent(in)        :: atom1 , atom2 
  type(atom)                   :: atoms(n_atoms)
  integer,intent(in)           :: n_atoms
  double precision,intent(in)  :: geometry(n_atoms,3)
  integer                      :: index1 , index2 

  integer                      :: i , j  , k 
  integer                      :: charge_atom
  double precision,parameter   :: pi = dacos(-1.d0)
  double precision             :: alpha , beta
  double precision             :: c1    , c2 
  double precision             :: p,mu
  double precision             :: Two_PIP , P_R
  double precision             :: x1 , x2 , y1 , y2 , z1 , z2 
  double precision             :: xp , yp , zp 
  double precision             :: X , Y , Z
  double precision             :: D_normal 
  double precision             :: xPC , yPC , zPC
  double precision             :: boys_function , R2PC
  double precision,intent(out) :: S_ss_normal

  x1 = r1(1) ; x2 = r2(1) 
  y1 = r1(2) ; y2 = r2(2)
  z1 = r1(3) ; z2 = r2(3)

  X            = (x1 - x2)
  Y            = (y1 - y2)
  Z            = (z1 - z2)

  D_normal     = (X*X+Y*Y+Z*Z)

  !-----------------------------------------------------------------!

  S_ss_normal = 0.d0
  do i = 1 , atom1%num_exponent_s
    alpha = atom1%exponent_s(i)
    c1    = atom1%coefficient_s(i,index1)
    do j = 1 , atom2%num_exponent_s
      beta = atom2%exponent_s(j)
      c2   = atom2%coefficient_s(j,index2)
        if (c1*c2 == 0.d0) cycle  
          p       = alpha + beta 
          P_R     = 1.d0 / p
          mu      = alpha*beta*p_R
          Two_PIP = 2.0d0*pi*p_R
          
          xp = (alpha*x1+beta*x2)*p_R
          yp = (alpha*y1+beta*y2)*p_R
          zp = (alpha*z1+beta*z2)*p_R

          do k = 1 , n_atoms

            xPC = xp - geometry(k,1) 
            yPC = yp - geometry(k,2)
            zPC = zp - geometry(k,3)

            R2PC = xPC*xPC + yPC*yPC + zPC*zPC

            charge_atom = (-1)*atoms(k)%charge

            S_ss_normal =  S_ss_normal +  c1 * c2 * charge_atom * Two_PIP * exp(-mu*D_normal) * boys_function(0,p*R2PC)

          end do 

    end do 
  end do



!-----------------------------------------------------------------!

end subroutine