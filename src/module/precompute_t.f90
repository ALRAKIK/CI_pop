module precomputed_bessel

      implicit none

      integer , parameter           :: Nmax_global = 5000 
      double precision, allocatable :: I_C_table_x(:, :)  ! (0:Nmax, 1:32 or 1:64)
      double precision, allocatable :: I_C_table_y(:, :)  ! (0:Nmax, 1:32 or 1:64)
      double precision, allocatable :: I_C_table_z(:, :)  ! (0:Nmax, 1:32 or 1:64)
      double precision, allocatable :: t_nodes(:)
      double precision, allocatable :: C_nodes(:)


      contains

      subroutine initialize_bessel_table_32_Lx()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_ax, inv_ax2
      integer                     :: i 
      integer                     :: npoint = 32
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_ax  = 1.d0/ax 
      inv_ax2 = inv_ax * inv_ax

      ! ============================================================
      ! Define the 32 fixed t values (from your printout)
      ! ============================================================
      t_nodes = [ &
          1.3699498041200081d-003, &
          7.2463831197837796d-003, &
          1.7934856587675952d-002, &
          3.3641909391797641d-002, &
          5.4673661692177175d-002, &
          8.1450751017493023d-002, &
          0.11452665180428827d0,   &
          0.15461298665154638d0,   &
          0.20261380375999602d0,   &
          0.25967132739808052d0,   &
          0.32722817813835642d0,   &
          0.40711168506416140d0,   &
          0.50164964969641823d0,   &
          0.61383070797638206d0,   &
          0.74753080276989492d0,   &
          0.90783685683233750d0,   &
          1.1015194993175779d0,    &
          1.3377375170288737d0,    &
          1.6291136741866559d0,    &
          1.9934231003752658d0,    &
          2.4563284147503617d0,    &
          3.0559715415986783d0,    &
          3.8510220208755785d0,    &
          4.9354978853491103d0,    &
          6.4677620014786665d0,    &
          8.7315920289792022d0,    &
          12.277357636459755d0,    &
          18.290342535134830d0,    &
          29.724828883934091d0,    &
          55.757345764735931d0,    &
          137.99987986694228d0,    &
          729.95375231456273d0     &
      ]
    
      ! ============================================================
      ! Compute the 32 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_ax2
      end do
    
      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_x(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 32 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 32 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_x(0:Nmax_global, i) = temp_array(0:Nmax_global)

      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_32_Lx

subroutine initialize_bessel_table_32_Ly()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_ay, inv_ay2
      integer                     :: i 
      integer                     :: npoint = 32
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_ay  = 1.d0/ay 
      inv_ay2 = inv_ay * inv_ay

      ! ============================================================
      ! Define the 32 fixed t values (from your printout)
      ! ============================================================
      t_nodes = [ &
          1.3699498041200081d-003, &
          7.2463831197837796d-003, &
          1.7934856587675952d-002, &
          3.3641909391797641d-002, &
          5.4673661692177175d-002, &
          8.1450751017493023d-002, &
          0.11452665180428827d0,   &
          0.15461298665154638d0,   &
          0.20261380375999602d0,   &
          0.25967132739808052d0,   &
          0.32722817813835642d0,   &
          0.40711168506416140d0,   &
          0.50164964969641823d0,   &
          0.61383070797638206d0,   &
          0.74753080276989492d0,   &
          0.90783685683233750d0,   &
          1.1015194993175779d0,    &
          1.3377375170288737d0,    &
          1.6291136741866559d0,    &
          1.9934231003752658d0,    &
          2.4563284147503617d0,    &
          3.0559715415986783d0,    &
          3.8510220208755785d0,    &
          4.9354978853491103d0,    &
          6.4677620014786665d0,    &
          8.7315920289792022d0,    &
          12.277357636459755d0,    &
          18.290342535134830d0,    &
          29.724828883934091d0,    &
          55.757345764735931d0,    &
          137.99987986694228d0,    &
          729.95375231456273d0     &
      ]
    
      ! ============================================================
      ! Compute the 32 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_ay2
      end do
    
      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_y(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 32 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 32 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_y(0:Nmax_global, i) = temp_array(0:Nmax_global)

      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_32_Ly

subroutine initialize_bessel_table_32_Lz()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_az, inv_az2
      integer                     :: i 
      integer                     :: npoint = 32
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_az  = 1.d0/az
      inv_az2 = inv_az * inv_az

      ! ============================================================
      ! Define the 32 fixed t values (from your printout)
      ! ============================================================
      t_nodes = [ &
          1.3699498041200081d-003, &
          7.2463831197837796d-003, &
          1.7934856587675952d-002, &
          3.3641909391797641d-002, &
          5.4673661692177175d-002, &
          8.1450751017493023d-002, &
          0.11452665180428827d0,   &
          0.15461298665154638d0,   &
          0.20261380375999602d0,   &
          0.25967132739808052d0,   &
          0.32722817813835642d0,   &
          0.40711168506416140d0,   &
          0.50164964969641823d0,   &
          0.61383070797638206d0,   &
          0.74753080276989492d0,   &
          0.90783685683233750d0,   &
          1.1015194993175779d0,    &
          1.3377375170288737d0,    &
          1.6291136741866559d0,    &
          1.9934231003752658d0,    &
          2.4563284147503617d0,    &
          3.0559715415986783d0,    &
          3.8510220208755785d0,    &
          4.9354978853491103d0,    &
          6.4677620014786665d0,    &
          8.7315920289792022d0,    &
          12.277357636459755d0,    &
          18.290342535134830d0,    &
          29.724828883934091d0,    &
          55.757345764735931d0,    &
          137.99987986694228d0,    &
          729.95375231456273d0     &
      ]
    
      ! ============================================================
      ! Compute the 32 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_az2
      end do
    
      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_z(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 32 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 32 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_z(0:Nmax_global, i) = temp_array(0:Nmax_global)

      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_32_Lz


subroutine initialize_bessel_table_64_Lx()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_ax, inv_ax2
      integer                     :: i 
      integer                     :: npoint = 64
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_ax  = 1.d0/ax 
      inv_ax2 = inv_ax * inv_ax

      ! ============================================================
      ! Define the 64 fixed t values (from your printout)
      ! ============================================================
      
      t_nodes = [ &
        0.00034758605080146475d0, &
        0.001833306865139011d0, &
        0.004513607056424145d0, &
        0.008401868305325152d0, &
        0.013516867918589854d0, &
        0.019883242288482272d0, &
        0.027532015873410566d0, &
        0.036500798092868027d0, &
        0.046834231406787760d0, &
        0.058584486267713690d0, &
        0.071811705480770233d0, &
        0.086584909928179909d0, &
        0.102982544212746070d0, &
        0.121093588425561030d0, &
        0.141018540955522880d0, &
        0.16287058918283662d0,  &
        0.18677716513660425d0,  &
        0.21288166838664668d0,  &
        0.24134523925221177d0,  &
        0.27234908283331816d0,  &
        0.30609733378807802d0,  &
        0.34281976554502269d0,  &
        0.38277585567787348d0,  &
        0.42625868077786250d0,  &
        0.47360026026799029d0,  &
        0.52517739833830801d0,  &
        0.58141877132733433d0,  &
        0.64281349707107416d0,  &
        0.70992138611004618d0,  &
        0.78338502252242681d0,  &
        0.86394484635210744d0,  &
        0.95245709764026898d0,  &
        1.0499160565630929d0,   &
        1.1574812955044149d0,   &
        1.2765115125383597d0,   &
        1.4086066704926512d0,   &
        1.5556611747519555d0,   &
        1.7199307096966907d0,   &
        1.9041185000802747d0,   &
        2.1114853261147752d0,   &
        2.3459932784832436d0,   &
        2.6124949762807241d0,   &
        2.9169846680519638d0,   &
        3.2669346956557788d0,   &
        3.6717582802069337d0,   &
        4.1434419966120615d0,   &
        4.6974453346717873d0,   &
        5.3539735399058257d0,   &
        6.1398439400094009d0,   &
        7.0912661074503607d0,   &
        8.2580755348143189d0,   &
        9.7103835183383502d0,   &
        11.549356589150186d0,   &
        13.925306373175895d0,   &
        17.069365350927502d0,   &
        21.351903724314528d0,   &
        27.396661230686686d0,   &
        36.321350554129388d0,   &
        50.293608330632686d0,   &
        73.981635836264417d0,   &
        119.02114668545737d0,   &
        221.55229453940080d0,   &
        545.46242040291190d0,   &
        2876.9854189896218d0 &
      ]
    
      ! ============================================================
      ! Compute the 64 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_ax2
      end do
    
      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_x(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 64 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 64 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_x(0:Nmax_global, i) = temp_array(0:Nmax_global)
        
      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_64_Lx

subroutine initialize_bessel_table_64_Ly()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_ay, inv_ay2
      integer                     :: i 
      integer                     :: npoint = 64
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_ay  = 1.d0/ay 
      inv_ay2 = inv_ay * inv_ay

      ! ============================================================
      ! Define the 64 fixed t values (from your printout)
      ! ============================================================
      
      t_nodes = [ &
        0.00034758605080146475d0, &
        0.001833306865139011d0, &
        0.004513607056424145d0, &
        0.008401868305325152d0, &
        0.013516867918589854d0, &
        0.019883242288482272d0, &
        0.027532015873410566d0, &
        0.036500798092868027d0, &
        0.046834231406787760d0, &
        0.058584486267713690d0, &
        0.071811705480770233d0, &
        0.086584909928179909d0, &
        0.102982544212746070d0, &
        0.121093588425561030d0, &
        0.141018540955522880d0, &
        0.16287058918283662d0,  &
        0.18677716513660425d0,  &
        0.21288166838664668d0,  &
        0.24134523925221177d0,  &
        0.27234908283331816d0,  &
        0.30609733378807802d0,  &
        0.34281976554502269d0,  &
        0.38277585567787348d0,  &
        0.42625868077786250d0,  &
        0.47360026026799029d0,  &
        0.52517739833830801d0,  &
        0.58141877132733433d0,  &
        0.64281349707107416d0,  &
        0.70992138611004618d0,  &
        0.78338502252242681d0,  &
        0.86394484635210744d0,  &
        0.95245709764026898d0,  &
        1.0499160565630929d0,   &
        1.1574812955044149d0,   &
        1.2765115125383597d0,   &
        1.4086066704926512d0,   &
        1.5556611747519555d0,   &
        1.7199307096966907d0,   &
        1.9041185000802747d0,   &
        2.1114853261147752d0,   &
        2.3459932784832436d0,   &
        2.6124949762807241d0,   &
        2.9169846680519638d0,   &
        3.2669346956557788d0,   &
        3.6717582802069337d0,   &
        4.1434419966120615d0,   &
        4.6974453346717873d0,   &
        5.3539735399058257d0,   &
        6.1398439400094009d0,   &
        7.0912661074503607d0,   &
        8.2580755348143189d0,   &
        9.7103835183383502d0,   &
        11.549356589150186d0,   &
        13.925306373175895d0,   &
        17.069365350927502d0,   &
        21.351903724314528d0,   &
        27.396661230686686d0,   &
        36.321350554129388d0,   &
        50.293608330632686d0,   &
        73.981635836264417d0,   &
        119.02114668545737d0,   &
        221.55229453940080d0,   &
        545.46242040291190d0,   &
        2876.9854189896218d0 &
      ]
    
      ! ============================================================
      ! Compute the 64 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_ay2
      end do
    
      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_y(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 64 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 64 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_y(0:Nmax_global, i) = temp_array(0:Nmax_global)
        
      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_64_Ly

subroutine initialize_bessel_table_64_Lz()

      use bessel_functions
      use torus_init

      implicit none
    
      double precision, parameter :: pi = dacos(-1.d0)
      double precision            :: inv_az, inv_az2
      integer                     :: i 
      integer                     :: npoint = 64
      double precision            :: temp_array(0:Nmax_global)  ! Temporary array for one C value

      ! =============================================================== !

      allocate(t_nodes(npoint))
      allocate(C_nodes(npoint))
    
      ! ============================================================
      ! Set your constant Lx
      ! ============================================================

      inv_az  = 1.d0/az 
      inv_az2 = inv_az * inv_az

      ! ============================================================
      ! Define the 64 fixed t values (from your printout)
      ! ============================================================
      
      t_nodes = [ &
        0.00034758605080146475d0, &
        0.001833306865139011d0, &
        0.004513607056424145d0, &
        0.008401868305325152d0, &
        0.013516867918589854d0, &
        0.019883242288482272d0, &
        0.027532015873410566d0, &
        0.036500798092868027d0, &
        0.046834231406787760d0, &
        0.058584486267713690d0, &
        0.071811705480770233d0, &
        0.086584909928179909d0, &
        0.102982544212746070d0, &
        0.121093588425561030d0, &
        0.141018540955522880d0, &
        0.16287058918283662d0,  &
        0.18677716513660425d0,  &
        0.21288166838664668d0,  &
        0.24134523925221177d0,  &
        0.27234908283331816d0,  &
        0.30609733378807802d0,  &
        0.34281976554502269d0,  &
        0.38277585567787348d0,  &
        0.42625868077786250d0,  &
        0.47360026026799029d0,  &
        0.52517739833830801d0,  &
        0.58141877132733433d0,  &
        0.64281349707107416d0,  &
        0.70992138611004618d0,  &
        0.78338502252242681d0,  &
        0.86394484635210744d0,  &
        0.95245709764026898d0,  &
        1.0499160565630929d0,   &
        1.1574812955044149d0,   &
        1.2765115125383597d0,   &
        1.4086066704926512d0,   &
        1.5556611747519555d0,   &
        1.7199307096966907d0,   &
        1.9041185000802747d0,   &
        2.1114853261147752d0,   &
        2.3459932784832436d0,   &
        2.6124949762807241d0,   &
        2.9169846680519638d0,   &
        3.2669346956557788d0,   &
        3.6717582802069337d0,   &
        4.1434419966120615d0,   &
        4.6974453346717873d0,   &
        5.3539735399058257d0,   &
        6.1398439400094009d0,   &
        7.0912661074503607d0,   &
        8.2580755348143189d0,   &
        9.7103835183383502d0,   &
        11.549356589150186d0,   &
        13.925306373175895d0,   &
        17.069365350927502d0,   &
        21.351903724314528d0,   &
        27.396661230686686d0,   &
        36.321350554129388d0,   &
        50.293608330632686d0,   &
        73.981635836264417d0,   &
        119.02114668545737d0,   &
        221.55229453940080d0,   &
        545.46242040291190d0,   &
        2876.9854189896218d0 &
      ]
    
      ! ============================================================
      ! Compute the 64 C values
      ! ============================================================
      do i = 1, npoint
          C_nodes(i) = 2.0d0 * (t_nodes(i) ** 2) * inv_az2
      end do

      ! ============================================================
      ! Allocate the big table
      ! ============================================================

      allocate(I_C_table_z(0:Nmax_global, npoint))
    
      ! ============================================================
      ! Precompute I_C for all 64 C values (THIS IS THE KEY!)
      ! ============================================================

      print *, "Precomputing Bessel functions for 64 quadrature points..."
    
      do i = 1, npoint
        
        print *, "Computing I_C for point ", i, " (C =", C_nodes(i), ")"
        
        call bessel_I_scaled_n(Nmax_global, C_nodes(i), temp_array)
        
        I_C_table_z(0:Nmax_global, i) = temp_array(0:Nmax_global)
        
      end do
    
      print *, "Precomputation COMPLETE!"

      deallocate(t_nodes)
      deallocate(C_nodes)
    
end subroutine initialize_bessel_table_64_Lz


      subroutine bessel_I_scaled_n(Nmax, C, I_array)

        use bessel_functions
        
        implicit none 

        integer, intent(in)           :: Nmax
        double precision, intent(in)  :: C
        double precision, intent(out) :: I_array(0:Nmax)

        integer                       :: i 

        do i = 0 , Nmax 
          I_array(i) = iv_scaled(i, C)
        end do

      end subroutine bessel_I_scaled_n


end module precomputed_bessel