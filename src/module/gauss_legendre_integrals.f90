module gauss_legendre_quadrature
      implicit none
      
      ! Precomputed nodes and weights for different accuracies
      ! These are symmetric: x_i = -x_{n+1-i}
      
contains

      ! ============================================================
      ! 16-point Gauss-Legendre (good for low accuracy, fastest)
      ! ============================================================
      subroutine gauss_legendre_16(f, a, b, result)
            interface
                  function f(t) result(y)
                        double precision, intent(in) :: t
                        double precision :: y
                  end function
            end interface
            double precision, intent(in) :: a, b
            double precision, intent(out) :: result
            
            ! Precomputed nodes and weights for 16 points
            double precision, parameter :: x(8) = [ &
                  -0.095012509837637440185, &  ! +/- x
                  -0.281603550779258913230, &
                  -0.458016777657227386342, &
                  -0.617876244402643748447, &
                  -0.755404408355003033895, &
                  -0.865631202387831743880, &
                  -0.944575023073232576077, &
                  -0.989400934991649932596  &
            ]
            double precision, parameter :: w(8) = [ &
                  0.189450610455068496285, &
                  0.182603415044923588867, &
                  0.169156519395002538189, &
                  0.149595988816576732081, &
                  0.124628971255533872052, &
                  0.095158511682492784810, &
                  0.062253523938647892863, &
                  0.027152459411754094852  &
            ]
            
            integer :: i
            double precision :: mid, half, t
            
            mid = (a + b) / 2.0d0
            half = (b - a) / 2.0d0
            
            result = 0.0d0
            
            ! Sum over symmetric pairs
            do i = 1, 8
                  t = mid + half * x(i)
                  result = result + w(i) * (f(t) + f(2.0d0*mid - t))
            end do
            
            result = result * half
            
      end subroutine gauss_legendre_16
      
      
      ! ============================================================
      ! 32-point Gauss-Legendre (balanced accuracy/speed)
      ! ============================================================
      subroutine gauss_legendre_32(f, a, b, result)
            interface
                  function f(t,i_point) result(y)
                        double precision, intent(in) :: t
                        integer, intent(in)          :: i_point
                        double precision             :: y
                  end function
            end interface
            double precision, intent(in) :: a, b
            double precision, intent(out) :: result
            
            ! Precomputed nodes and weights for 32 points
            double precision :: x(16), w(16)
            integer          :: i
            double precision :: mid, half, t

            integer          :: i_left, i_right
            
            ! Nodes (positive values only)
            x = [ &
                  0.048307665687738316235, &
                  0.144471961582796493485, &
                  0.239287362252137074545, &
                  0.331868602282127649780, &
                  0.421351276130635345364, &
                  0.506899908932229390024, &
                  0.587715757240762329041, &
                  0.663044266930215200975, &
                  0.732182118740289680387, &
                  0.794483795967942406963, &
                  0.849367613732569970133, &
                  0.896321155766052123965, &
                  0.934906075937739689171, &
                  0.964762255587506430773, &
                  0.985611511545268335400, &
                  0.997263861849481563545  &
            ]
            
            w = [ &
                  0.096540088514727800567, &
                  0.095638720079274859419, &
                  0.093844399080804565639, &
                  0.091173878695763884713, &
                  0.087652093004403811143, &
                  0.083311924226946755222, &
                  0.078193895787070306472, &
                  0.072345794108848506225, &
                  0.065822222776361846838, &
                  0.058684093478535547145, &
                  0.050998059262376176196, &
                  0.042835898022226680657, &
                  0.034273862913021433103, &
                  0.025392065309262059456, &
                  0.016274394730905670605, &
                  0.007018610009470096600  &
            ]
            
            mid = (a + b) / 2.0d0
            half = (b - a) / 2.0d0
            
            result = 0.0d0
            
            ! do i = 1, 16
            !       t = mid + half * x(i)
            !       result = result + w(i) * (f(t) + f(2.0d0*mid - t))
            ! end do

            do i = 1, 16
                t = mid + half * x(i)
                i_left = 16 + i
                i_right = 17 - i
                result = result + w(i) * (f(t, i_left) + f(2.0d0*mid - t, i_right))
            end do
            
            result = result * half
            
      end subroutine gauss_legendre_32
      
      
      ! ============================================================
      ! 64-point Gauss-Legendre (high accuracy, ~15-16 digit accuracy)
      ! ============================================================

      subroutine gauss_legendre_64(f, a, b, result)
      
            interface
            function f(t,i_point) result(y)
                  double precision, intent(in) :: t
                  integer, intent(in)          :: i_point
                  double precision             :: y
            end function
      end interface

      double precision, intent(in)  :: a, b
      double precision, intent(out) :: result
      
      double precision :: x(32), w(32)
      integer          :: i
      double precision :: mid, half, t
      double precision :: u
      integer          :: i_left, i_right
      
      ! Nodes (positive x values, symmetric about 0)
      x = [ &
            0.024350292663424432509, &
            0.072993121787799039450, &
            0.121462819296120553470, &
            0.169644420423992818037, &
            0.217423643740007084149, &
            0.264687162208767416374, &
            0.311322871990210956157, &
            0.357220158337668115950, &
            0.402270157963991603695, &
            0.446366017253464087984, &
            0.489403145707052957478, &
            0.531279464019894545658, &
            0.571895646202634034283, &
            0.611155355172393250249, &
            0.648965471254657339857, &
            0.685236313054233242564, &
            0.719881850171610826848, &
            0.752819907260531896612, &
            0.783972358943341407610, &
            0.813265315122797559741, &
            0.840629296252580362751, &
            0.865999398154092819760, &
            0.889315445995114105853, &
            0.910522137078502805756, &
            0.929569172131939575821, &
            0.946411374858402816062, &
            0.961008799652053718918, &
            0.973326827789910963741, &
            0.983336253884625956931, &
            0.991013371476744320739, &
            0.996340116771955279346, &
            0.999305041735772139456  &
      ]
      
      ! Weights for the corresponding nodes
      w = [ &
            0.048690957009139720383, &
            0.048575467441503426934, &
            0.048344762234802957211, &
            0.047999388596458307728, &
            0.047540165017830302004, &
            0.046968182816210017325, &
            0.046284796581314417287, &
            0.045491627927418144479, &
            0.044590558163756563060, &
            0.043583724529323453377, &
            0.042473515123653589607, &
            0.041262563242623528610, &
            0.039953741132720341386, &
            0.038550153178615629128, &
            0.037055128540240046040, &
            0.035472213256882383810, &
            0.033805161837141609391, &
            0.032057928354851553585, &
            0.030234657072402478867, &
            0.028339672614283483455, &
            0.026377469715054658671, &
            0.024352702568710873338, &
            0.022270173808383254159, &
            0.020134823153530209372, &
            0.017951715775697343085, &
            0.015726030476024719322, &
            0.013463047896718642598, &
            0.011168139460131128805, &
            0.008846759826363947723, &
            0.006504457968978362856, &
            0.004147033260562467635, &
            0.001783280721696432947  &
      ]
      
      mid = (a + b) / 2.0d0
      half = (b - a) / 2.0d0
      
      result = 0.0d0
      
      ! Sum over symmetric pairs: use positive x and its negative counterpart
      do i = 1, 32
            t = mid + half * x(i)
            i_left  = 32 + i
            i_right = 33 - i
            result = result + w(i) * (f(t,i_left) + f(2.0d0*mid - t,i_right))
      end do
      
      result = result * half
      
      end subroutine gauss_legendre_64
      
      
      ! ============================================================
      ! Automatic selection based on desired accuracy
      ! ============================================================
      ! subroutine gauss_legendre_auto(f, a, b, result, accuracy)
      !       interface
      !             function f(t) result(y)
      !                   double precision, intent(in) :: t
      !                   double precision :: y
      !             end function
      !       end interface
      !       double precision, intent(in) :: a, b
      !       double precision, intent(out) :: result
      !       integer, intent(in) :: accuracy  ! 1=fast, 2=normal, 3=high
            
      !       select case(accuracy)
      !       case(1)  ! Fast (16 points, ~4-5 digit accuracy)
      !             call gauss_legendre_16(f, a, b, result)
      !       case(2)  ! Normal (32 points, ~8-9 digit accuracy)
      !             call gauss_legendre_32(f, a, b, result)
      !       case(3)  ! High (64 points, ~12-13 digit accuracy)
      !             call gauss_legendre_64(f, a, b, result)
      !       case default
      !             call gauss_legendre_32(f, a, b, result)
      !       end select
            
      ! end subroutine gauss_legendre_auto
      
end module gauss_legendre_quadrature