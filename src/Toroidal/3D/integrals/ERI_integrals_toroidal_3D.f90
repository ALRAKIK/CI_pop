subroutine ERI_integral_4_function_toroidal_3D(one,two,three,four,value)
      
      use torus_init    
      use classification_ERI
      use HeavisideModule

      implicit none 

      !-----------------------------------------------------------------!

      type(ERI_function),intent(in)   :: one , two , three , four 

      double precision,intent(out)    :: value 

      integer                         ::  i , j   , k  , l
      character(LEN=2)                :: o1 , o2 , o3 , o4 
      double precision,parameter      :: pi     = 3.14159265358979323846D00
      
      double precision                :: xa  , ya , za 
      double precision                :: xb  , yb , zb
      double precision                :: xc  , yc , zc
      double precision                :: xd  , yd , zd
      double precision                :: xAB , xCD 
      double precision                :: yAB , yCD
      double precision                :: zAB , zCD
      double precision                :: xp  , xq 
      double precision                :: yp  , yq
      double precision                :: zp  , zq
      double precision                :: c1  , c2 , c3 , c4 
      double precision                :: alpha , beta , gamma , delta 
      double precision                :: const  
      double precision                :: value_s 
      double precision                :: der 
      double precision                :: mu_x, mu_y , mu_z
      double precision                :: nu_x, nu_y , nu_z
      double precision                :: mu  , nu 
      double precision                :: xpA , xpB , xqC , xqD , phix
      double precision                :: ypA , ypB , yqC , yqD , phiy
      double precision                :: zpA , zpB , zqC , zqD , phiz 
      double precision                :: const_x , const_y , const_z 
      double precision                :: inv_ax2, inv_ay2, inv_az2
      integer                         :: pattern_id, encode_orbital_pattern


      inv_ax2 = 1.0d0/(ax*ax)
      inv_ay2 = 1.0d0/(ay*ay) 
      inv_az2 = 1.0d0/(az*az)

      !-----------------------------------------------------------------!

      value = 0.d0 
            
      xa =   one%x ; ya =   one%y ; za =   one%z 
      xb =   two%x ; yb =   two%y ; zb =   two%z
      xc = three%x ; yc = three%y ; zc = three%z
      xd =  four%x ; yd =  four%y ; zd =  four%z


      xAB = xa - xb ; xCD = xc - xd
      yAB = ya - yb ; yCD = yc - yd
      zAB = za - zb ; zCD = zc - zd 

      do i = 1 , size  (one%exponent)
        alpha = one%exponent(i)
        c1    = one%coefficient(i)
        o1    = one%orbital
        do j = 1 , size  (two%exponent)
          beta  = two%exponent(j)
          c2    = two%coefficient(j)
          o2    = two%orbital

          mu_x  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ax*(XAB)))
          mu_y  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(ay*(YAB)))
          mu_z  = dsqrt(alpha**2+beta**2+2.d0*alpha*beta*cos(az*(ZAB)))

          call bary_center_toroidal(alpha,beta,xa,xb,xp)
          call bary_center_toroidal(alpha,beta,ya,yb,yp)
          call bary_center_toroidal(alpha,beta,za,zb,zp)

          mu = alpha+beta 
          
          do k = 1 , size(three%exponent)
            gamma = three%exponent(k)
            c3    = three%coefficient(k)
            o3    = three%orbital
            do l = 1 , size (four%exponent)
              delta = four%exponent(l)
              c4    = four%coefficient(l)
              o4    = four%orbital


              nu_x  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ax*(XCD)))
              nu_y  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(ay*(YCD)))
              nu_z  = dsqrt(gamma**2+delta**2+2.d0*gamma*delta*dcos(az*(ZCD)))

              call bary_center_toroidal(gamma,delta,xc,xd,xq)
              call bary_center_toroidal(gamma,delta,yc,yd,yq)
              call bary_center_toroidal(gamma,delta,zc,zd,zq)
               
              nu    = gamma + delta

              const   = (c1*c2*c3*c4) * 2.d0 /dsqrt(pi) * (Lx*Lx) * (Ly*Ly) * (Lz*Lz)  

              const_x = dexp(2.d0*( mu_x + nu_x - mu - nu )*inv_ax2)
              const_y = dexp(2.d0*( mu_y + nu_y - mu - nu )*inv_ay2)
              const_z = dexp(2.d0*( mu_z + nu_z - mu - nu )*inv_az2)

              if ( (const_x*const_y*const_z < 1.d-30 )  ) cycle

              xpA     = ax*(xp - xa) ; ypA     = ay*(yp - ya) ; zpA     = az*(zp - za)
              xpB     = ax*(xp - xb) ; ypB     = ay*(yp - yb) ; zpB     = az*(zp - zb) 
              xqC     = ax*(xq - xc) ; yqC     = ay*(yq - yc) ; zqC     = az*(zq - zc)
              xqD     = ax*(xq - xd) ; yqD     = ay*(yq - yd) ; zqD     = az*(zq - zd)
              phix    = ax*(xp - xq) ; phiy    = ay*(yp - yq) ; phiz    = az*(zp - zq)

              pattern_id = encode_orbital_pattern(o1, o2, o3, o4)
 
              call integrate_ERI_3D(pattern_id,mu,nu,mu_x,nu_x,phix,const_x,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                                     mu_y,nu_y,phiy,const_y,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                                     mu_z,nu_z,phiz,const_z,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,value_s)

              value  = value    + const * value_s

            end do
          end do 
        end do 
      end do 

end subroutine ERI_integral_4_function_toroidal_3D

subroutine integrate_ERI_3D(pattern_id,p,q,p_x,q_x,phix,const_x,xpA,xpB,xqC,xqD,xa,xb,xc,xd,xp,xq,&
                                           p_y,q_y,phiy,const_y,ypA,ypB,yqC,yqD,ya,yb,yc,yd,yp,yq,&
                                           p_z,q_z,phiz,const_z,zpA,zpB,zqC,zqD,za,zb,zc,zd,zp,zq,result)
      
      use quadpack , only : dqagi
      use iso_c_binding
      use torus_init
      use gsl_bessel_mod

      use, intrinsic :: ieee_arithmetic

      implicit none

      ! Input parameters
      double precision, intent(in)       :: p_x , q_x 
      double precision, intent(in)       :: p_y , q_y
      double precision, intent(in)       :: p_z , q_z 
      double precision, intent(in)       :: p   , q 
      double precision, intent(in)       :: phix , phiy , phiz
      double precision, intent(in)       :: xpA , xpB , xqC , xqD
      double precision, intent(in)       :: ypA , ypB , yqC , yqD
      double precision, intent(in)       :: zpA , zpB , zqC , zqD
      double precision, intent(in)       :: xa , xb , xc ,xd ,xp ,xq 
      double precision, intent(in)       :: ya , yb , yc ,yd ,yp ,yq
      double precision, intent(in)       :: za , zb , zc ,zd ,zp ,zq
      double precision, intent(in)       :: const_x , const_y , const_z
      integer         , intent(in)       :: pattern_id
      

      ! Output parameters

      double precision, intent(out)      :: result
    
      ! Local variables

      double precision,parameter         :: epsabs = 1.0e-12 , epsrel = 1.0e-10
      integer,parameter                  :: inf = 1 
      double precision,parameter         :: bound = 0.0d0
      integer, parameter                 :: limit = 50
      integer, parameter                 :: lenw = limit*4
      integer                            :: ier, iwork(limit), last, neval
      double precision                   :: abserr, work(lenw)
      integer                            :: Nmax
      double precision,parameter         :: pi   = 3.14159265358979323846D00
      double precision,parameter         :: pi2  = pi * pi


      double precision                   :: inv_ax  , inv_ay  , inv_az 
      double precision                   :: inv_ax2 , inv_ay2 , inv_az2 

      inv_ax  = 1.d0/ax 
      inv_ay  = 1.d0/ay 
      inv_az  = 1.d0/az 

      inv_ax2 = inv_ax * inv_ax 
      inv_ay2 = inv_ay * inv_ay 
      inv_az2 = inv_az * inv_az 
      
      call dqagi(f_decay , bound, inf, epsabs, epsrel, result, abserr, neval, ier,Limit,Lenw,Last,Iwork,Work)

      if (ier /= 0) then
        write(*,'(A,I8,A)') 'Error code = ', ier
      end if

      contains

      function f_decay(u) result(fu)
        double precision, intent(in) :: u
        double precision             :: fu

          if (u < 1.d-30) then 
            fu = 0.d0 
            return
          else 
            fu =  S(u)
          end if 

      end function f_decay

      double precision function S(t) result(sum)

      use gsl_bessel_mod

      implicit none
      double precision, intent(in)         :: t
      double precision                     :: AAx,  BBx, CCx
      double precision                     :: AAy,  BBy, CCy
      double precision                     :: AAz,  BBz, CCz
      double precision                     :: tol  = 1D-30
      COMPLEX(KIND=KIND(1.0D0)), PARAMETER :: I_dp = (0.0D0, 1.0D0)
      integer                              :: n  
      COMPLEX(KIND=KIND(1.0D0))            :: termAn , termBn
      double precision                     :: termC 
      double precision                     :: sum1, sum2, sum3
      COMPLEX(KIND=KIND(1.0D0))            :: term

      double precision                     :: cxpa    , sxpa
      double precision                     :: cypa    , sypa
      double precision                     :: czpa    , szpa

      double precision                     :: cxpb    , sxpb
      double precision                     :: cypb    , sypb
      double precision                     :: czpb    , szpb

      double precision                     :: cxqc    , sxqc
      double precision                     :: cyqc    , syqc
      double precision                     :: czqc    , szqc

      double precision                     :: cxqd    , sxqd
      double precision                     :: cyqd    , syqd
      double precision                     :: czqd    , szqd

      double precision                     :: c2xpab  , s2xpab
      double precision                     :: c2ypab  , s2ypab
      double precision                     :: c2zpab  , s2zpab

      double precision                     :: c2xqcd  , s2xqcd
      double precision                     :: c2yqcd  , s2yqcd
      double precision                     :: c2zqcd  , s2zqcd

      double precision , allocatable       :: bessi_AAx(:), bessi_BBx(:), prefactor_x(:)
      double precision , allocatable       :: bessi_AAy(:), bessi_BBy(:), prefactor_y(:)
      double precision , allocatable       :: bessi_AAz(:), bessi_BBz(:), prefactor_z(:)

      AAx   = inv_ax2 * 2.d0 * p_x
      BBx   = inv_ax2 * 2.d0 * q_x
      CCx   = inv_ax2 * 2.d0 * t * t 
      

      AAy   = inv_ay2 * 2.d0 * p_y
      BBy   = inv_ay2 * 2.d0 * q_y
      CCy   = inv_ay2 * 2.d0 * t * t 

      AAz   = inv_az2 * 2.d0 * p_z 
      BBz   = inv_az2 * 2.d0 * q_z 
      CCz   = inv_az2 * 2.d0 * t * t

      AAx   = max(AAx, 1.0d-30)
      BBx   = max(BBx, 1.0d-30)

      AAy   = max(AAy, 1.0d-30)
      BBy   = max(BBy, 1.0d-30)

      AAz   = max(AAz, 1.0d-30)
      BBz   = max(BBz, 1.0d-30)

      cxpa  = dcos(xpA) ; cypa  = dcos(ypA) ; czpa  = dcos(zpA) 
      sxpa  = dsin(xpA) ; sypa  = dsin(ypA) ; szpa  = dsin(zpA) 

      cxpb  = dcos(xpB) ; cypb  = dcos(ypB) ; czpb  = dcos(zpB) 
      sxpb  = dsin(xpB) ; sypb  = dsin(ypB) ; szpb  = dsin(zpB) 
      
      
      cxqc  = dcos(xqC) ; cyqc  = dcos(yqC) ; czqc  = dcos(zqC) 
      sxqc  = dsin(xqC) ; syqc  = dsin(yqC) ; szqc  = dsin(zqC) 

      cxqd  = dcos(xqD) ; cyqd  = dcos(yqD) ; czqd  = dcos(zqD)
      sxqd  = dsin(xqD) ; syqd  = dsin(yqD) ; szqd  = dsin(zqD)

      c2xpab = dcos(ax*(2.d0*xp-xA-xB)) ; c2ypab = dcos(ay*(2.d0*yp-yA-yB)) ; c2zpab = dcos(az*(2.d0*zp-zA-zB))
      s2xpab = dsin(ax*(2.d0*xp-xA-xB)) ; s2ypab = dsin(ay*(2.d0*yp-yA-yB)) ; s2zpab = dsin(az*(2.d0*zp-zA-zB))

      c2xqcd = dcos(ax*(2.d0*xq-xC-xD)) ; c2yqcd = dcos(ay*(2.d0*yq-yC-yD)) ; c2zqcd = dcos(az*(2.d0*zq-zC-zD)) 
      s2xqcd = dsin(ax*(2.d0*xq-xC-xD)) ; s2yqcd = dsin(ay*(2.d0*yq-yC-yD)) ; s2zqcd = dsin(az*(2.d0*zq-zC-zD)) 

      select case(pattern_id)
        
      ! case (0000) ! | s   s   s   s    ( 1 ) 

      ! sum1      = 0.d0
      ! sum2      = 0.d0
      ! sum3      = 0.d0
      ! sum       = 0.d0

      ! n         = 0
      ! sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
      ! Nmax      = floor(max(AAx,BBx,CCx)) + 50 
      ! do n      = 1 , Nmax
      !   termAn  = bessi_scaled(n, AAx)
      !   termBn  = bessi_scaled(n, BBx)
      !   termc   = bessi_scaled(n, CCX) 
      !   term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
      !   if (abs(term) < tol ) exit
      !     sum2 = sum2 + 2.d0 * real(term)
      ! end do

      ! n         = 0
      ! sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
      ! Nmax      = floor(max(AAy,BBy,CCy)) + 50 
      ! do n      = 1 , Nmax
      !   termAn  = bessi_scaled(n, AAy)
      !   termBn  = bessi_scaled(n, BBy)
      !   termc   = bessi_scaled(n, CCY) 
      !   term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
      !   if (abs(term) < tol ) exit
      !     sum2 = sum2 + 2.d0 * real(term)
      ! end do

      ! n         = 0
      ! sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
      ! Nmax      = floor(max(AAz,BBz,CCz)) + 50 
      ! do n      = 1 , Nmax
      !   termAn  = bessi_scaled(n, AAz)
      !   termBn  = bessi_scaled(n, BBz)
      !   termc   = bessi_scaled(n, CCZ) 
      !   term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
      !   if (abs(term) < tol ) exit
      !     sum3 = sum3 + 2.d0 * real(term)
      ! end do

      case (0000) ! | s   s   s   s    ( 1 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0001) ! | s   s   s   px   ( 2 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0002) ! | s   s   s   py   ( 3 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0003) ! | s   s   s   pz   ( 4 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0010) ! | s   s   px  s    ( 5 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0011) ! | s   s   px  px   ( 6 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0012) ! | s   s   px  py   ( 7 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0013) ! | s   s   px  pz   ( 8 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0020) ! | s   s   py  s    ( 9 ) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0021) ! | s   s   py  px   ( 10) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0022) ! | s   s   py  py   ( 11) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0023) ! | s   s   py  pz   ( 12) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0030) ! | s   s   pz  s    ( 13) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0031) ! | s   s   pz  px   ( 14) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0032) ! | s   s   pz  py   ( 15) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0033) ! | s   s   pz  pz   ( 16) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0100) ! | s   px  s   s    ( 17) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0101) ! | s   px  s   px   ( 18) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0102) ! | s   px  s   py   ( 19) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0103) ! | s   px  s   pz   ( 20) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0110) ! | s   px  px  s    ( 21) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0111) ! | s   px  px  px   ( 22) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0112) ! | s   px  px  py   ( 23) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0113) ! | s   px  px  pz   ( 24) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0120) ! | s   px  py  s    ( 25) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0121) ! | s   px  py  px   ( 26) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0122) ! | s   px  py  py   ( 27) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0123) ! | s   px  py  pz   ( 28) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0130) ! | s   px  pz  s    ( 29) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0131) ! | s   px  pz  px   ( 30) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0132) ! | s   px  pz  py   ( 31) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0133) ! | s   px  pz  pz   ( 32) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0200) ! | s   py  s   s    ( 33) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0201) ! | s   py  s   px   ( 34) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0202) ! | s   py  s   py   ( 35) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0203) ! | s   py  s   pz   ( 36) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0210) ! | s   py  px  s    ( 37) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0211) ! | s   py  px  px   ( 38) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0212) ! | s   py  px  py   ( 39) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0213) ! | s   py  px  pz   ( 40) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0220) ! | s   py  py  s    ( 41) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0221) ! | s   py  py  px   ( 42) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0222) ! | s   py  py  py   ( 43) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0223) ! | s   py  py  pz   ( 44) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0230) ! | s   py  pz  s    ( 45) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0231) ! | s   py  pz  px   ( 46) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0232) ! | s   py  pz  py   ( 47) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0233) ! | s   py  pz  pz   ( 48) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0300) ! | s   pz  s   s    ( 49) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0301) ! | s   pz  s   px   ( 50) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0302) ! | s   pz  s   py   ( 51) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0303) ! | s   pz  s   pz   ( 52) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0310) ! | s   pz  px  s    ( 53) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0311) ! | s   pz  px  px   ( 54) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0312) ! | s   pz  px  py   ( 55) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0313) ! | s   pz  px  pz   ( 56) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0320) ! | s   pz  py  s    ( 57) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0321) ! | s   pz  py  px   ( 58) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0322) ! | s   pz  py  py   ( 59) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0323) ! | s   pz  py  pz   ( 60) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0330) ! | s   pz  pz  s    ( 61) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0331) ! | s   pz  pz  px   ( 62) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0332) ! | s   pz  pz  py   ( 63) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (0333) ! | s   pz  pz  pz   ( 64) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1000) ! | px  s   s   s    ( 65) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1001) ! | px  s   s   px   ( 66) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1002) ! | px  s   s   py   ( 67) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1003) ! | px  s   s   pz   ( 68) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1010) ! | px  s   px  s    ( 69) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1011) ! | px  s   px  px   ( 70) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1012) ! | px  s   px  py   ( 71) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1013) ! | px  s   px  pz   ( 72) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1020) ! | px  s   py  s    ( 73) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1021) ! | px  s   py  px   ( 74) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1022) ! | px  s   py  py   ( 75) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1023) ! | px  s   py  pz   ( 76) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1030) ! | px  s   pz  s    ( 77) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1031) ! | px  s   pz  px   ( 78) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1032) ! | px  s   pz  py   ( 79) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1033) ! | px  s   pz  pz   ( 80) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1100) ! | px  px  s   s    ( 81) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1101) ! | px  px  s   px   ( 82) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1102) ! | px  px  s   py   ( 83) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1103) ! | px  px  s   pz   ( 84) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1110) ! | px  px  px  s    ( 85) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1111) ! | px  px  px  px   ( 86) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1112) ! | px  px  px  py   ( 87) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1113) ! | px  px  px  pz   ( 88) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1120) ! | px  px  py  s    ( 89) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1121) ! | px  px  py  px   ( 90) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1122) ! | px  px  py  py   ( 91) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1123) ! | px  px  py  pz   ( 92) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1130) ! | px  px  pz  s    ( 93) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1131) ! | px  px  pz  px   ( 94) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1132) ! | px  px  pz  py   ( 95) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1133) ! | px  px  pz  pz   ( 96) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax2 * (cxpa * cxpb * bessi_scaled(n,AAx)-c2xpab*(0.25d0*(bessi_scaled(n-2,AAx)+2.d0*bessi_scaled(n,AAx)+bessi_scaled(n+2,AAx)))+I_dp/AAx*n * s2xpab * (0.5d0*(bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1200) ! | px  py  s   s    ( 97) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1201) ! | px  py  s   px   ( 98) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1202) ! | px  py  s   py   ( 99) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1203) ! | px  py  s   pz   ( 100) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1210) ! | px  py  px  s    ( 101) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1211) ! | px  py  px  px   ( 102) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1212) ! | px  py  px  py   ( 103) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1213) ! | px  py  px  pz   ( 104) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1220) ! | px  py  py  s    ( 105) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1221) ! | px  py  py  px   ( 106) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1222) ! | px  py  py  py   ( 107) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1223) ! | px  py  py  pz   ( 108) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1230) ! | px  py  pz  s    ( 109) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1231) ! | px  py  pz  px   ( 110) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1232) ! | px  py  pz  py   ( 111) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1233) ! | px  py  pz  pz   ( 112) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1300) ! | px  pz  s   s    ( 113) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1301) ! | px  pz  s   px   ( 114) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1302) ! | px  pz  s   py   ( 115) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1303) ! | px  pz  s   pz   ( 116) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1310) ! | px  pz  px  s    ( 117) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1311) ! | px  pz  px  px   ( 118) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1312) ! | px  pz  px  py   ( 119) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1313) ! | px  pz  px  pz   ( 120) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1320) ! | px  pz  py  s    ( 121) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1321) ! | px  pz  py  px   ( 122) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1322) ! | px  pz  py  py   ( 123) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1323) ! | px  pz  py  pz   ( 124) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1330) ! | px  pz  pz  s    ( 125) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1331) ! | px  pz  pz  px   ( 126) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1332) ! | px  pz  pz  py   ( 127) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (1333) ! | px  pz  pz  pz   ( 128) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpa * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpa * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2000) ! | py  s   s   s    ( 129) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2001) ! | py  s   s   px   ( 130) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2002) ! | py  s   s   py   ( 131) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2003) ! | py  s   s   pz   ( 132) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2010) ! | py  s   px  s    ( 133) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2011) ! | py  s   px  px   ( 134) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2012) ! | py  s   px  py   ( 135) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2013) ! | py  s   px  pz   ( 136) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2020) ! | py  s   py  s    ( 137) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2021) ! | py  s   py  px   ( 138) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2022) ! | py  s   py  py   ( 139) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2023) ! | py  s   py  pz   ( 140) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2030) ! | py  s   pz  s    ( 141) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2031) ! | py  s   pz  px   ( 142) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2032) ! | py  s   pz  py   ( 143) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2033) ! | py  s   pz  pz   ( 144) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2100) ! | py  px  s   s    ( 145) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2101) ! | py  px  s   px   ( 146) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2102) ! | py  px  s   py   ( 147) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2103) ! | py  px  s   pz   ( 148) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2110) ! | py  px  px  s    ( 149) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2111) ! | py  px  px  px   ( 150) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2112) ! | py  px  px  py   ( 151) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2113) ! | py  px  px  pz   ( 152) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2120) ! | py  px  py  s    ( 153) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2121) ! | py  px  py  px   ( 154) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2122) ! | py  px  py  py   ( 155) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2123) ! | py  px  py  pz   ( 156) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2130) ! | py  px  pz  s    ( 157) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2131) ! | py  px  pz  px   ( 158) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2132) ! | py  px  pz  py   ( 159) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2133) ! | py  px  pz  pz   ( 160) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2200) ! | py  py  s   s    ( 161) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2201) ! | py  py  s   px   ( 162) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2202) ! | py  py  s   py   ( 163) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2203) ! | py  py  s   pz   ( 164) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2210) ! | py  py  px  s    ( 165) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2211) ! | py  py  px  px   ( 166) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2212) ! | py  py  px  py   ( 167) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2213) ! | py  py  px  pz   ( 168) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2220) ! | py  py  py  s    ( 169) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2221) ! | py  py  py  px   ( 170) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2222) ! | py  py  py  py   ( 171) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2223) ! | py  py  py  pz   ( 172) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2230) ! | py  py  pz  s    ( 173) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2231) ! | py  py  pz  px   ( 174) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2232) ! | py  py  pz  py   ( 175) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2233) ! | py  py  pz  pz   ( 176) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay2 * (cypa*cypb*bessi_scaled(n,AAy)- c2ypab * (0.25d0*(bessi_scaled(n-2,AAy)+2.d0*bessi_scaled(n,AAy)+bessi_scaled(n+2,AAy)))+I_dp/AAy*n * s2ypab * (0.5d0*(bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = bessi_scaled(n, AAz) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAz)
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2300) ! | py  pz  s   s    ( 177) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2301) ! | py  pz  s   px   ( 178) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2302) ! | py  pz  s   py   ( 179) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2303) ! | py  pz  s   pz   ( 180) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2310) ! | py  pz  px  s    ( 181) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2311) ! | py  pz  px  px   ( 182) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2312) ! | py  pz  px  py   ( 183) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2313) ! | py  pz  px  pz   ( 184) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2320) ! | py  pz  py  s    ( 185) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2321) ! | py  pz  py  px   ( 186) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2322) ! | py  pz  py  py   ( 187) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2323) ! | py  pz  py  pz   ( 188) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2330) ! | py  pz  pz  s    ( 189) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2331) ! | py  pz  pz  px   ( 190) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2332) ! | py  pz  pz  py   ( 191) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (2333) ! | py  pz  pz  pz   ( 192) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypa * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypa * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpb *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpb * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3000) ! | pz  s   s   s    ( 193) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3001) ! | pz  s   s   px   ( 194) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3002) ! | pz  s   s   py   ( 195) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3003) ! | pz  s   s   pz   ( 196) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3010) ! | pz  s   px  s    ( 197) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3011) ! | pz  s   px  px   ( 198) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3012) ! | pz  s   px  py   ( 199) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3013) ! | pz  s   px  pz   ( 200) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3020) ! | pz  s   py  s    ( 201) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3021) ! | pz  s   py  px   ( 202) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3022) ! | pz  s   py  py   ( 203) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3023) ! | pz  s   py  pz   ( 204) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3030) ! | pz  s   pz  s    ( 205) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3031) ! | pz  s   pz  px   ( 206) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3032) ! | pz  s   pz  py   ( 207) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3033) ! | pz  s   pz  pz   ( 208) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3100) ! | pz  px  s   s    ( 209) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3101) ! | pz  px  s   px   ( 210) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3102) ! | pz  px  s   py   ( 211) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3103) ! | pz  px  s   pz   ( 212) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3110) ! | pz  px  px  s    ( 213) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3111) ! | pz  px  px  px   ( 214) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3112) ! | pz  px  px  py   ( 215) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3113) ! | pz  px  px  pz   ( 216) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3120) ! | pz  px  py  s    ( 217) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3121) ! | pz  px  py  px   ( 218) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3122) ! | pz  px  py  py   ( 219) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3123) ! | pz  px  py  pz   ( 220) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3130) ! | pz  px  pz  s    ( 221) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3131) ! | pz  px  pz  px   ( 222) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3132) ! | pz  px  pz  py   ( 223) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3133) ! | pz  px  pz  pz   ( 224) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx))) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ax * 0.5d0 * (I_dp * cxpb * (bessi_scaled(n-1,AAx) - bessi_scaled(n+1,AAx)) + sxpb * (bessi_scaled(n-1,AAx)+bessi_scaled(n+1,AAx)))
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3200) ! | pz  py  s   s    ( 225) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3201) ! | pz  py  s   px   ( 226) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3202) ! | pz  py  s   py   ( 227) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3203) ! | pz  py  s   pz   ( 228) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3210) ! | pz  py  px  s    ( 229) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3211) ! | pz  py  px  px   ( 230) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3212) ! | pz  py  px  py   ( 231) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3213) ! | pz  py  px  pz   ( 232) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3220) ! | pz  py  py  s    ( 233) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3221) ! | pz  py  py  px   ( 234) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3222) ! | pz  py  py  py   ( 235) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3223) ! | pz  py  py  pz   ( 236) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3230) ! | pz  py  pz  s    ( 237) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3231) ! | pz  py  pz  px   ( 238) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3232) ! | pz  py  pz  py   ( 239) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3233) ! | pz  py  pz  pz   ( 240) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy))) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = inv_ay * 0.5d0 * (I_dp * cypb * (bessi_scaled(n-1,AAy)-bessi_scaled(n+1,AAy)) + sypb * (bessi_scaled(n-1,AAy)+bessi_scaled(n+1,AAy)))
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az * 0.5d0 * (I_dp * czpa *(bessi_scaled(n-1,AAz)-bessi_scaled(n+1,AAz))+ szpa * (bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3300) ! | pz  pz  s   s    ( 241) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3301) ! | pz  pz  s   px   ( 242) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3302) ! | pz  pz  s   py   ( 243) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3303) ! | pz  pz  s   pz   ( 244) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3310) ! | pz  pz  px  s    ( 245) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3311) ! | pz  pz  px  px   ( 246) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax2 * (cxqc*cxqd*bessi_scaled(n,BBx)-c2xqcd*(0.25d0*(bessi_scaled(n-2,BBx)+2.d0*bessi_scaled(n,BBx)+bessi_scaled(n+2,BBx)))-I_dp/BBx*n * s2xqcd * (0.5d0*(bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3312) ! | pz  pz  px  py   ( 247) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3313) ! | pz  pz  px  pz   ( 248) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqc * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqc * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3320) ! | pz  pz  py  s    ( 249) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3321) ! | pz  pz  py  px   ( 250) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3322) ! | pz  pz  py  py   ( 251) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay2 * (cyqc*cyqd*bessi_scaled(n,BBy)- c2yqcd * (0.25d0*(bessi_scaled(n-2,BBy)+2.d0*bessi_scaled(n,BBy)+bessi_scaled(n+2,BBy)))-I_dp/BBy*n * s2yqcd * (0.5d0*(bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * bessi_scaled(n, BBz) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = bessi_scaled(n, BBz)
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3323) ! | pz  pz  py  pz   ( 252) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqc * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqc * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqd * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqd * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3330) ! | pz  pz  pz  s    ( 253) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3331) ! | pz  pz  pz  px   ( 254) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx))) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = inv_ax * 0.5d0 * (-I_dp * cxqd * (bessi_scaled(n-1,BBx)-bessi_scaled(n+1,BBx)) + sxqd * (bessi_scaled(n-1,BBx)+bessi_scaled(n+1,BBx)))
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3332) ! | pz  pz  pz  py   ( 255) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy))) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = inv_ay * 0.5d0 * (-I_dp * cyqd * (bessi_scaled(n-1,BBy)-bessi_scaled(n+1,BBy)) + syqd * (bessi_scaled(n-1,BBy)+bessi_scaled(n+1,BBy)))
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az * 0.5d0 * (-I_dp * czqc * (bessi_scaled(n-1,BBz)-bessi_scaled(n+1,BBz)) + szqc * (bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do


case (3333) ! | pz  pz  pz  pz   ( 256) 

sum1      = 0.d0
sum2      = 0.d0
sum3      = 0.d0
sum       = 0.d0

n         = 0
sum1      = bessi_scaled(n, AAx) * bessi_scaled(n, BBx) * bessi_scaled(n, CCx) 
Nmax      = floor(max(AAx,BBx,CCx)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAx)
  termBn  = bessi_scaled(n, BBx)
  termc   = bessi_scaled(n, CCX) 
  term    = exp(I_dp*dble(n)*phix) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum1 = sum1 + 2.d0 * real(term)
end do

n         = 0
sum2      = bessi_scaled(n, AAy) * bessi_scaled(n, BBy) * bessi_scaled(n, CCY) 
Nmax      = floor(max(AAy,BBy,CCy)) + 50 
do n      = 1 , Nmax
  termAn  = bessi_scaled(n, AAy)
  termBn  = bessi_scaled(n, BBy)
  termc   = bessi_scaled(n, CCY) 
  term    = exp(I_dp*dble(n)*phiy) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum2 = sum2 + 2.d0 * real(term)
end do

n         = 0
sum3      = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz)))) * inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz)))) * bessi_scaled(n, CCZ) 
Nmax      = floor(max(AAz,BBz,CCz)) + 50 
do n      = 1 , Nmax
  termAn  = inv_az2 * (czpa*czpb*bessi_scaled(n,AAz)-c2zpab*(0.25d0*(bessi_scaled(n-2,AAz)+2.d0*bessi_scaled(n,AAz)+bessi_scaled(n+2,AAz)))+I_dp/AAz*n * s2zpab * (0.5d0*(bessi_scaled(n-1,AAz)+bessi_scaled(n+1,AAz))))
  termBn  = inv_az2 * (czqc*czqd*bessi_scaled(n,BBz)-c2zqcd*(0.25d0*(bessi_scaled(n-2,BBz)+2.d0*bessi_scaled(n,BBz)+bessi_scaled(n+2,BBz)))-I_dp/BBz*n*s2zqcd*(0.5d0*(bessi_scaled(n-1,BBz)+bessi_scaled(n+1,BBz))))
  termc   = bessi_scaled(n, CCZ) 
  term    = exp(I_dp*dble(n)*phiz) * termC * termAn * termBn
  if (abs(term) < tol ) exit
    sum3 = sum3 + 2.d0 * real(term)
end do




      case default
        sum = 0.0d0
      end select


      sum = sum1 * sum2 * sum3 * const_x * const_y * const_z

      end function S

end subroutine integrate_ERI_3D