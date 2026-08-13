subroutine Hartree_Fock(nBas,nO,ERI_size,S,T,V,Hc,ERI,X,Enuc,EHF,e,c,coeff_SAO &
                       ,n_alpha, n_beta )

      use keywords
      use constants_module
      implicit none

      ! input ! 

      integer,intent(in)            :: nBas , ERI_size
  
      integer,intent(in)            ::                          nO
      double precision,intent(in)   ::                 S(nBas,nBas)
      double precision,intent(in)   ::                 T(nBas,nBas)
      double precision,intent(in)   ::                 V(nBas,nBas)
      double precision,intent(in)   ::                Hc(nBas,nBas)
      double precision,intent(in)   ::                 X(nBas,nBas)
      double precision,intent(in)   :: ERI(ERI_size,nBas,nBas,nBas)

      integer        ,intent(in)    :: n_alpha , n_beta 

      double precision,intent(in)   :: ENuc


      ! local ! 

      ! output ! 

      double precision,intent(out)  :: EHF 
      double precision,intent(out)  :: e(nBas)
      double precision,intent(out)  :: c(nBas,nBas)
      complex(dpc),    intent(out)  :: coeff_SAO(nBas,nBas)


      if (c_SAO) then 
        call RHF_SAO(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,coeff_SAO)
        go to 100
      end if 

      if (c_UHF) then 
        call UHF(nBas,c_details,n_alpha,n_beta,S,T,V,Hc,ERI,X,Enuc,EHF,e,c)
        go to 100
      end if

      call RHF(nBas,nO,S,T,V,Hc,ERI,X,Enuc,EHF,e,c)


100 continue

end subroutine 