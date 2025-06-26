subroutine FCI(Hc,ERI,nbas,E_nuc)

      use files

      implicit none 

      integer ,intent(in)           :: nbas
      double precision, intent(in)  :: Hc(nbas,nbas)
      double precision, intent(in)  :: ERI(nbas,nbas,nbas,nbas)
      double precision, intent(in)  :: E_nuc
      double precision              :: e1 , e2 , ehf , delta , efci , etrip

      e1  = Hc(1,1)+ERI(1,1,1,1)
      e2  = Hc(2,2)+2.d0*ERI(1,1,2,2)-ERI(1,2,1,2)

      ehf = 2.d0*Hc(1,1)+ERI(1,1,1,1)

      delta=2.d0*(e2-e1)+ERI(1,1,1,1)+ERI(2,2,2,2)-4.d0*ERI(1,1,2,2)+2.d0*ERI(1,2,1,2)

      delta = delta/2.d0

      efci = ehf + delta - dsqrt(delta*delta+ERI(1,2,1,2)*ERI(1,2,1,2))

      etrip = Hc(1,1) + Hc(2,2) + ERI(1,1,2,2) - ERI(1,2,1,2)


      call header("FCI Calculation",-1)

      write(outfile,'(a,E16.10,4x,E16.10)') "EPS              " , e1 , e2
      write(outfile,*) "" 
      write(outfile,'(a,E16.10)')        "EHF              " , ehf + E_nuc
      write(outfile,*) ""
      write(outfile,'(a,E16.10)')        "EFCI(Ground)     " , efci + E_nuc
      write(outfile,*) ""
      write(outfile,'(a,E16.10)')        "E_corr           " , delta - dsqrt(delta*delta+ERI(1,2,1,2)*ERI(1,2,1,2))
      write(outfile,*) ""
      write(outfile,'(a,E16.10)')        "EFCI(Triplet)    " , etrip + E_nuc
      write(outfile,*) ""


end subroutine FCI