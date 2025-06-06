subroutine FCI(Hc,ERI,nbas)

      use files

      implicit none 

      integer ,intent(in)           :: nbas
      double precision, intent(in)  :: Hc(nbas,nbas)
      double precision, intent(in)  :: ERI(nbas,nbas,nbas,nbas)


      double precision              :: e1 , e2 , ehf , delta , efci


      e1  = HC(1,1)+ERI(1,1,1,1)
      e2  = HC(2,2)+2.d0*ERI(1,1,2,2)-ERI(1,2,1,2)

      ehf = 2.d0*HC(1,1)+ERI(1,1,1,1)

      delta=2.d0*(e2-e1)+ERI(1,1,1,1)+ERI(2,2,2,2)-4.d0*ERI(1,1,2,2)+2.d0*ERI(1,2,1,2)

      delta = delta/2.d0

      efci = ehf + delta - dsqrt(delta*delta+4.d0*ERI(1,2,1,2)*ERI(1,2,1,2))

      call header("FCI Calculation",-1)

      write(outfile,*) "EPS   :" , e1 , e2
      write(outfile,*) "" 
      write(outfile,*) "EHF   :" , ehf
      write(outfile,*) ""
      write(outfile,*) "EFCI  :" , efci

end subroutine FCI