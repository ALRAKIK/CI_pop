subroutine print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,Ex,EHF)

      ! Print one- and two-electron energies and other stuff for RHF calculation

      use files 
      use keywords
      implicit none

      integer,intent(in)                 :: nBas,nO
      double precision,intent(in)        :: e(nBas),c(nBas,nBas),ENuc,ET,EV,EJ,Ex,EHF

      double precision                   :: e_sorted(nBas)
      integer                            :: HOMO,LUMO , i , o
      double precision                   :: Gap
      double precision,parameter         :: autoev = 27.21138602d0

      ! HOMO and LUMO

      HOMO = nO
      LUMO = HOMO + 1
      
      if (nBas > nO) then
        e_sorted = e
        call sort_array(e_sorted, nBas)
        Gap = e_sorted(nO+1) - e_sorted(nO)
      else
        Gap = 0.d0
      endif


      ! Dump results

      write(outfile,*)
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32)')           ' Summary              '
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' One-electron energy      ',ET + EV
      write(outfile,'(A36,1X,F24.16)')      ' Kinetic      energy      ',ET
      write(outfile,'(A36,1X,F24.16)')      ' Potential    energy      ',EV
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Two-electron energy      ',EJ + Ex
      write(outfile,'(A36,1X,F24.16)')      ' Coulomb      energy      ',EJ
      write(outfile,'(A36,1X,F24.16)')      ' Exchange     energy      ',Ex
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Electronic   energy      ',EHF
      write(outfile,'(A36,1X,F24.16)')      ' Nuclear   repulsion      ',ENuc
      write(outfile,'(A36,1X,F24.16)')      ' Hartree-Fock energy      ',EHF + ENuc
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F24.16)')      ' Potential energy        :',(EV+EJ+Ex+ENuc)
      write(outfile,'(A36,1X,F24.16)')      ' Kinetic   energy        :', ET
      write(outfile,'(A36,1X,F24.16)')      ' Virial    theorem (-V/T):',-(EV+EJ+Ex+ENuc)/ET      
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO      energy (au):',e_sorted(HOMO)
      write(outfile,'(A36,1X,F13.6)')       ' HF LUMO      energy (au):',e_sorted(LUMO)
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO-LUMO gap    (au):',Gap
      write(outfile,'(A36,1X,F13.6)')       ' HF HOMO-LUMO gap    (eV):',Gap*autoev
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,*)

      ! Print results

      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') ' Hartree-Fock orbital energies         '
      write(outfile,'(A50)') '---------------------------------------'
      call matout(nBas,1,e)
      write(outfile,*)

      if (c_Orbitals) then 

      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"
      write(orbfile,'(A)') "|                                 The molecular orbitals                                  |"
      write(orbfile,'(A)') "|-----------------------------------------------------------------------------------------|"

      write(orbfile,'(a)') ""
      write(orbfile,'(A)') "|---------------------------------- The orbital coeff c ----------------------------------|"
      write(orbfile,'(a)') ""
      write(orbfile,'(17x,1000(i3,15x))') (i,i=1,nBas)
      do i = 1 , nBas
        write(orbfile,'(i3,6x,1000(f16.10,2x))') i ,  (c(i,o),o=1,nBas)
      end do 
      write(orbfile,'(a)') ""

      Flush(orbfile)
    
      end if

end subroutine print_RHF