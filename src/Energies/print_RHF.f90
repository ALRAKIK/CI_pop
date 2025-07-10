subroutine print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,Ex,EHF)

      ! Print one- and two-electron energies and other stuff for RHF calculation

      use files 
      implicit none

      integer,intent(in)                 :: nBas,nO
      double precision,intent(in)        :: e(nBas),c(nBas,nBas),ENuc,ET,EV,EJ,Ex,EHF

      integer                            :: HOMO,LUMO
      double precision                   :: Gap
      double precision,parameter         :: autoev = 27.21138602d0

      ! HOMO and LUMO

      HOMO = nO
      LUMO = HOMO + 1
      
      Gap = e(LUMO) - e(HOMO)


      ! Dump results

      write(outfile,*)
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32)')           ' Summary              '
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
      write(outfile,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
      write(outfile,'(A32,1X,F16.10)') ' Potential    energy  ',EV
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + Ex
      write(outfile,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
      write(outfile,'(A32,1X,F16.10)') ' Exchange     energy  ',Ex
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A32,1X,F16.10)') ' Electronic   energy  ',EHF
      write(outfile,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
      write(outfile,'(A32,1X,F16.10)') ' Hartree-Fock energy  ',EHF + ENuc
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,F13.6)')     ' HF HOMO      energy (au):',e(HOMO)
      write(outfile,'(A36,F13.6)')     ' HF LUMO      energy (au):',e(LUMO)
      write(outfile,'(A36,F13.6)')     ' HF HOMO-LUMO gap    (au):',Gap
      write(outfile,'(A36,F13.6)')     ' HF HOMO-LUMO gap    (eV):',Gap*autoev
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,'(A36,F13.6)')     ' Potential energy        :',(EV+EJ+Ex+ENuc)
      write(outfile,'(A36,F13.6)')     ' Kinetic   energy        :', ET
      write(outfile,'(A36,F13.6)')     ' Virial    theorem (-V/T):',-(EV+EJ+Ex+ENuc)/ET
      write(outfile,'(A50)')           '---------------------------------------'
      write(outfile,*)

      ! Print results

      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') 'Hartree-Fock orbital coefficients      '
      write(outfile,'(A50)') '---------------------------------------'
      call matout(nBas,nBas,C)
      write(outfile,*)
      write(outfile,'(A50)') '---------------------------------------'
      write(outfile,'(A50)') ' Hartree-Fock orbital energies         '
      write(outfile,'(A50)') '---------------------------------------'
      call matout(nBas,1,e)
      write(outfile,*)

end subroutine print_RHF