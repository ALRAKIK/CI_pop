subroutine Title()

      use files 

      implicit none 

      character(len=8)  :: date
      character(len=10) :: time

      ! --------------------------------------------------------------- !

      write(outfile,*) 
      write(outfile,*) 

      write(outfile,"(8x,a)") "__| |___________________________________________________| |__"
      write(outfile,"(8x,a)") "__   ___________________________________________________   __"
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "  | |  ______     __        ______   ______     ______  | |  "
      write(outfile,"(8x,a)") "  | | /\  ___\   /\ \      /\  == \ /\  __ \   /\  == \ | |  "
      write(outfile,"(8x,a)") "  | | \ \ \____  \ \ \     \ \  _-/ \ \ \/\ \  \ \  _-/ | |  "
      write(outfile,"(8x,a)") "  | |  \ \_____\  \ \_\     \ \_\    \ \_____\  \ \_\   | |  "
      write(outfile,"(8x,a)") "  | |   \/_____/   \/_/      \/_/     \/_____/   \/_/   | |  "
      write(outfile,"(8x,a)") "  | |                                                   | |  "
      write(outfile,"(8x,a)") "__| |___________________________________________________| |__"
      write(outfile,"(8x,a)") "__   ___________________________________________________   __"
      write(outfile,"(8x,a)") "  | |                                                   | |  "

      write(outfile,'(a)')
      write(outfile,'(a)')
      write(outfile,'(a)')  "      CI_POP: Calculate Integrals, is a program to calculate   the" 
      write(outfile,'(a)')  "              molecular integrals using Clifford Toroidal Gaussian"
      write(outfile,'(a)')  "              functions for advanced quantum chemistry computations"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ", "Amer Alrakik        Email  : alrakikamer@gmail.com"
      write(outfile,'(2a)') "      ", "                    Github : https://github.com/ALRAKIK/CI_pop"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ", "Contributors        Arjan Berger"
      write(outfile,'(2a)') "      ", "                    Stefano Evangelisti"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  "                          Version: 0.0.3 | Quantum Chemistry | 2026"
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  ""


      call date_and_time(date, time)

      write(outfile,'(6x,10a)') "Executing Date : ", date(7:8), "/", date(5:6), "/", date(1:4)
      write(outfile,'(6x,10a)') "Executing time : ", time(1:2), ":", time(3:4), ":", time(5:6)
      write(outfile,'(a)')  ""

      FLUSH(outfile)

      ! --------------------------------------------------------------- !

end subroutine Title