subroutine get_charge(number_of_atoms,type,charge)

      implicit none
      
      integer               :: i 
      integer               :: number_of_atoms
      character(len=2)      :: type(number_of_atoms)
      integer               :: charge(number_of_atoms)
      

      do i = 1 , number_of_atoms
        if (type(i) == "H" ) charge(i) = 1
        if (type(i) == "He") charge(i) = 2
        if (type(i) == "Li") charge(i) = 3
        if (type(i) == "Be") charge(i) = 4
        if (type(i) == "B" ) charge(i) = 5
        if (type(i) == "C" ) charge(i) = 6
        if (type(i) == "N" ) charge(i) = 7
        if (type(i) == "O" ) charge(i) = 8
        if (type(i) == "F" ) charge(i) = 9
        if (type(i) == "Ne") charge(i) = 10
        if (type(i) == "Na") charge(i) = 11
        if (type(i) == "Mg") charge(i) = 12
        if (type(i) == "Al") charge(i) = 13
        if (type(i) == "Si") charge(i) = 14
        if (type(i) == "P" ) charge(i) = 15
        if (type(i) == "S" ) charge(i) = 16
        if (type(i) == "Cl") charge(i) = 17
        if (type(i) == "Ar") charge(i) = 18
        if (type(i) == "K" ) charge(i) = 19
        if (type(i) == "Ca") charge(i) = 20
        if (type(i) == "Sc") charge(i) = 21
        if (type(i) == "Ti") charge(i) = 22
        if (type(i) == "V" ) charge(i) = 23
        if (type(i) == "Cr") charge(i) = 24
        if (type(i) == "Mn") charge(i) = 25
        if (type(i) == "Fe") charge(i) = 26
        if (type(i) == "Co") charge(i) = 27
        if (type(i) == "Ni") charge(i) = 28
        if (type(i) == "Cu") charge(i) = 29
        if (type(i) == "Zn") charge(i) = 30
        if (type(i) == "Ga") charge(i) = 31
        if (type(i) == "Ge") charge(i) = 32
        if (type(i) == "As") charge(i) = 33
        if (type(i) == "Se") charge(i) = 34
        if (type(i) == "Br") charge(i) = 35
        if (type(i) == "Kr") charge(i) = 36
        if (type(i) == "GO") charge(i) = 0
    end do






end subroutine