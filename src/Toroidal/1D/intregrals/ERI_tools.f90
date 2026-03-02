integer function encode_orbital_pattern(o1, o2, o3, o4)

      implicit none
      character(len=*), intent(in) :: o1, o2, o3, o4
      integer                      :: code1, code2, code3, code4
    
      ! Convert orbital types to numeric codes
      ! s-orbitals: s=0
      ! p-orbitals: px=1, py=2, pz=3
      ! d-orbitals: dxy=4, dxz=5, dyz=6, dx2=7, dy2=8, dz2=9

      ! First orbital

      code1 = 0 ; code2 = 0 ; code3 = 0 ; code4 = 0 

      ! s orbital
      if (trim(adjustl(o1)) == "s") then
          code1 = 0
      ! p orbitals
      else if (trim(adjustl(o1)) == "px") then
          code1 = 1
      else if (trim(adjustl(o1)) == "py") then
          code1 = 2
      else if (trim(adjustl(o1)) == "pz") then
          code1 = 3
      ! d orbitals
      else if (trim(adjustl(o1)) == "dxy") then
          code1 = 4
      else if (trim(adjustl(o1)) == "dxz") then
          code1 = 5
      else if (trim(adjustl(o1)) == "dyz") then
          code1 = 6
      else if (trim(adjustl(o1)) == "dx2") then
          code1 = 7
      else if (trim(adjustl(o1)) == "dy2") then
          code1 = 8
      else if (trim(adjustl(o1)) == "dz2") then
          code1 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if
      
      ! s orbital
      if (trim(adjustl(o2)) == "s") then
          code2 = 0
      ! p orbitals
      else if (trim(adjustl(o2)) == "px") then
          code2 = 1
      else if (trim(adjustl(o2)) == "py") then
          code2 = 2
      else if (trim(adjustl(o2)) == "pz") then
          code2 = 3
      ! d orbitals
      else if (trim(adjustl(o2)) == "dxy") then
          code2 = 4
      else if (trim(adjustl(o2)) == "dxz") then
          code2 = 5
      else if (trim(adjustl(o2)) == "dyz") then
          code2 = 6
      else if (trim(adjustl(o2)) == "dx2") then
          code2 = 7
      else if (trim(adjustl(o2)) == "dy2") then
          code2 = 8
      else if (trim(adjustl(o2)) == "dz2") then
          code2 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if

      ! s orbital
      if (trim(adjustl(o3)) == "s") then
          code3 = 0
      ! p orbitals
      else if (trim(adjustl(o3)) == "px") then
          code3 = 1
      else if (trim(adjustl(o3)) == "py") then
          code3 = 2
      else if (trim(adjustl(o3)) == "pz") then
          code3 = 3
      ! d orbitals
      else if (trim(adjustl(o3)) == "dxy") then
          code3 = 4
      else if (trim(adjustl(o3)) == "dxz") then
          code3 = 5
      else if (trim(adjustl(o3)) == "dyz") then
          code3 = 6
      else if (trim(adjustl(o3)) == "dx2") then
          code3 = 7
      else if (trim(adjustl(o3)) == "dy2") then
          code3 = 8
      else if (trim(adjustl(o3)) == "dz2") then
          code3 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if

      ! s orbital
      if (trim(adjustl(o4)) == "s") then
          code4 = 0
      ! p orbitals
      else if (trim(adjustl(o4)) == "px") then
          code4 = 1
      else if (trim(adjustl(o4)) == "py") then
          code4 = 2
      else if (trim(adjustl(o4)) == "pz") then
          code4 = 3
      ! d orbitals
      else if (trim(adjustl(o4)) == "dxy") then
          code4 = 4
      else if (trim(adjustl(o4)) == "dxz") then
          code4 = 5
      else if (trim(adjustl(o4)) == "dyz") then
          code4 = 6
      else if (trim(adjustl(o4)) == "dx2") then
          code4 = 7
      else if (trim(adjustl(o4)) == "dy2") then
          code4 = 8
      else if (trim(adjustl(o4)) == "dz2") then
          code4 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if
    
      ! Create unique pattern ID

      encode_orbital_pattern = code1*1000 + code2*100 + code3*10 + code4

end function


integer function encode_orbital_pattern_AO(o1, o2)
      implicit none
      character(len=*), intent(in) :: o1, o2
      integer                      :: code1, code2
    
      ! Convert orbital types to numeric codes
      ! s-orbitals: s=0
      ! p-orbitals: px=1, py=2, pz=3
      ! d-orbitals: dxy=4, dxz=5, dyz=6, dx2=7, dy2=8, dz2=9

      ! First orbital

      code1 = 0
      code2 = 0

      ! s orbital
      if (trim(adjustl(o1)) == "s") then
          code1 = 0
      ! p orbitals
      else if (trim(adjustl(o1)) == "px") then
          code1 = 1
      else if (trim(adjustl(o1)) == "py") then
          code1 = 2
      else if (trim(adjustl(o1)) == "pz") then
          code1 = 3
      ! d orbitals
      else if (trim(adjustl(o1)) == "dxy") then
          code1 = 4
      else if (trim(adjustl(o1)) == "dxz") then
          code1 = 5
      else if (trim(adjustl(o1)) == "dyz") then
          code1 = 6
      else if (trim(adjustl(o1)) == "dx2") then
          code1 = 7
      else if (trim(adjustl(o1)) == "dy2") then
          code1 = 8
      else if (trim(adjustl(o1)) == "dz2") then
          code1 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if
      
      
      ! Second orbital
      if (trim(adjustl(o2)) == "s") then
          code2 = 0
      else if (trim(adjustl(o2)) == "px") then
          code2 = 1
      else if (trim(adjustl(o2)) == "py") then
          code2 = 2
      else if (trim(adjustl(o2)) == "pz") then
          code2 = 3
      else if (trim(adjustl(o2)) == "dxy") then
          code2 = 4
      else if (trim(adjustl(o2)) == "dxz") then
          code2 = 5
      else if (trim(adjustl(o2)) == "dyz") then
          code2 = 6
      else if (trim(adjustl(o2)) == "dx2") then
          code2 = 7
      else if (trim(adjustl(o2)) == "dy2") then
          code2 = 8
      else if (trim(adjustl(o2)) == "dz2") then
          code2 = 9
      else 
        Print*, "The user have an orbital not implemented yet"
        stop
      end if
      
      ! Create unique pattern ID
      ! Using multiplication by 100 to accommodate two-digit codes
      encode_orbital_pattern_AO = code1*10 + code2

end function