integer function count_px_orbitals(o1, o2, o3, o4)
    implicit none
    character(len=*), intent(in) :: o1, o2, o3, o4
    
    count_px_orbitals = 0
    if (o1 == "px") count_px_orbitals = count_px_orbitals + 1
    if (o2 == "px") count_px_orbitals = count_px_orbitals + 1
    if (o3 == "px") count_px_orbitals = count_px_orbitals + 1
    if (o4 == "px") count_px_orbitals = count_px_orbitals + 1
end function count_px_orbitals