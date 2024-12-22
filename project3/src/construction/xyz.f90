program xyz_file

implicit none

contains

! Subroutine to create the format of the xyz file
subroutine write_xyz(step, Natoms, coord, kinetic_enrgy, potential_energ
y, total_energy)
implicit none
integer, intent(in) :: step, Natoms
double precision, intent(in) :: coord(Natoms, 3), kinetic_enrgy, potential_energy, total_energy
integer :: i

! Unit 10 bc later in the main we will open a file 
write(10,*) Natoms
write(10,*) "STEP: ", step, "T: ", kinetic_energy, "V: ", potential_energy, "E: ", total_energy

do i = 1, Natoms
    write(10, '(A,3F12.6)') "Ar", coord(i, 1), coord(i, 2), coord(i, 3)
end do

end subroutine write_xyz

end program xyz_file
