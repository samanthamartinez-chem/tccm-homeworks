program input 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DECLARATION AREA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer(8), parameter :: n_max = 10000
real(8) :: positions(n_max, 3), cell_size
real(8), parameter :: mass_Ar = 39.95 ! Ar mass in a.u.
character(len=100) :: filename
integer :: natoms, i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAIN CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prompt the user to introduce the number of Ar atoms in the cell
! The number of atoms must be in the range between 1 and n_max
write(6, *) "Introduce the number of Ar atoms in the cell"
read(5, *) natoms

! Warning if the introduced number of particles is not valid
if (natoms .le. 0 .or. natoms .gt. n_max) then
    write(6, *) "The number of particles must be between 1 and ", n_max
    write(6, *) "Introduce a valid number of particles"
    read(5, *) natoms
end if

! Prompt the user for the simulation cubic box size
write(6, *) "Eneter the size of the cubic simulation box in nm:"
read(5, *) cell_size

! Prompt the user for the name of the input 
write(6, *) "Enter the name of the input (recommended: inp.txt)"
read(5, *) filename

! Call subroutines for random number generation
call random_seed() ! Initialize the random number generation
! Random numbers are between 0 and 1
call random_number(positions) ! Random positions (x, y, z) for each particles

! Rescale of the positions to fit them in the simulation box
! Normalized positions rescaled using the cell_size chosen by the user
positions = positions * cell_size

! Open a file to write the input data
open(16, file=filename, status="replace", action="write")
write(16, *) natoms

!  Write in the file the x, y, z coordinates + the mass of the atom
do i = 1, natoms
    write(16, '(3F12.6, F8.3)') (positions(i, j), j = 1, 3), mass_Ar
end do

close(16)

end program input
