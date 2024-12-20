!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: coord(:,:)
double precision, allocatable :: mass(:)
! Variables for the LJ calculation
double precision, allocatable :: distance(:,:)
double precision :: potential_energy
double precision, parameter :: epsilon=0.0661d0 ! in j/mol
double precision, parameter :: sigma=0.3345d0 ! in nm
character(len=100)::input_file
integer :: Natoms,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Execution area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Define the name of the file to read
input_file="inp.txt"

!Call the function to determine the number of atoms
Natoms=read_Natoms(input_file)

!Call the subroutine that gets the coordinates and the mass of the atoms
call read_molecule(input_file,Natoms,coord,mass)
write(*,*) mass(2),coord(1,1),coord(22,2)

! Alocates memory for the distance matrix 
call allocate_distance(Natoms, distance)

! Computes the distance and LJ potential energy
call compute_distance_potential(Natoms, coord, distance, sigma, epsilon, potential_energy)

! To confirm everything ok so far
write(6,*) "The total LJ potential energy is: ", potential_energy


! Deallocate all we allocated 
deallocate(coord)
deallocate(mass)
deallocate(distance)

contains

function get_dimension(input_file)
        implicit none
        character(len=100),intent(in)::input_file
        character(len=100)::line
        integer :: pos_x
        integer,intent(out) :: n,m
        open(unit=10, file=input_file, status="old", action="read")
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) line
        pos_x=index(line,'x')
        if (pos_x>0) then
                read(line(35:pos_x-1),*) n
                read(line(pos_x+1:),*) m
        end if 
        close(unit=10)
        end function get_dimension



!subroutine read_matrix(input_file,n,m,row,column,matrixvalue)
!        implicit none
!        character(len=100),intent(in)::input_file
!        double precision,intent(out),allocatable::row(:),column(:),matrixvalue(:)
!        integer,intent(out) :: n,m
!        open(unit=10, file=input_file, status="old", action="read")
        
!        read(10,*)

        
!        integer function read_Natoms(input_file) result(Natoms)
!        !This function opens the file, takes the first line and save it in Natoms
!        implicit none
        character(len=*), intent(in)::input_file
        open(unit=10, file=input_file, status="old", action="read")
        read(10,*) Natoms
        close(unit=10)
end function read_Natoms

subroutine read_molecule(input_file,Natoms,coord,mass)
        !This subroutine allocates the matrix coord where the coordinates of each atom are saved
        !Aldo allocates the vector where the masses of the atoms are saved
        implicit none
        character(len=100),intent(in)::input_file
        integer,intent(in)::Natoms
        integer :: n,i_stat,i,m
        double precision,intent(out),allocatable::coord(:,:)
        !character(len=2),intent(out),allocatable::mass(:)
        double precision,intent(out),allocatable::mass(:)
        n=Natoms
        m=3
        
        allocate(coord(n,m),stat=i_stat)
        allocate(mass(n),stat=i_stat)
        if (i_stat /= 0) then
                print *, "Memory allocation failed"
                stop
        end if
        open(unit=10, file=input_file, status="old", action="read")
        read(10,*)
        do i=1,n
                read(10,*) coord(i,1), coord(i,2), coord(i,3),mass(i)
        end do
        close(10)
end subroutine read_molecule

! Subroutine to allocate the distance matrix

subroutine allocate_distance(Natoms, distance)
    implicit none
    integer, intent(in) :: Natoms
    double precision, allocatable :: distance(:,:)
    integer :: i_stat

    allocate(distance(Natoms, Natoms), stat=i_stat)
        if (i_stat /= 0) then
            write(6,*) "Memory allocation for distance failed"
            stop
        end if
end subroutine allocate_distance
! Subroutine to compute the distances and the LJ potential energy

subroutine compute_distance_potential(Natoms, coord, distance, sigma, epsilon, potential_energy)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coord(Natoms, 3)
    double precision, intent(in) :: sigma, epsilon
    double precision, intent(out) :: distance(Natoms, Natoms)
    double precision, intent(out) :: potential_energy
    integer :: i, j
    double precision :: r

    ! Start Counter for V at 0
    potential_energy = 0.0d0
    do i=1, Natoms
        do j = i+1, Natoms
            r = sqrt(sum((coord(i,:) - coord(j,:))**2))
            distance(i, j) = r
            distance(j, i) = r
            potential_energy = potential_energy + lj_potential(r, sigma, epsilon)
        end do
    end do
end subroutine compute_distance_potential

double precision function lj_potential(r, sigma, epsilon)
    implicit none
    double precision, intent(in) :: r, sigma, epsilon
    lj_potential = 4.0d0 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
end function lj_potential
    






!
END PROGRAM readfile
