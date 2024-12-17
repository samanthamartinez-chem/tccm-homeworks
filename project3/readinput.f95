!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: coord(:,:)
double precision, allocatable :: mass(:)
character(len=100)::input_file
integer :: Natoms,n
!
!----------Execution area-----------
!Define the name of the file to read
input_file="inp.txt"
!Call the function to determine the number of atoms
Natoms=read_Natoms(input_file)
!Call the subroutine that gets the coordinates and the mass of the atoms
call read_molecule(input_file,Natoms,coord,mass)
!We dealocate the matrix and vector
write(*,*) mass(2),coord(1,1),coord(22,2)

deallocate(coord)
deallocate(mass)

contains

integer function read_Natoms(input_file) result(Natoms)
        !This function opens the file, takes the first line and save it in Natoms
        implicit none
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
        read(*,*)
        do i=1,n
                read(*,*) coord(i,1), coord(i,2), coord(i,3),mass(i)
        end do
end subroutine read_molecule

!
END PROGRAM readfile
