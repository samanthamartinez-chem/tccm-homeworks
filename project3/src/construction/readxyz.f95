!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: A(:,:)
!character(len=1),dimension(100) :: A_atoms
character(len=2),allocatable :: A_atoms(:)
integer :: i_stat,i,m,n
!
!----------Execution area-----------
m=3
!Read matrix
read(*,*) n

allocate(A(m,n),stat=i_stat)
allocate(A_atoms(n),stat=i_stat)
if (i_stat /= 0) then
        print *, "Memory allocation failed!"
        stop
end if

read(*,*)
!Save the atoms in a vector and the coordinates in a matrix
do i=1,n
        read(*,*) A_atoms(i), A(1,i), A(2,i), A(3,i)
end do

write(*,*) A(1,1),A(2,1) 

deallocate(A)
deallocate(A_atoms)
!
END PROGRAM readfile
