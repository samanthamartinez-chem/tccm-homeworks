program test_acceleration

implicit none

contains

subroutine compute_acc(Natoms, coord, mass, distance, acceleration, sigma, epsilon)
! call compute_acc(Natoms, coord, mass, distance, acceleration, sigma, epsilon)
implicit none
integer, intent(in) :: Natoms
double precision, intent(in) :: coord(Natoms, 3)
double precision, intent(in) :: mass(Natoms)
double precision, intent(in) :: distance(Natoms, Natoms)
double precision, intent(out) :: acceleration(Natoms, 3)
double precision, intent(in) :: sigma, epsilon
integer :: i, j
double precision :: r, U_deriv

acceleration = 0.0d0
 
do i = 1, Natoms
    do j = 1, i-1
        r = distance(i, j)
        U_deriv = (24.0d0 * epsilon / r) * (2.0d0 * (sigma / r)**12 - (sigma/r) **6)
        acceleration(i, 1) = acceleration(i, 1) - (1.0d0 / mass(i)) * U_deriv * (coord(i, 1) - coord(j, 1)) / r
        acceleration(i, 2) = acceleration(i, 2) - (1.0d0 / mass(i)) * U_deriv * (coord(i, 2) - coord(j, 2)) / r
        acceleration(i, 3) = acceleration(i, 3) - (1.0d0 / mass(i)) * U_deriv * (coord(i, 3) - coord(j, 3)) / r
    
        acceleration(j, 1) = acceleration(j, 1) + (1.0d0 / mass(j)) * U_deriv * (coord(i, 1) - coord(j, 1)) / r
        acceleration(j, 2) = acceleration(j, 2) + (1.0d0 / mass(i)) * U_deriv * (coord(i, 2) - coord(j, 2)) / r
        acceleration(j, 3) = acceleration(j, 3) + (1.0d0 / mass(i)) * U_deriv * (coord(i, 3) - coord(j, 3)) / r
    end do
end do

end subroutine compute_acc

subroutine allocate_acceleration(Natoms, acceleration)
    implicit none
    integer, intent(in) :: Natoms
    double precision, allocatable, intent(out) :: acceleration(:,:)
    integer :: i_stat

    allocate(acceleration(Natoms, 3), stat=i_stat)
        if (i_stat /= 0) then
            write(6,*) "Memory allocation for acceleration failed"
            stop
        end if
end subroutine allocate_acceleration

end program
