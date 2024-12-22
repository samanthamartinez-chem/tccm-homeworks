program test_verlet

implicit none

contains

! subroutine to update coordintes using eq. 7
subroutine position_update(Natoms, coord, velocity, acceleration, dt)
implicit none
integer, intent(in) :: Natoms
double precision, intent(inout) :: coord(Natoms, 3)
double precision, intent(in) :: velocity(Natoms, 3)
double precision, intent(in) :: acceleration(Natoms, 3)
double precision, intent(in) :: dt
integer :: i, j

do i = 1, Natoms
    do j = 1, 3
        coord(i, j) = coord(i, j) + velocity(i, j) * dt + 0.5d0 * acceleration(i, j) * dt**2
    end do
end do
end subroutine position_update

! subroutine to update the velocity (which depends on the current acceleration, a_n)
subroutine update_velocity_current(Natoms, velocity, acceleration, dt)
implicit none
integer, intent(in) :: Natoms
double precision, intent(inout) :: velocity(Natoms, 3)
double precision, intent(in) :: acceleration(Natoms, 3)
double precision, intent(in) :: dt
integer :: i, j

do i = 1, Natoms
    do j = 1, 3
        velocity(i, j) = velocity(i, j) + 0.5d0 * acceleration(i, j) * dt
    end do
end do
end subroutine update_velocity_current

! subroutine to update the velocity (Finalize the update using the a_{n+1})
! takes as an input the acceleration computed from the updated positions
subroutine update_velocity_final(Natoms, velocity, acceleration_new, dt)
implicit none
integer, intent(in) :: Natoms
double precision, intent(inout) :: velocity(Natoms, 3)
double precision, intent(in) :: acceleration_new(Natoms, 3)
double precision, intent(in) :: dt
integer :: i, j

do i = 1, Natoms
    do j = 1, 3
        velocity(i, j) = velocity(i, j) + 0.5d0 * acceleration_new(i, j) * dt
    end do
end do

end subroutine update_velocity_final

end program test_verlet

