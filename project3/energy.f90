program energy_test
implicit none
integer, parameter :: Natoms = 5
double precision, allocatable :: velocity(:,:)
double precision :: mass(Natoms)
double precision :: kinetic_energy, potential_energy, total_energy
integer :: i

call allocate_velocities(Natoms, velocity)

do i = 1, Natoms
    mass(i) = 1.0d0  ! Todas las masas son 1.0 como ejemplo
    velocity(i, 1) = 0.1d0 * i  ! Velocidad en x
    velocity(i, 2) = 0.2d0 * i  ! Velocidad en y
    velocity(i, 3) = 0.3d0 * i  ! Velocidad en z
end do

potential_energy = 10.0d0 
kinetic_energy = compute_T(Natoms, velocity, mass)
total_energy = compute_E(kinetic_energy, potential_energy)

write(6,*) "The total kinetic energy (T) is: ", kinetic_energy
write(6,*) "The total potential energy (V) is: ", potential_energy
write(6,*) "The total energy (E = T + V) is: ", total_energy

deallocate(velocity)

contains

! Function to compute the T
double precision function compute_T(Natoms, velocity, mass)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: velocity(Natoms, 3)
    double precision, intent(in) :: mass(Natoms)
    integer :: i
    double precision :: vx, vy, vz
    ! Start count of T with 0
    compute_T = 0
    do i = 1, Natoms
        vx = velocity(i, 1)
        vy = velocity(i, 2)
        vz = velocity(i, 3)
        compute_T = compute_T + 0.5d0 * mass(i) * (vx**2 + vy**2 + vz**2)
     end do
end function compute_T

! Function to comput E = V + T

double precision function compute_E(T, potential_energy)
    implicit none
    double precision, intent(in) :: T, potential_energy
    compute_E = T + potential_energy
end function compute_E

! Subroutine to allocat velocities

subroutine allocate_velocities(Natoms, velocity)
    implicit none
    integer, intent(in) :: Natoms
    double precision, allocatable, intent(out) :: velocity(:,:)
    integer :: i_stat

    allocate(velocity(Natoms, 3), stat=i_stat)
        if (i_stat /= 0) then
            write(6,*) "Memory allocation for velocity failed"
            stop
        end if
    velocity = 0.0d0
end subroutine allocate_velocities

end program energy_test
