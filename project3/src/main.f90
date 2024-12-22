!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! A Molecular Dynamics Code !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program molecular_dynamics

! **************************************************************************** !
! **************************** Declaration Area  ***************************** !
! **************************************************************************** !

implicit none

double precision, allocatable :: coord(:,:), velocity(:,:), acceleration(:,:)
double precision, allocatable :: acceleration_new(:,:), distance(:,:)
double precision, allocatable :: mass(:)
double precision :: potential_energy, kinetic_energy, total_energy
double precision, parameter :: epsilon = 0.0661d0 ! in j/mol
double precision, parameter :: sigma = 0.3345d0 ! in nm
double precision, parameter :: dt = 0.2d0
character(len=100)::input_file
integer :: Natoms, step, num_steps

! **************************************************************************** !
! ***************************** Execution Area  ****************************** !
! **************************************************************************** !

! We welcome the user
write(6, *) "*****************************************************"
write(6, *) "*            Welcome to PedASMD Simulation Code!    *"
write(6, *) "*   This program simulates particle dynamics using  *"
write(6, *) "*       the Lennard-Jones potential in 3D space.    *"
write(6, *) "*****************************************************"
write(6, *)
write(6, *) "Please follow the prompts to set up your simulation."
write(6, *)

! Ask for the name of the file to read

write(6, *) "Provide your input file name: "
read(5, *) input_file

! Ask the user to introduce the number of steps
write(6, *) "Introduce the number of steps of the simulation: "
read(5, *) num_steps

if (num_steps .le. 0) then
    write(6, *) "The number of steps must be larger than 0."
    stop
end if

!Call the function to determine the number of atoms
Natoms=read_Natoms(input_file)

!Call the subroutine that gets the coordinates and the mass of the atoms
call read_molecule(input_file,Natoms,coord,mass)

! Alocates memory for the distance matrix 
call allocate_distance(Natoms, distance)

! Computes the distance and LJ potential energy
call compute_distance_potential(Natoms, coord, distance, sigma, epsilon, potential_energy)

! Allocates velocity
call allocate_velocities(Natoms, velocity)

! Allocate accelerations
call allocate_acceleration(Natoms, acceleration)
call allocate_acceleration(Natoms, acceleration_new)

! Compute acceleration
call compute_acc(Natoms, coord, mass, distance, acceleration, sigma, epsilon)


open(unit=16, file='trajectory.xyz', status='replace', action='write')
do step = 1, num_steps
    call position_update(Natoms, coord, velocity, acceleration, dt)
    call update_velocity_current(Natoms, velocity, acceleration, dt)
    call compute_distance_potential(Natoms, coord, distance, sigma, epsilon, potential_energy)
    call compute_acc(Natoms, coord, mass, distance, acceleration_new, sigma, epsilon)
    call update_velocity_final(Natoms, velocity, acceleration_new, dt)
    kinetic_energy = compute_T(Natoms, velocity, mass)
    total_energy = compute_E(kinetic_energy, potential_energy)
    !if (mod(step, 1) == 0) then
       ! write(6, *) "STEP: ", step, "V: ", potential_energy, "T: ", kinetic_energy, "E: ", total_energy
    !end if
    if (mod(step, 10) == 0) then
        call write_xyz(step, Natoms, coord, kinetic_energy, potential_energy, total_energy)
    end if
    acceleration = acceleration_new
end do
close(16)

! For the user if simulation finished successfully
write(6, *) "Simulation completed successfully!"
write(6, *) "Trajectory file saved as 'trajectory.xyz'."

! Deallocate all we allocated 
deallocate(coord)
deallocate(mass)
deallocate(distance)
deallocate(velocity)
deallocate(acceleration)
deallocate(acceleration_new)


contains

! **************************************************************************** !
! *************************** Memory Allocation ****************************** !
! **************************************************************************** !

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



! **************************************************************************** !
! ****************************** Subroutines ********************************* !
! **************************************************************************** !

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


! Subroutine to create the format of the xyz file
subroutine write_xyz(step, Natoms, coord, kinetic_energy, potential_energy, total_energy)
implicit none
integer, intent(in) :: step, Natoms
double precision, intent(in) :: coord(Natoms, 3), kinetic_energy, potential_energy, total_energy
integer :: i

! Unit 10 bc later in the main we will open a file 
write(16,*) Natoms
write(16,*) "STEP: ", step, "T: ", kinetic_energy, "V: ", potential_energy, "E: ", total_energy

do i = 1, Natoms
    write(16, '(A,3F16.6)') "Ar", coord(i, 1), coord(i, 2), coord(i, 3)
end do

end subroutine write_xyz


! **************************************************************************** !
! ******************************* Functions ********************************** !
! **************************************************************************** !

integer function read_Natoms(input_file) result(Natoms)
        !This function opens the file, takes the first line and save it in Natoms
        implicit none
        character(len=*), intent(in)::input_file
        open(unit=10, file=input_file, status="old", action="read")
        read(10,*) Natoms
        close(unit=10)
end function read_Natoms


double precision function lj_potential(r, sigma, epsilon)
    implicit none
    double precision, intent(in) :: r, sigma, epsilon
    lj_potential = 4.0d0 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
end function lj_potential

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



end program molecular_dynamics
