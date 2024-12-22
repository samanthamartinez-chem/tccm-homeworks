do step = 1, num_steps
    call position_update(Natoms, coord, velocity, acceleration, dt)
    call update_velocity_current(Natoms, velocity, acceleration, dt)
    call compute_distance_potential(Natoms, coord, distance, sigma, epsilon, potential_energy)
    call compute_acc(Natoms, coord, mass, distance, acceleration_new, sigma, epsilon)
    call update_velocity_final(Natoms, velocity, acceleration_new, dt)
    kinetic_energy = compute_T(Natoms, velocity, mass)
    total_energy = compute_E(kinetic_energy, potential_energy)
    if (mod(step, 1) == 0) then
        write(6, *) "STEP: ", step, "V: ", potential_energy, "T: ", kinetic_energy, "E: ", total_energy
    end if
    acceleration = acceleration_new
end do

