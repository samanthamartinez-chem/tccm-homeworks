# PedASMD: Molecular Dynamics Simulation Code

Welcome to PedASMD, a Fortran-based molecular dynamics simulation program designed to simulate particle dynamics using the Lennard-Jones potential in 3D space.

## Features
- Generates input files with random atomic positions in a cubic box.
- Computes forces, accelerations, and potential energies using the Lennard-Jones potential.
- Implements the Verlet integration algorithm for particle motion.
- Outputs the simulation trajectory in XYZ format.

## How to Use
For details on how to compile and run the program, please refer to `INSTALL.md`.

## Directories
- `src` contains the main code `main.f90` and the input generator code `input.f90`.
- `src/construction` contains different pieces of the main code, as we built it.
- `tests` contains an input example and a trajectory example obtained with our code.
- `inputs` contains the input generator code and some examples.
