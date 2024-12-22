# Installation and Execution Instructions

## Compiler
- Fortran compiler
- Shell environment

## Compilation
- Navigate to `src` directory
- First, you need to create the input file, by default Ar atoms
  Run the commands:
$ gfortran -o input_generator input.f90
$ ./input_generator
Follow the on-screen prompts to define the system (number of atoms,
simulation box size, and input file name)

- Second, you run the MD simulation for your system defined in the
input file. Run the commands:
$ gfortran -o md_simulation main.f90
$ ./md_simulation
Follow the on-screen prompts to define the number of steps and
the trajectory output file name. The simulation trajectory will
be saved in the output file name, output_file.xyz.
