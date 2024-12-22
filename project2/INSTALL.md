#Installation and Excecution Instructions

##Compiler
- Fortran compiler
- Shell environment

##Preparation of the files
In order to compile the program you need the input file of the
2 matrixes that will be multiplied (located in the test directory).
It is noteworthy that the input files of the 125x125 matrices were
modified to have the same format as the 25x25 matrices (5 lines
were added in the begining of the file).
These files must be in the same directory as the script 
matrixmultiplication.f95 (located in the scr directory).
Once you have the script and the two input files you have to
write in the script the name of the two input files in the variables
input_file1 and input_file2.

##Compilation
Run the commands:
$ gfortran -o mult multiplicationmatrix.f95 -lopenblas
$ ./mult

##Results
The script prints:
-The fillinmultiplicationmatrix.f95 -lopenblasg degree of the 2 initial matrices
-The filling degree of the resulting matrix
-The number of multiplications that the sparse multiplication
executes
-The time taken for the sparse multiplication
-The time taken for the dense multiplication
-The time taken for the multiplication using DGEMM


