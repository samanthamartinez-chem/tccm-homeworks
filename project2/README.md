############################################################################
######################## DIRECTORY STRUCTURE ###############################
############################################################################

This directory contains the following:
-A LICENSE file
-An AUTHORS file
-An INSTALL.md where the instructions to compile the program are described
-A tests directory containing some of the outputs obtained with our script
-A data directory with the matrices used in the project
-A doc directory containing the report (in pdf) where the computational
efficiency of the multiplication methods is discussed.
-A scr directory that contains the source file of our program

############################################################################
############################# PROJECT 2 ####################################
############################################################################

Our program is coded in fortran and consists on several subroutines:
1)get_dimension: gets the sparse matrix dimensions and number of values from
the input file
2)read_matrix: saves the input values in three 1D arrays containing the rows,
columns and values of the sparse matrix
3)dense_matrix: saves the sparse matrix in dense format
4)multiply_sparse: does the multiplication in sparse format of 2 matrices and
save the result in three 1D arrays (rows, columns and values)
5)multiply_dense: multiplies the 2 matrices in the conventional way

In our program, first 2 matrices are read from the sparse format, converted
to dense format and then multiplied in sparse and dense format. Also, they are
multiplied using the BLAS routine. The execution time is monitored for each
multiplication and the filling degree of each matrix is also calculated.
More details of the code are written as comments in the script.
Finally, the profiling of each type of multiplication was analyzed using
gprof.

Our program was used to multiply five 25x25 and six 125x125 matrices. 
In total 61 multiplications were done and the results were used
to discuss the computational efficiency of each multiplication method.
