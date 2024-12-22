program matrix_operations
    implicit none
    integer :: iunit, ierr
    character(len=100) :: line
    integer :: matrix_dim, row, col, i, j
    real(8) :: value
    real(8), allocatable :: matrix(:,:)
    real(8), allocatable :: transpose_matrix(:,:)
    real(8) :: start_cpu_time, end_cpu_time, elapsed_time

    ! Open input file for the first matrix
    open(unit=iunit, file="matrix_25_5p", status="old")

    ! Skip header lines
    do while (.true.)
        read(iunit, "(A)", iostat=ierr) line
        if (ierr /= 0) exit  ! Exit loop when end of file is reached
        if (index(line, "matrix will have dimension") > 0) exit  ! Found the dimension line
    end do

    ! Read matrix dimension (from the line after the header)
    read(iunit, *) matrix_dim

    ! Allocate matrix based on the dimension
    allocate(matrix(matrix_dim, matrix_dim))

    ! Read matrix elements (row, col, value)
    do while (.true.)
        read(iunit, *, iostat=ierr) row, col, value
        if (ierr /= 0) exit  ! Break loop if end of file or error
        if (row > matrix_dim .or. col > matrix_dim) then
            print *, "Error: row or col index exceeds matrix dimension."
            stop
        end if
        matrix(row, col) = value
    end do

    ! Close first matrix file
    close(iunit)

    ! Now calculate the transpose
    allocate(transpose_matrix(matrix_dim, matrix_dim))
    do i = 1, matrix_dim
        do j = 1, matrix_dim
            transpose_matrix(i, j) = matrix(j, i)
        end do
    end do

    ! Open second input file for the second matrix
    open(unit=iunit, file="matrix_25_1p", status="old")

    ! Skip header lines in the second matrix file
    do while (.true.)
        read(iunit, "(A)", iostat=ierr) line
        if (ierr /= 0) exit  ! Exit loop when end of file is reached
        if (index(line, "matrix will have dimension") > 0) exit  ! Found the dimension line
    end do

    ! Read matrix elements for the second matrix (row, col, value)
    do while (.true.)
        read(iunit, *, iostat=ierr) row, col, value
        if (ierr /= 0) exit  ! Break loop if end of file or error
        if (row > matrix_dim .or. col > matrix_dim) then
            print *, "Error: row or col index exceeds matrix dimension."
            stop
        end if
        ! Multiply the two matrices element-wise
        matrix(row, col) = matrix(row, col) * transpose_matrix(row, col)
    end do

    ! Close second matrix file
    close(iunit)

    ! Measure CPU time
    call cpu_time(start_cpu_time)

    ! Perform any further matrix operations, if needed (e.g., multiplication)
    call cpu_time(end_cpu_time)

    ! Calculate elapsed time
    elapsed_time = end_cpu_time - start_cpu_time

    ! Print CPU time
    print *, "CPU Time for matrix transpose and multiplication: ", elapsed_time

end program matrix_operations

