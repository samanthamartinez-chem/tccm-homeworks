program matrix_operations
    implicit none
    integer :: iunit1, iunit2, i, j, k
    real(8) :: cpu_start, cpu_end
    real(8) :: random_fraction, scalefactor
    integer :: matrix_dim_1, matrix_dim_2, num_elements_1, num_elements_2
    integer, dimension(:), allocatable :: row_index_1, col_index_1, row_index_2, col_index_2
    real(8), dimension(:), allocatable :: matrix_values_1, matrix_values_2
    real(8), dimension(:,:), allocatable :: matrix_product
    character(len=255) :: line

    ! Open the input files for reading
    open(unit=iunit1, file="matrix_25_5p", status="old")
    open(unit=iunit2, file="matrix_25_1p", status="old")

    ! Read the first matrix
    call read_sparse_matrix(iunit1, matrix_dim_1, row_index_1, col_index_1, matrix_values_1, num_elements_1)

    ! Read the second matrix
    call read_sparse_matrix(iunit2, matrix_dim_2, row_index_2, col_index_2, matrix_values_2, num_elements_2)

    ! Check if matrix dimensions match for multiplication
    if (matrix_dim_1 /= matrix_dim_2) then
        print *, "Matrix dimensions do not match for multiplication"
        stop
    end if

    ! Allocate memory for the product matrix
    allocate(matrix_product(matrix_dim_1, matrix_dim_2))
    matrix_product = 0.0d0

    ! Start CPU time measurement
    call cpu_time(cpu_start)

    ! Perform matrix multiplication (sparse format)
    do i = 1, num_elements_1
        do j = 1, num_elements_2
            if (col_index_1(i) == row_index_2(j)) then
                matrix_product(row_index_1(i), col_index_2(j)) = matrix_product(row_index_1(i), col_index_2(j)) + matrix_values_1(i) * matrix_values_2(j)
            end if
        end do
    end do

    ! Stop CPU time measurement
    call cpu_time(cpu_end)

    ! Print the resulting matrix product (for checking purposes)
    print *, "Matrix multiplication result:"
    do i = 1, matrix_dim_1
        print *, matrix_product(i, :)
    end do

    ! Print the CPU time for the whole process
    print *, "CPU Time for matrix operations (s): ", cpu_end - cpu_start

    ! Close the files
    close(iunit1)
    close(iunit2)

contains

    ! Subroutine to read a sparse matrix from a file
    subroutine read_sparse_matrix(iunit, matrix_dim, row_index, col_index, matrix_values, num_elements)
        implicit none
        integer, intent(in) :: iunit
        integer, intent(out) :: matrix_dim
        integer, dimension(:), allocatable, intent(out) :: row_index, col_index
        real(8), dimension(:), allocatable, intent(out) :: matrix_values
        integer, intent(out) :: num_elements
        integer :: row, col, i
        real(8) :: value

        ! Read matrix dimension
        read(iunit,*) "matrix will have dimension", matrix_dim

        ! Initialize variables
        num_elements = 0

        ! First, count the number of non-zero elements
        do
            read(iunit, *, iostat=iunit) row, col, value
            if (iunit /= 0) exit
            num_elements = num_elements + 1
        end do

        ! Allocate arrays to store row, col indices and values
        allocate(row_index(num_elements), col_index(num_elements), matrix_values(num_elements))

        ! Rewind file to start reading the data again
        rewind(iunit)
        read(iunit,*)  ! Skip the first line (header)

        ! Populate arrays with non-zero elements
        do i = 1, num_elements
            read(iunit,*) row_index(i), col_index(i), matrix_values(i)
        end do
    end subroutine read_sparse_matrix

end program matrix_operations

