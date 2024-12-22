!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: rowA(:), colA(:), valA(:)
double precision, allocatable :: rowB(:), colB(:), valB(:)
double precision, allocatable :: rowR(:), colR(:), valR(:)
character(len=100) :: input_fileA, input_fileB
integer :: i, nA, mA, totvalA, nB, mB, totvalB
double precision :: start_time, end_time, cpu_time_taken

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Execution area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the name of the file to read
input_fileA = "matA"

! Call the subroutine that gets the coordinates and the mass of the atoms
call get_dimension(input_fileA, nA, mA, totvalA)
call read_matrix(input_fileA, totvalA, rowA, colA, valA)

input_fileB = "matB"
call get_dimension(input_fileB, nB, mB, totvalB)
call read_matrix(input_fileB, totvalB, rowB, colB, valB)

! Start measuring CPU time
call cpu_time(start_time)

! Perform matrix multiplication
call multiply_matrix(rowA, colA, valA, totvalA, rowB, colB, valB, totvalB, rowR, colR, valR)

! End measuring CPU time
call cpu_time(end_time)

! Calculate and print the CPU time taken for multiplication
cpu_time_taken = end_time - start_time
print *, "CPU time for matrix multiplication: ", cpu_time_taken, " seconds"

! Output the results
do i = 1, min(8, size(rowR))
    write(*,*) rowR(i), colR(i), valR(i)
end do

! Deallocate arrays
deallocate(rowA)
deallocate(colA)
deallocate(valA)
deallocate(rowB)
deallocate(colB)
deallocate(valB)
deallocate(rowR)
deallocate(colR)
deallocate(valR)

contains

subroutine get_dimension(input_file, n, m, totval)
    implicit none
    character(len=100), intent(in) :: input_file
    integer, intent(out) :: n, m, totval
    integer :: i, iostat
    character(len=100) :: dummy

    totval = 0
    open(unit=10, file=input_file, status="old", action="read")
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *) dummy, dummy, dummy, dummy, n, dummy, m
    do
        read(10, *, iostat=iostat)
        if (iostat /= 0) exit
        totval = totval + 1
    end do
    close(unit=10)
end subroutine get_dimension

subroutine read_matrix(input_file, totval, row, col, val)
    implicit none
    character(len=100), intent(in) :: input_file
    double precision, intent(out), allocatable :: row(:), col(:), val(:)
    integer, intent(out) :: totval
    integer :: i, i_stat

    allocate(row(totval), stat=i_stat)
    allocate(col(totval), stat=i_stat)
    allocate(val(totval), stat=i_stat)
    if (i_stat /= 0) then
        print *, "Memory allocation failed"
        stop
    end if

    open(unit=10, file=input_file, status="old", action="read")
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    read(10, *)
    do i = 1, totval
        read(10, *) row(i), col(i), val(i)
    end do
    close(10)
end subroutine read_matrix

subroutine multiply_matrix(rowA, colA, valA, totvalA, rowB, colB, valB, totvalB, rowR, colR, valR)
    implicit none
    double precision, intent(in), allocatable :: rowA(:), colA(:), valA(:)
    double precision, intent(in), allocatable :: rowB(:), colB(:), valB(:)
    double precision, intent(out), allocatable :: rowR(:), colR(:), valR(:)
    integer, intent(in) :: totvalA, totvalB
    integer :: i, j, k, n
    logical :: same

    ! Allocate memory for the result matrix
    allocate(rowR(totvalA + totvalB))
    allocate(colR(totvalA + totvalB))
    allocate(valR(totvalA + totvalB))

    n = 1
    do i = 1, totvalA
        do j = 1, totvalB
            if (colA(i) == rowB(j)) then
                same = .false.
                ! Check if the entry already exists in the result matrix
                do k = 1, n
                    if (rowR(k) == rowA(i) .and. colR(k) == colB(j)) then
                        valR(k) = valR(k) + valA(i) * valB(j)
                        same = .true.
                    end if
                end do
                ! If not found, add new entry
                if (.not. same) then
                    rowR(n) = rowA(i)
                    colR(n) = colB(j)
                    valR(n) = valA(i) * valB(j)
                    n = n + 1
                end if
            end if
        end do
    end do
end subroutine multiply_matrix

END PROGRAM readfile

