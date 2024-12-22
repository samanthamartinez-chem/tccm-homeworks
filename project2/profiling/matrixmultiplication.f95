!----------Declaration area-----------
PROGRAM matrixmultiplication
!
IMPLICIT NONE
!
double precision, allocatable :: rowA(:),colA(:),valA(:)
double precision, allocatable :: rowB(:),colB(:),valB(:)
double precision, allocatable :: rowR(:),colR(:),valR(:)
double precision, allocatable :: denseA(:,:),denseB(:,:),denseR(:,:),denseC(:,:)
double precision,allocatable::finalC(:,:)
character(len=100)::input_fileA,input_fileB
integer :: i,nA,mA,totvalA,nB,mB,totvalB,totvalR,nC,mC,nmult
real(8)::isparse,fsparse,tsparse,idense,fdense,tdense,fillA,fillB,fillR,iblas,fblas,tblas
double precision::alpha,beta
alpha=1.0d0
beta=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Execution area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The program does a matrix multiplication of a matrix A and matrix B in
!sparse format, dense format and with the BLAS function

!In this part we will define the matrix A:
!In the subroutine get_dimention we get from the input file the nA=number 
!of rows, mA=number of columns and totvalA= number of values different from
!zero
!In the subroutine read_matrix we use these values to save the values in 
!three 1D arrays: rowA= row index, colA=column index and valA=value
!Also, with the subroutine dense_matrix we write this matrix in dense format
!In addition, we print the filling degree of the matrix

input_fileA="MATRIX_125_1p"
call get_dimension(input_fileA,nA,mA,totvalA)
call read_matrix(input_fileA,totvalA,rowA,colA,valA)
call dense_matrix(nA,mA,rowA,colA,valA,totvalA,denseA)
fillA=real(totvalA)/(real(nA)*real(mA))
write(*,*) 'The filling degree of A is ',fillA

!We do the same to define the matrix B

input_fileB="MATRIX_125_1p"
call get_dimension(input_fileB,nB,mB,totvalB)
call read_matrix(input_fileB,totvalB,rowB,colB,valB)
call dense_matrix(nB,mB,rowB,colB,valB,totvalB,denseB)
fillB=real(totvalB)/(real(nB)*real(mB))
write(*,*) 'The filling degree of B is ',fillB

!Now that we have the matrixA and matrixB in sparse and dense format we can
!multiply them. First we check that the two matrixes can be multiplied by
!checking that the mA(number of columns of A) is the same than nB(number of
!rows of B). Then we multiply in sparse format with the subroutine
!multiply_sparse. Then, we multiply them in dense format with the subroutine
!multiply_dense and finally we multiply them in dense format using the
!function BLAS. Also, the time of each multiplication is recorded and printed
!The filling degree of the resulting multiplied matrix is printed
!Also, for the sparse multiplication, the number of multiplications in 
!comparisson to the theoretical maximum is printed

if (mA==nB) then
        call cpu_time(isparse)
        call multiply_sparse(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR,totvalR,nmult)
        call cpu_time(fsparse)
        tsparse=fsparse-isparse
        fillR=real(totvalR)/(real(nA)*real(mB))
        write(*,*) 'The filling degree of the resulting matrix is ',fillR
        write(*,*) 'The number of multiplications in sparse format is ',nmult,' in comparisson to the theoretical maximum ',nA**3
        write(*,*) 'The time taken for the sparse multiplication is ',tsparse,' seconds'
        call cpu_time(idense)
        call multiply_dense(nA,mA,denseA,mB,denseB,nC,mC,denseC)
        call cpu_time(fdense)
        tdense=fdense-idense
        write(*,*) 'The time taken for the dense multiplication is ',tdense,' seconds'
        allocate(finalC(nA,mB))
        call cpu_time(iblas)
        call DGEMM('N','N',nA,mB,mA,alpha,denseA,nA,denseB,nB,beta,finalC,nA)
        call cpu_time(fblas)
        tblas=fblas-iblas
        write(*,*) 'The time taken for the multiplication using BLAS is ',tblas,' seconds'
else
        write(*,*) 'Multiplication not possible'
end if

!All of the matrixes used in this program are allocated, so the final step is
!to deallocate them

deallocate(rowA,colA,valA,denseA)
deallocate(rowB,colB,valB,denseB)
deallocate(rowR,colR,valR)
deallocate(denseC,finalC)

contains

!This subroutine gets the dimension of the matrix in sparse format. It opens the input file and
!then it reads the dimension that is in the 5th line of the file and saves it in the variables
!n(number of rows) and m(number of columns)

subroutine get_dimension(input_file,n,m,totval)
        implicit none
        character(len=100),intent(in)::input_file
        character(len=100)::dummy
        integer ::iostat
        integer,intent(out) :: n,m,totval
        totval=0
        open(unit=10, file=input_file, status="old", action="read")
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) dummy,dummy,dummy,dummy,n,dummy,m
        do 
                read(10,'(A)',iostat=iostat)
                if (iostat/=0)exit 
                totval=totval+1
        end do
        close(unit=10)
end subroutine get_dimension


!This subroutine saves the values of the input file in three 1D arrays.
!It starts saving the data from the 6th line of the file and saves the
!first column as row(rows), the second as col(columns) and the third as
!val(the values of the matrix).

subroutine read_matrix(input_file,totval,row,col,val)
        implicit none
        character(len=100),intent(in)::input_file
        double precision,intent(out),allocatable::row(:),col(:),val(:)
        integer,intent(out) :: totval
        double precision::temp
        integer :: i,i_stat

        allocate(row(totval),stat=i_stat)
        allocate(col(totval),stat=i_stat)
        allocate(val(totval),stat=i_stat)
        if (i_stat /= 0) then
                print *, "Memory allocation failed"
                stop
        end if
        open(unit=10, file=input_file, status="old", action="read")       
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        do i=1,totval
                read(10,*) row(i), col(i), val(i)
        end do
        close(10)
        !do i=1,totval
        !        if (col(i)<row(i)) then
        !                write(*,*) 'yes'
        !                temp=row(i)
        !                row(i)=col(i)
        !                col(i)=temp
        !        end if
        !end do        
        end subroutine read_matrix
        

!This subroutine multiplies the matrix in sparse format. For each i row index of matrix A it looks for
!any row index of matrix B that is the same as the i column index of matrix A. Then, the values of 
!matrix A and B that correspond to those indices are multiplied and saved in the same row index as 
!the row index of matrix A and the same column index as the column index of matrix B. These values are
!saved in three 1D arrays, rowR, colR and valR (row, column and values). 
!Given that it is possible that more than one multiplication contributes to a certain row/column index,
!there is also a loop that checks if there is already a value saved for a certain row/column value and
!sum the multiplication to this value.
!Finally, the total number of multiplications and the number of values different from zero are saved
!in variables (nmult and totvalR respectively).

subroutine multiply_sparse(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR,totvalR,nmult)
        implicit none
        double precision,intent(in),allocatable::rowA(:),colA(:),valA(:)
        double precision,intent(in),allocatable::rowB(:),colB(:),valB(:)
        double precision,intent(out),allocatable::rowR(:),colR(:),valR(:)
        integer :: n,i,j,k,i_stat
        integer,intent(out)::totvalR,nmult
        integer,intent(in)::totvalA,totvalB
        logical::same
        n=1
        nmult=0
        allocate(rowR((totvalA+totvalB)*2))
        allocate(colR((totvalA+totvalB)*2))
        allocate(valR((totvalA+totvalB)*2))
        do i=1,totvalA
                do j=1,totvalB
                        if (colA(i)==rowB(j)) then
                                same=.false.
                                if (n>1) then
                                        do k=1,n
                                                if (rowR(k)==rowA(i) .and. colR(k)==colB(j)) then
                                                        valR(k)=valR(k)+valA(i)*valB(j)
                                                        nmult=nmult+1
                                                        same=.true.
                                                end if
                                        end do
                                        if (.not. same) then
                                                rowR(n)=rowA(i)
                                                colR(n)=colB(j)
                                                valR(n)=valA(i)*valB(j)
                                                n=n+1
                                                nmult=nmult+1
                                        end if
                                else
                                        rowR(n)=rowA(i)
                                        colR(n)=colB(j)
                                        valR(n)=valA(i)*valB(j)
                                        n=n+1
                                        nmult=nmult+1
                                end if       
                        end if
                end do
         end do
         totvalR=n-1
         end subroutine multiply_sparse

!This subroutine converts the sparse matrix into a dense format. It creates a 
!matrix with the desired dimension and put every value as zero. Then it 
!checks for the indexes of rows and columns that exist in the 1D vectors and
!saves it values.

subroutine dense_matrix(n,m,row,col,val,totval,dense)
        implicit none
        double precision,intent(in),allocatable::row(:),col(:),val(:)
        integer,intent(in)::n,m,totval
        integer::i,j,k
        double precision,intent(out),allocatable::dense(:,:)
        allocate(dense(n,m))
        do i=1,n
                do j=1,m
                        dense(i,j)=0.000000000
                        do k=1,totval 
                                if (i==row(k) .and. j==col(k)) then
                                        dense(i,j)=val(k)
                                        exit
                                end if
                        end do
                end do
        end do
        end subroutine dense_matrix        

!This subroutine multiply the matrices in dense format. It reads the matrix
!in dense format and multiply the matrix A rows with the matrix B columns.

subroutine multiply_dense(nA,mA,denseA,mB,denseB,nC,mC,denseC)
        implicit none        
        double precision,intent(in),allocatable::denseA(:,:),denseB(:,:)
        double precision,intent(out),allocatable::denseC(:,:)
        integer,intent(in)::nA,mA,mB
        integer,intent(out)::nC,mC
        integer::i,j,k
        nC=nA
        mC=mB
        allocate(denseC(nC,mC))
        do i=1,nC
                do j=1,mC
                        denseC(i,j)=0
                        do k=1,mA
                                denseC(i,j)=denseC(i,j)+denseA(i,k)*denseB(k,j)
                        end do
                end do
        end do
        end subroutine multiply_dense

END PROGRAM matrixmultiplication
