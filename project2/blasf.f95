!----------Declaration area-----------
PROGRAM matrixmultiplication
!
IMPLICIT NONE
!
double precision, allocatable :: rowA(:),colA(:),valA(:)
double precision, allocatable :: rowB(:),colB(:),valB(:)
double precision, allocatable :: rowR(:),colR(:),valR(:)
double precision, allocatable :: denseA(:,:),denseB(:,:),denseR(:,:),denseC(:,:)
character(len=100)::input_fileA,input_fileB
integer :: i,nA,mA,totvalA,nB,mB,totvalB,totvalR,nC,mC
real(8)::isparse,fsparse,tsparse
double precision,allocatable::finalC(:,:)
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

input_fileA="matA"
call get_dimension(input_fileA,nA,mA,totvalA)
call read_matrix(input_fileA,totvalA,rowA,colA,valA)
call dense_matrix(nA,mA,rowA,colA,valA,totvalA,denseA)

!We do the same to define the matrix B

input_fileB="matB"
call get_dimension(input_fileB,nB,mB,totvalB)
call read_matrix(input_fileB,totvalB,rowB,colB,valB)
call dense_matrix(nB,mB,rowB,colB,valB,totvalB,denseB)

!Now that we have the matrixA and matrixB in sparse and dense format we can
!multiply them. First we check that the two matrixes can be multiply by
!checking that the mA(number of columns of A) is the same than nB(number of
!rows of B). Then we multiply in sparse format with the subroutine
!multiply_sparse. Then, we multiply them in dense format with the subroutine
!multiply_dense and finally we multiply them in dense format using the
!function BLAS
!Also, the time 

allocate(finalC(nA,mB))
call DGEMM('N','N',nA,mB,mA,alpha,denseA,nA,denseB,nB,beta,finalC,nA)

!write(*,*) tsparse

!All of the matrixes used in this program are allocated, so the final step is
!to deallocate them

!write(*,*) totvalR
!do i=1,totvalR
!        write(*,*) rowR(i),colR(i),valR(i)
!end do
!write(*,*) denseA,denseB
write(*,*) finalC

deallocate(rowA,colA,valA,denseA)
deallocate(rowB,colB,valB,denseB)
!deallocate(rowR,colR,valR)
!deallocate(denseC,finalC)
deallocate(finalC)

contains

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


subroutine read_matrix(input_file,totval,row,col,val)
        implicit none
        character(len=100),intent(in)::input_file
        double precision,intent(out),allocatable::row(:),col(:),val(:)
        integer,intent(out) :: totval
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
        !                row(i)=col(i)
        !                col(i)=row(i)
        !        end if
        !end do        
        end subroutine read_matrix
        
!

subroutine multiply_sparse(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR,totvalR)
        implicit none
        double precision,intent(in),allocatable::rowA(:),colA(:),valA(:)
        double precision,intent(in),allocatable::rowB(:),colB(:),valB(:)
        double precision,intent(out),allocatable::rowR(:),colR(:),valR(:)
        integer :: n,i,j,k,i_stat,totvalR
        integer,intent(in)::totvalA,totvalB
        logical::same
        n=1
        allocate(rowR(totvalA+totvalB))
        allocate(colR(totvalA+totvalB))
        allocate(valR(totvalA+totvalB))
        do i=1,totvalA
                do j=1,totvalB
                        if (colA(i)==rowB(j)) then
                                same=.false.
                                if (n>1) then
                                        do k=1,n
                                                if (rowR(k)==rowA(i) .and. colR(k)==colB(j)) then
                                                        valR(k)=valR(k)+valA(i)*valB(j)
                                                        same=.true.
                                                end if
                                        end do
                                        if (.not. same) then
                                                rowR(n)=rowA(i)
                                                colR(n)=colB(j)
                                                valR(n)=valA(i)*valB(j)
                                                n=n+1
                                        end if
                                else
                                        rowR(n)=rowA(i)
                                        colR(n)=colB(j)
                                        valR(n)=valA(i)*valB(j)
                                        n=n+1
                                end if       
                        end if
                end do
         end do
         totvalR=n-1
         end subroutine multiply_sparse
         
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
