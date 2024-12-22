!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: rowA(:),colA(:),valA(:)
double precision, allocatable :: rowB(:),colB(:),valB(:)
double precision, allocatable :: rowR(:),colR(:),valR(:)
character(len=100)::input_fileA,input_fileB
integer :: i,nA,mA,totvalA,nB,mB,totvalB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Execution area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Define the name of the file to read
input_fileA="matA"

!Call the subroutine that gets the coordinates and the mass of the atoms
call get_dimension(input_fileA,nA,mA,totvalA)
call read_matrix(input_fileA,totvalA,rowA,colA,valA)

input_fileB="matB"
call get_dimension(input_fileB,nB,mB,totvalB)
call read_matrix(input_fileB,totvalB,rowB,colB,valB)

call multiply_matrix(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR)
do i=1,8
        write(*,*) rowR(i),colR(i),valR(i)
end do

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
        !        if (col(i)>row(i)) then
        !                row(i)=col(i)
        !                col(i)=row(i)
        !        end if
        !end do        
        end subroutine read_matrix
        
!

subroutine multiply_matrix(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR)
        implicit none
        double precision,intent(in),allocatable::rowA(:),colA(:),valA(:)
        double precision,intent(in),allocatable::rowB(:),colB(:),valB(:)
        double precision,intent(out),allocatable::rowR(:),colR(:),valR(:)
        integer :: n,i,j,k,i_stat
        integer,intent(in)::totvalA,totvalB
        logical::same
        write(*,*) totvalA
        n=1
        allocate(rowR(totvalA+totvalB))
        allocate(colR(totvalA+totvalB))
        allocate(valR(totvalA+totvalB))
        !rowR(n)=0.00000000
        !colR(n)=0.00000000
        !valR(n)=0.00000000
        !same=.false.
        do i=1,totvalA
                do j=1,totvalB
                        if (colA(i)==rowB(j)) then
                                write(*,*) 'matches',rowA(i),colA(i),rowB(j),colB(j)
                                same=.false.
                                if (n>1) then
                                        do k=1,n
                                                if (rowR(k)==rowA(i) .and. colR(k)==colB(j)) then
                                                        valR(k)=valR(k)+valA(i)*valB(j)
                                                        write(*,*) 'la suma es',valR(k)
                                                        same=.true.
                                                end if
                                        end do
                                        if (.not. same) then
                                                rowR(n)=rowA(i)
                                                colR(n)=colB(j)
                                                valR(n)=valA(i)*valB(j)
                                                write(*,*) n,'val',valR(n)
                                                n=n+1
                            !                    same=.false.
                                        end if
                                else
                                        rowR(n)=rowA(i)
                                        colR(n)=colB(j)
                                        valR(n)=valA(i)*valB(j)
                                        write(*,*) '1st val',valR(n)
                                        n=n+1
                                end if       
                        end if
                end do
         end do
         write(*,*) n
         end subroutine multiply_matrix
         
END PROGRAM readfile
