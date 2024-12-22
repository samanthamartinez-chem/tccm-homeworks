!----------Declaration area-----------
PROGRAM readfile
!
IMPLICIT NONE
!
double precision, allocatable :: rowA(:),colA(:),valA(:)
double precision, allocatable :: rowB(:),colB(:),valB(:)
double precision, allocatable :: rowR(:),colR(:),valR(:)
double precision, allocatable :: denseA(:,:),denseB(:,:),denseR(:,:),denseC(:,:)
double precision,allocatable::finalC(:,:)
character(len=100)::input_fileA,input_fileB
integer :: i,nA,mA,totvalA,nB,mB,totvalB,totvalR,nC,mC
real::alpha,beta
alpha=1.0d0
beta=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Execution area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Define the name of the file to read
input_fileA="matA"

!Call the subroutine that gets the coordinates and the mass of the atoms
call get_dimension(input_fileA,nA,mA,totvalA)
call read_matrix(input_fileA,totvalA,rowA,colA,valA)
call dense_matrix(nA,mA,rowA,colA,valA,totvalA,denseA)
!write(*,*) 'A'
!do i=1,25
!        write(*,*) denseA(i,1),denseA(i,2),denseA(i,3),denseA(i,4)
!end do

input_fileB="matB"
call get_dimension(input_fileB,nB,mB,totvalB)
call read_matrix(input_fileB,totvalB,rowB,colB,valB)
call dense_matrix(nB,mB,rowB,colB,valB,totvalB,denseB)
!write(*,*) 'B'
!do i=1,25
!        write(*,*) denseB(i,1),denseB(i,2),denseB(i,3),denseB(i,4)
!end do

!if (mA==nB) then
!        call multiply_matrix(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR,totvalR)
    !    write(*,*) 'Sparsemultiplication'
     !   do i=1,totvalR
     !           write(*,*) rowR(i),colR(i),valR(i)
      !  end do 
!else
!        write(*,*) 'Multiplication not possible'
!end if

!call dense_matrix(mA,nB,rowR,colR,valR,totvalR,denseR)
!write(*,*) 'R'
!do i=1,25
!        write(*,*) denseR(i,1),denseR(i,2),denseR(i,3),denseR(i,4)
!end do

!call multiply_dense(nA,mA,denseA,mB,denseB,nC,mC,denseC)
!write(*,*) 'C'
!do i=1,25
!        write(*,*) denseC(i,1),denseC(i,2),denseC(i,3),denseC(i,4)
!end do

allocate(finalC(nA,mB))
call DGEMM('T','T',nA,mB,mA,alpha,denseA,nA,denseB,nB,beta,finalC,nA)
write(*,*) 'finalC'
do i=1,4
        write(*,*) finalC(i,1),finalC(i,2),finalC(i,3),finalC(i,4)
end do

write(*,*) finalC(1,1)

deallocate(rowA)
deallocate(colA)
deallocate(valA)
deallocate(rowB)
deallocate(colB)
deallocate(valB)
!deallocate(rowR)
!deallocate(colR)
!deallocate(valR)
!deallocate(denseA)
!deallocate(denseB)
!deallocate(denseR)
!deallocate(denseC)
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

subroutine multiply_matrix(rowA,colA,valA,totvalA,rowB,colB,valB,totvalB,rowR,colR,valR,totvalR)
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
         totvalR=n
         end subroutine multiply_matrix
         
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
        !write(*,*) dense(5,10),dense(8,13),dense(18,25)
        !do i=1,25
        !        write(*,*) dense(i,1),dense(i,2),dense(i,3),dense(i,4)
        !end do
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

END PROGRAM readfile
