subroutine herm_diag(matr,nn,eigvs,job,info)
	complex*8, dimension(:,:), intent(inout) :: matr
	integer*4,intent(out) :: info
	integer*4,intent(in) :: nn
	character*1, intent(inout) :: job
	integer*4 :: lwork
	real*4, dimension(:),allocatable,intent(inout) :: eigvs
	complex*8, dimension(:), allocatable :: work
	real*4, dimension(:), allocatable :: rwork
	! optimal lwork
	lwork=-1
	allocate(work(1))
	allocate(rwork(max(1, 3*nn-2)))
	if(.not.allocated(eigvs))allocate(eigvs(nn))
	if((job /= 'N').OR.(job /= 'V'))then
		job = 'N'
	end if
	call cheev('V','U',nn,matr,nn,eigvs,work,lwork,rwork,info)
	lwork = int(real(work(1)))
	deallocate(work)
	deallocate(rwork)


	! actual diag
	allocate(work(max(1,lwork)))
	allocate(rwork(max(1, 3*nn-2)))
	call cheev('V','U',nn,matr,nn,eigvs,work,lwork,rwork,info)

	deallocate(work)
	deallocate(rwork)
end subroutine
program hubbard
implicit none
integer*4 :: L, ii,jj, N, next, prev, x
real*4 :: hopt
real*4, dimension(:,:), allocatable :: H
L=10
N=2**L
allocate(H(L,L))
H=0
! loop through the possible states
do ii=1,N
	do jj=1,N
		
	end do
end do



end program
