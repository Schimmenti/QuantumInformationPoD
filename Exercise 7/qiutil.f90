module qiutil
implicit none
contains
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

subroutine write_vec_col(x,sz,fileUnit)
	real*4, dimension(:) :: x
	integer*4 :: sz, fileUnit
	integer*4 :: ii
	do ii=1,sz
		write(fileUnit,"(G0,A)")x(ii)
	end do
end subroutine

subroutine write_complex_vec_col(x,sz,fileUnit)
	complex*8, dimension(:) :: x
	integer*4 :: sz, fileUnit
	integer*4 :: ii
	do ii=1,sz
		write(fileUnit,"(G0,A,G0)")real(x(ii)),',',aimag(x(ii))
	end do
end subroutine

subroutine wrm_r4(x, nRows, nCols, fileUnit)
	real*4, dimension(:,:), intent(in) :: x
	integer*4, intent(in) :: nRows, nCols, fileUnit
	integer*4 :: ii,jj
	do ii=1, nRows
		do jj=1,nCols-1
			write(fileUnit,"(G0,A)", advance="no")x(ii,jj),','
		end do
		write(fileUnit, "(G0)")x(ii,nCols)
	end do
end subroutine

subroutine wrm_c8(x, nRows, nCols, fileUnit)
	complex*8, dimension(:,:), intent(in) :: x
	integer*4, intent(in) :: nRows, nCols, fileUnit
	integer*4 :: ii,jj
	do ii=1, nRows
		do jj=1,nCols-1
			write(fileUnit,"(A,G0,A,G0,A,A)", advance="no")'(',real(x(ii,jj)),',',aimag(x(ii,jj)),')','	'
		end do
		write(fileUnit,"(A,G0,A,G0,A)")'(',real(x(ii,nCols)),',',aimag(x(ii,nCols)),')'
	end do
end subroutine


function shift_difference(x,sz)result(s)
 	real*4, dimension(:) :: x
	integer*4 :: sz, shift

	real*4, dimension(:), allocatable :: s
	integer*4 :: ii

	allocate(s(sz))
	s(1)=0
	do ii=2,sz
		s(ii)=x(ii)-x(ii-1)
	end do
end function

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

character(len=20) function strr(num)
!   "Convert an integer to string."
    real, intent(in) :: num
    write (strr, "(G0)") num
    strr = adjustl(strr)
end function strr

end module
