module heutil
implicit none
contains
! create a pair of standard gaussian distributed real
function box_muller_r4(u,v)result(z)
	real*4 :: u,v
	real*4, dimension(2) :: z
	real*4 :: rad,c,s
	rad=sqrt(-2*log(u))
	c=cos(2*acos(-1.0)*v)
	s=sin(2*acos(-1.0)*v)
	z(1)=rad*c
	z(2)=rad*s
end function
! create a complex with real and imaginary part normally distributed
function box_muller_c8(u,v)result(c)
	real*4 :: u,v
	real*4, dimension(2) :: z
	complex*8 :: c
	z = box_muller_r4(u,v)
	c= complex(z(1),z(2))
end function
! creates a random hermitian matrix with normal entries
function rnd_hem(mSize)result(hem)
	integer*4 :: mSize
	integer*4 :: ii,jj, sz
	complex*8, dimension(:,:), allocatable :: hem
	real*4, dimension(2) :: temp
	if(mSize <= 0)then
		sz = 0
	else
		sz = mSize
	end if
	allocate(hem(sz,sz))
	do jj=1,sz
		temp = box_muller_r4(RAND(0),RAND(0))
		hem(jj,jj)=complex(temp(1),0)
		do ii=1,jj-1
			hem(ii,jj)=box_muller_c8(RAND(0),RAND(0))
			hem(jj,ii)=conjg(hem(ii,jj))
		end do
	end do
end function
! creates a random diagonal hermitian matrix with normal entries
function rnd_hem_diag(mSize)result(hem)
	integer*4 :: mSize
	integer*4 :: ii,jj, sz
	complex*8, dimension(:,:), allocatable :: hem
	real*4, dimension(2) :: temp
	if(mSize <= 0)then
		sz = 0
	else
		sz = mSize
	end if
	allocate(hem(sz,sz))
	hem = complex(0.0,0.0)
	do jj=1,sz
		temp = box_muller_r4(RAND(0),RAND(0))
		hem(jj,jj)=complex(temp(1),0)
	end do
end function

! froma sorted array creates an histogram
function sorted_hist(oArray,aSize,nBins,left,right,h)result(counts)
	real*4, dimension(:) :: oArray
	real*4, intent(in) :: left,right
	integer*4 :: aSize,nBins,ii,jj,total
	real*4 :: delta,temp,limit
	integer*4, dimension(:),allocatable :: counts
	real*4, intent(out) :: h
	logical :: flag
	flag = .TRUE.
	delta = oArray(aSize)-oArray(1)
	h = (right-left)/nBins
	allocate(counts(nBins))
	counts = 0
	limit = left+h
	jj=1
	total=0
	do ii=1,nBins
		flag = .TRUE.
		do while((jj <= aSize).AND.flag)
			temp = oArray(jj)
			if(temp > limit)then
				limit=limit+h
				flag=.FALSE.
			else
				counts(ii)=counts(ii)+1
				total=total+1
				jj=jj+1
			end if
		end do
	end do
end function

subroutine sort_merge(arr,left,center,right)
	real*4, dimension(:), intent(inout) :: arr
	real*4, dimension(:), allocatable :: b
	integer*4, intent(in) :: left,center,right
	real*4 :: temp
	integer*4 :: ii,jj,kk
	
	if(right-left==1)then
		if(arr(left)>arr(right))then
			temp = arr(left)
			arr(left)=arr(right)
			arr(right)=temp
		end if
	else
		allocate(b(right-left+1))
		ii = left
		jj = center+1
		kk = 1
		do while((ii <= center).AND.(jj <= right))
			if(arr(ii)<=arr(jj))then
				b(kk)=arr(ii)
				ii=ii+1
			else
				b(kk)=arr(jj)
				jj=jj+1
			end if
			kk=kk+1
		end do
		do while(ii <= center)
			b(kk)=arr(ii)
			ii=ii+1
			kk=kk+1
		end do
		do while(jj <= right)
			b(kk)=arr(jj)
			jj=jj+1
			kk=kk+1
		end do
		do kk=left,right
			arr(kk)=b(kk-left+1)
		end do
		deallocate(b)
	end if
end subroutine

recursive subroutine sort_split(arr,left,right)
	real*4, dimension(:), intent(inout) :: arr
	integer*4, intent(in) :: left, right
	integer*4 :: center
	if(left < right)then
		center=(left+right)/2
		call sort_split(arr,left,center)
		call sort_split(arr,center+1,right)
		call sort_merge(arr,left,center,right)
	end if
end subroutine

! in place array ordering
subroutine mergesort(arr,sz)
	real*4, dimension(:), intent(inout) :: arr
	integer*4, intent(in) :: sz
	call sort_split(arr,1,sz)
end subroutine

! hermitian matrix diagonalization
subroutine hediagonalize(matr,nn,eigvs,info)
	complex*8, dimension(:,:), intent(inout) :: matr
	integer*4,intent(out) :: info
	integer*4,intent(in) :: nn
	integer*4 :: lwork
	real*4, dimension(:),allocatable,intent(inout) :: eigvs
	complex*8, dimension(:), allocatable :: work
	real*4, dimension(:), allocatable :: rwork
	! optimal lwork
	lwork=-1
	allocate(work(1))
	allocate(rwork(max(1, 3*nn-2)))
	if(.not.allocated(eigvs))allocate(eigvs(nn))
	call cheev('N','U',nn,matr,nn,eigvs,work,lwork,rwork,info)
	lwork = int(real(work(1)))
	deallocate(work)
	deallocate(rwork)


	! actual diag
	allocate(work(max(1,lwork)))
	allocate(rwork(max(1, 3*nn-2)))
	call cheev('N','U',nn,matr,nn,eigvs,work,lwork,rwork,info)

	deallocate(work)
	deallocate(rwork)
end subroutine

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

! computes a local mean
! the mean at position i is given by the mean from i-steps to i+steps
function localmean(x, sz,steps)result(xm)
	real*4, dimension(:) :: x
	integer*4 :: sz, steps
	real*4, dimension(:), allocatable :: xm
	integer*4 :: ii,jj,ll,rl
	allocate(xm(sz))
	do ii=1,sz
		ll = max(1,ii-steps)
		rl = min(sz,ii+steps)
		xm(ii)=sum(x(ll:rl))/(rl-ll+1)
	end do
end function

function nn_minmax_ratio(x, sz)result(r)
	real*4, dimension(:),intent(in) :: x
	integer*4, intent(in) :: sz
	integer*4 :: ii
	real*4, dimension(:), allocatable :: r
	allocate(r(sz-1))
	do ii=2,sz
		r(ii-1)=min(x(ii), x(ii-1))/max(x(ii),x(ii-1))
	end do
end function

end module

program eigenproblem
use heutil
integer*4 ii,rr,nn,reps,nBins, stat,info,mode, avgMode
character(10) :: nnchar
complex*8, dimension(:,:), allocatable :: hem
real*4, dimension(:), allocatable :: evs,sp,ratios, lAvgs
real*4 :: avgSp,binWidth,cutoff
integer*4, dimension(:), allocatable :: counts
character(:), allocatable :: fname, fname2
cutoff=3.5

! CMD ARGS ARE:
! MATRIX SIZE
! NUMBER OF MATRICES TO GENERATE
! NUMBER OF BINS
! MODE (HERMITIAN==0, HERMITIAN DIAGONAL)
! LOCAL==0, GLOBAL==1

if(COMMAND_ARGUMENT_COUNT()  < 1)then
	nn = 100
else
	call GET_COMMAND_ARGUMENT(1, nnchar)
	read(nnchar,FMT=*, IOSTAT=stat)nn
	if((stat.NE.0).OR.(nn <= 0))then
		nn = 100
	end if
	if(COMMAND_ARGUMENT_COUNT()>1)then
		call GET_COMMAND_ARGUMENT(2, nnchar)
		read(nnchar,FMT=*, IOSTAT=stat)reps
		if((stat.NE.0).OR.(reps <= 0))then
			reps = 1
		end if
	end if

	if(COMMAND_ARGUMENT_COUNT()>2)then
		call GET_COMMAND_ARGUMENT(3, nnchar)
		read(nnchar,FMT=*, IOSTAT=stat)nBins
		if((stat.NE.0).OR.(nBins <= 0))then
			nBins = 20
		end if
	end if	

	if(COMMAND_ARGUMENT_COUNT()>3)then
		call GET_COMMAND_ARGUMENT(4, nnchar)
		read(nnchar,FMT=*, IOSTAT=stat)mode
		if((stat.NE.0).OR.(mode < 0))then
			mode = 0
		end if
	end if	

	if(COMMAND_ARGUMENT_COUNT()>4)then
		call GET_COMMAND_ARGUMENT(5, nnchar)
		read(nnchar,FMT=*, IOSTAT=stat)avgMode
		if((stat.NE.0).OR.(avgMode < 0))then
			avgMode = 0
		end if
	end if	
endif



fname = "hh_"//trim(str(nn))//"_"//trim(str(reps))//"_"//trim(str(nBins))//"_"//trim(str(mode))//"_"//trim(str(avgMode))//".dat"
fname2 = "hr_"//trim(str(nn))//"_"//trim(str(reps))//"_"//trim(str(nBins))//"_"//trim(str(mode))//"_"//trim(str(avgMode))//".dat"



open(unit=42, file=fname,action="write",status="replace")
open(unit=420, file=fname2,action="write",status="replace")
rr=1
write(42,"(A,G0)")"#cutoff=",cutoff
do while(rr <=reps)
	write(*,"(F6.2,A)")(100.0*rr)/reps,"%"
	if(mode==0)then
		hem = rnd_hem(nn)
	else
		hem = rnd_hem_diag(nn)
	end if
	call hediagonalize(hem,nn,evs,info)
	
	if(info==0)then
!		Diagonalization completed.
		allocate(sp(nn-1))

		do ii=2,nn
			sp(ii-1)=evs(ii)-evs(ii-1)
		end do

!		Writes the ratios
		ratios = nn_minmax_ratio(sp, nn-1)
		lAvgs = localmean(sp, nn-1, int(nn/10.0))
		do ii=1,nn-3
			write(420,"(G0)",advance="no")ratios(ii)
			write(420,"(A)",advance="no")","
		end do
		write(420,"(G0)",advance="no")ratios(nn-2)
		write(420,*)



		if(avgMode==0)then
			avgSp = sum(sp)/(nn-1)
			sp = sp/avgSp
		else
			sp = sp/lAvgs
		end if

		
!		Writes the histogram
		call mergesort(sp,nn-1)
		counts = sorted_hist(sp, nn-1, nBins,0.0,cutoff,binWidth)
		write(42,"(G0)",advance="no")binWidth
		write(42,"(A)",advance="no")","
		write(42,"(I0)",advance="no")nBins
		write(42,"(A)",advance="no")","
		do ii=1,nBins-1
			write(42,"(I0)",advance="no")counts(ii)
			write(42,"(A)",advance="no")","
		end do
		write(42,"(I0)",advance="no")counts(nBins)
		write(42,*)
		print*,real(sum(counts))/(nn-1)

		




		
		deallocate(counts)
		deallocate(hem)
		deallocate(evs)
		deallocate(sp)
		rr=rr+1
	else
!		print*,"Error in diagonalization."
	end if

	
end do
close(42)
close(420)
end program
