module heutil
implicit none
contains
function rnd_hem(mSize)result(hem)
	integer*4 :: mSize
	integer*4 :: ii,jj, sz
	complex*8, dimension(:,:), allocatable :: hem
	
	if(mSize <= 0)then
		sz = 0
	else
		sz = mSize
	end if
	allocate(hem(sz,sz))
	do jj=1,sz
		do ii=1,jj-1
		hem(ii,jj)=complex(RAND(0),RAND(0))
		hem(jj,ii)=conjg(hem(ii,jj))
		end do
		hem(ii,ii)=complex(RAND(0),0)
	end do
end function



end module

program eigenproblem
use heutil
integer*4 ii,nn, stat,info,lwork
character(10) :: nnchar
complex*8, dimension(:,:), allocatable :: hem
complex*8, dimension(:), allocatable :: work
real*4, dimension(:), allocatable :: evs,rwork
if(COMMAND_ARGUMENT_COUNT() < 1)then
	nn = 100
else
! get the only command line argument given
	read(nnchar,FMT=*, IOSTAT=stat)nn
! if stat is different from zero we fall back into invalid state and we use default matrix size
	if((stat.NE.0).OR.(nn <= 0))then
		nn = 100
	end if
	
endif

lwork=-1
allocate(evs(nn))
allocate(work(max(1,lwork)))
allocate(rwork(max(1, 3*nn-2)))
hem = rnd_hem(nn)
call cheev('N','U',nn,hem,nn,evs,work,lwork,rwork,info)
print*,"Procedure finished with:",info
print*,"First 10 Eigenvalues:",evs
end program
