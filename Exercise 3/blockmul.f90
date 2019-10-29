module matstuff
implicit none
contains
function RANDMATF(N)result(x)
	implicit none
	integer*4 :: N
	real*4, dimension(:,:), allocatable :: x
	integer*4 :: ii,jj,kk

	if(N <= 0)then
		allocate(x(0,0))
	else
		allocate(x(N,N))
		do jj=1,N
			do ii=1,N
				x(ii,jj)=RAND(0)
			end do
		end do
	end if
end function

function blockmulf(a,b,nn,bsz)result(c)
	implicit none
	integer*4  :: nn, bsz
	real*4, dimension(:,:) :: a, b
	real*4, dimension(:,:), allocatable :: c
	integer*4 :: bi,bj,bk,ii,jj,kk
	real*4 :: temp
	if(((nn>0).AND.(bsz>0)).AND.(mod(nn,bsz)==0))then
		allocate(c(nn,nn))
		c = 0E0
		do bj=1,nn,bsz
			do bk=1,nn,bsz
				do bi=1,nn,bsz
					do jj=0,bsz-1
						do kk=0,bsz-1
						temp = b(bk+kk,bj+jj)
							do ii=0,bsz-1
								c(bi+ii,bj+jj)=c(bi+ii,bj+jj)+a(bi+ii,bk+kk)*temp
							end do
						end do
					end do
				end do
			end do
		end do
	else
		allocate(c(0,0))
	end if
end function
end module

program blockmul
use matstuff
implicit none
integer*4 :: sz, bsz, stat
real*4, dimension(:,:), allocatable :: mat1, mat2, res
real*4 :: t0,t1
character(len=10) :: inChars
if(COMMAND_ARGUMENT_COUNT() == 2)then
	call GET_COMMAND_ARGUMENT(1, inChars)
	read(inChars,*,IOSTAT=stat)sz
	if(stat /= 0)then
		sz = 100
	end if
	call GET_COMMAND_ARGUMENT(2, inChars)
	read(inChars,*,IOSTAT=stat)bsz
	if(stat == 0)then
		if(mod(sz,bsz)/=0)then
			bsz = sz
		end if
	else
		bsz = sz
	end if
else
	sz = 100
	bsz = 100
end if

mat1 = RANDMATF(sz)
mat2 = RANDMATF(sz)

print*,sz,bsz

call cpu_time(t0)
res = blockmulf(mat1,mat2,sz,10)
call cpu_time(t1)
print*, t1-t0

call cpu_time(t0)
res = matmul(mat1,mat2)
call cpu_time(t1)
print*, t1-t0

deallocate(mat1)
deallocate(mat2)
deallocate(res)

end program
