program deposition
implicit none
integer :: N,M
integer :: tt, ii, jj, jp, jn, lX, rX, cX
integer,dimension(:), allocatable :: X
real :: p,r
N = 100
M = 1000000
p = 0.5
allocate(X(N))
do tt=1,M
	jj = mod(INT(RAND(0)*N), N) + 1
	jp = mod(jj-1,N) + 1
	jn = mod(jj+1,N) + 1		
	lX=X(jp)
	rX=X(jn)
	cX=X(jj)+1
	if((cX > lX).AND.(cX > rX)) then
		r = RAND(0)
		if(r.LE.p)then
			X(jp)=lX+1
		else
			X(jn)=rX+1
		end if	
	elseif(cX > lX) then
		X(jp)=lX+1
	elseif(cX > rX) then
		X(jn)=rX+1
	else
		X(jj)=cX
	end if
	print*,"Estimate is:",real(sum(X))/(N),"+/-", sqrt(real(sum((X-real(sum(X))/N)**2))/N)
end do
deallocate(X)
end program


