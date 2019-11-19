module quantumutil
implicit none
contains
function laplmatr(N, pbc)result(A)
	integer*4 :: N
	logical :: pbc


	integer*4 :: ii
	complex*8, dimension(:,:), allocatable :: A

	allocate(A(N,N))
	A=0
	if(pbc)then
		A(1,1)=-2
		A(1,2)=1
		A(1,N)=1

		A(N,N)=-2
		A(N,N-1)=1
		A(N,1)=1	
	else
		A(1,1)=-2
		A(1,2)=1

		A(N,N)=-2
		A(N,N-1)=1
	end if
	do ii=2,N-1
		A(ii,ii)=-2
		A(ii,ii+1)=1
		A(ii,ii-1)=1
	end do
	return
end function

function hosc_potmatr(eps, L, m,omega)result(Vpot)
	! lattice size
	real*4 :: eps, m,omega, x_ii
! number of discretized points and loop index
	integer*4 :: L
	integer*4 :: N,ii
	complex*8, dimension(:,:), allocatable :: Vpot
	N=2*L+1
	allocate(Vpot(N,N))
	Vpot=0.0
	do ii=1,N
		x_ii=eps*(ii-L-1)
		Vpot(ii,ii)=(0.5*m)*((omega*x_ii)**2)
	end do
	return
end function

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

program shotproblem
use quantumutil
implicit none
integer*4 :: L,N,info,ii, argCount, stat
real*4 :: eps,hbar,m,omega,kinTyp,x0,mdiff
real*4, dimension(:), allocatable :: eig, res,diff
complex*8, dimension(:,:), allocatable :: H
character(10) :: argChar
character(:), allocatable :: fname, fname2
argCount=COMMAND_ARGUMENT_COUNT()

L=500
eps=1E-4
omega=1.0E4
m=1.0
hbar=1.0

if(argCount >= 1)then
	call GET_COMMAND_ARGUMENT(1, argChar)
	read(argChar,FMT=*, IOSTAT=stat)L
	if((stat.NE.0).OR.(L <= 0))then
		L = 500
	end if
	if(argCount > 1)then
		call GET_COMMAND_ARGUMENT(2, argChar)
		read(argChar,FMT=*, IOSTAT=stat)eps
		if((stat.NE.0).OR.(eps <= 0.0))then
			eps = 1E-4
		end if
	end if

	if(argCount > 2)then
		call GET_COMMAND_ARGUMENT(3, argChar)
		read(argChar,FMT=*, IOSTAT=stat)omega
		if((stat.NE.0).OR.(omega < 0.0))then
			omega = 1.0E4
		end if
	end if

	if(argCount > 3)then
		call GET_COMMAND_ARGUMENT(4, argChar)
		read(argChar,FMT=*, IOSTAT=stat)m
		if((stat.NE.0).OR.(m <= 0.0))then
			m = 1.0
		end if
	end if

	if(argCount > 4)then
		call GET_COMMAND_ARGUMENT(5, argChar)
		read(argChar,FMT=*, IOSTAT=stat)hbar
		if((stat.NE.0).OR.(hbar <= 0.0))then
			hbar = 1.0
		end if
	end if
end if
N=2*L+1
x0=sqrt(hbar/(m*omega))
kinTyp=(hbar*hbar)/(2*m*eps*eps)
H=laplmatr(N, .FALSE.)*(-kinTyp)+hosc_potmatr(eps, L, m,omega)
!print*,"Right x-cutoff:",eps*L
!print*,"Typical length:", x0
!print*,"Spacing:",eps
!print*,"Maximum order for kinetic energy:", kinTyp
!print*,"Minimum order for kinetic energy:", kinTyp/(L*L)
!print*,"Typical potential energy:", (0.5*m*omega*omega*x0*x0)
allocate(eig(N))
call hediagonalize(H,N,eig,info)
if(info.NE.0)then
	print*,'DIAGONALIZATION FAILED'
	stop
end if


fname="energies_"//trim(str(L))//"_"//trim(strr(eps))//"_"//trim(strr(omega))//"_"//trim(strr(m))//"_"//trim(strr(hbar))//".dat"
fname2="states_"//trim(str(L))//"_"//trim(strr(eps))//"_"//trim(strr(omega))//"_"//trim(strr(m))//"_"//trim(strr(hbar))//".dat"

open(unit=42, file=fname,action="write",status="replace")
write(42, "(A,I0,A,G0,A,G0,A,G0,A,G0)")"#",L,",",eps,",",omega,",",m,",",hbar
call write_vec_col(eig,N,42)
close(42)
open(unit=420, file=fname2,action="write",status="replace")
do ii=1,int(N/100)
	call write_complex_vec_col(H(:,ii),N,420)
end do
close(420)

print*,fname
print*,fname2

deallocate(eig)
deallocate(H)
end program
