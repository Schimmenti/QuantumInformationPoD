module qiutil
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'
include 'fftw3l.f03'
type grid1
	real*8 :: xmin
	real*8 :: xmax
	integer*4 :: sz
	real*8 :: a
end type grid1
contains

function create_grid1d(minPos, maxPos, nPoints)result(grd)
	real*8, intent(in) :: minPos, maxPos	
	integer*4, intent(in) :: nPoints
	type(grid1) :: grd
	grd%xmin = minPos
	grd%xmax = maxPos
	grd%sz = nPoints
	grd%a = (maxPos-minPos)/nPoints
	return
end function

subroutine print_grid1(unit,grd)
	integer*4, intent(in) :: unit 
	type(grid1), intent(in) :: grd
	write(unit, '(G0,A,G0,A,G0,A,I0)')grd%xmin,',',grd%xmax,',',grd%a,',',grd%sz
end subroutine

function laplmatr(N, pbc)result(A)
	integer*4 :: N
	logical :: pbc


	integer*4 :: ii
	complex*16, dimension(:,:), allocatable :: A

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

function hosc_potmatr(m, omega, grd)result(Vpot)
	type(grid1), intent(in) :: grd
	real*8 :: x_ii, m, omega
	complex*16, dimension(:,:), allocatable :: Vpot
	integer*4 :: ii
	allocate(Vpot(grd%sz,grd%sz))
	Vpot=0.0
	x_ii = grd%xmin
	do ii=1,grd%sz
		Vpot(ii,ii)=(0.5*m)*((omega*x_ii)**2)
		x_ii = x_ii + grd%a
	end do
	return
end function

subroutine propagate_hov_shifted(psi,grd,tau,hbar,m,omega,t,Tmax)
! args
	complex*16, dimension(:), intent(inout) :: psi
	type(grid1), intent(in) :: grd
	real*8, intent(in) :: tau,hbar,m,omega, t, Tmax
! dummy vars
	complex*16 :: factor
	real*8 :: x_sh
	integer*4 :: idx
	factor = complex(0,-tau*0.5*m*omega*omega/hbar)
	x_sh = grd%xmin-t/Tmax
	do idx=1,grd%sz
		psi(idx)=psi(idx)*exp(factor*x_sh**2)
		x_sh = x_sh + grd%a
	end do
end subroutine

function propogate_wf_hov_shifted(psi0,grd,tGrd,hbar,m,omega,Tmax)result(psi)
! args
	complex*16, dimension(:), intent(in) :: psi0
	type(grid1), intent(in) :: grd, tGrd
	real*8, intent(in) :: hbar,m,omega,Tmax
! dummy vars
	complex*16, dimension(:,:), allocatable :: psi
	complex*16, dimension(:), allocatable :: temp,temp2
	integer*4 :: idx,jdx
	complex*16 :: factor
	real*8 :: t
	type(C_PTR) :: planf,planb
! since fortran is column major better to fill the matrix column by column
	allocate(psi(grd%sz, tGrd%sz))
	allocate(temp(grd%sz))
	allocate(temp2(grd%sz))
	factor = complex(0, -0.5*tGrd%a*4*acos(-1.0)/(hbar*m*(grd%a*grd%sz)**2))
	psi(:,1)=psi0

	call dfftw_plan_dft_1d(planf,grd%sz,temp,temp2,FFTW_FORWARD,FFTW_ESTIMATE);
	call dfftw_plan_dft_1d(planb,grd%sz,temp2,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
	
	temp = psi0
	t = 0.0
	do idx=2, tGrd%sz
! potential part, first half
		call propagate_hov_shifted(temp,grd,tGrd%a/2,hbar,m,omega,t,Tmax)
! momentum space
		call dfftw_execute_dft(planf,temp,temp2)
! multiply by kinetic term
		do jdx=0,grd%sz/2-1
			temp2(jdx+1)=temp2(jdx+1)*exp(factor*jdx*jdx)
		end do

		do jdx=grd%sz/2, grd%sz-1
			temp2(jdx+1)=temp2(jdx+1)*exp(factor*(grd%sz-jdx)**2)
		end do
! position space	
		call dfftw_execute_dft(planb, temp2, temp)
! potential part, first half
		call propagate_hov_shifted(temp,grd,tGrd%a/2,hbar,m,omega,t,Tmax)
		psi(:,idx)=temp/grd%sz
		temp = psi(:,idx)
		t = t + tGrd%a
	end do
	call dfftw_destroy_plan(planf)
	call dfftw_destroy_plan(planb)
	deallocate(temp)
end function

subroutine herm_diag(matr,nn,eigvs,job,info)
	complex*16, dimension(:,:), intent(inout) :: matr
	integer*4,intent(out) :: info
	integer*4,intent(in) :: nn
	character*1, intent(in) :: job
	integer*4 :: lwork
	real*8, dimension(:),allocatable,intent(inout) :: eigvs
	complex*16, dimension(:), allocatable :: work
	real*8, dimension(:), allocatable :: rwork
	! optimal lwork
	lwork=-1
	allocate(work(1))
	allocate(rwork(max(1, 3*nn-2)))
	if(.not.allocated(eigvs))then
		allocate(eigvs(nn))
	end if
	call zheev(job,'U',nn,matr,nn,eigvs,work,lwork,rwork,info)
	lwork = int(real(work(1)))
	deallocate(work)
	deallocate(rwork)


	! actual diag
	allocate(work(max(1,lwork)))
	allocate(rwork(max(1, 3*nn-2)))
	call zheev(job,'U',nn,matr,nn,eigvs,work,lwork,rwork,info)

	deallocate(work)
	deallocate(rwork)
end subroutine


function ho_diag(grd, hbar, m, omega)result(H)
	type(grid1), intent(in) :: grd
	real*8, intent(in) :: hbar,m,omega
	complex*16, dimension(:,:), allocatable :: H
	real*8, dimension(:), allocatable :: eig
	integer*4 :: nn
	real*8 :: kinTyp, x0

	! diag params
	integer*4::info




	nn=grd%sz
	allocate(eig(nn))
	
	! quantum stuff

	x0=sqrt(hbar/(m*omega))
	kinTyp=(hbar*hbar)/(2*m*grd%a*grd%a)
	H=laplmatr(nn, .FALSE.)*(-kinTyp)+hosc_potmatr(m,omega, grd)

	! diag

	call herm_diag(H, nn, eig, 'V', info)
end function

subroutine write_vec_col(x,sz,fileUnit)
	real*8, dimension(:) :: x
	integer*4 :: sz, fileUnit
	integer*4 :: ii
	do ii=1,sz
		write(fileUnit,"(G0,A)")x(ii)
	end do
end subroutine

subroutine write_vec_row(x,sz,fileUnit)
	real*8, dimension(:) :: x
	integer*4 :: sz, fileUnit
	integer*4 :: ii
	do ii=1,sz-1
		write(fileUnit,"(G0,A)",advance='no')x(ii),' '
	end do
	write(fileUnit,"(G0,A)")x(sz)
end subroutine

subroutine write_complex_vec_col(x,sz,fileUnit)
	complex*16, dimension(:) :: x
	integer*4 :: sz, fileUnit
	integer*4 :: ii
	do ii=1,sz
		write(fileUnit,"(G0,A,G0)")real(x(ii)),',',aimag(x(ii))
	end do
end subroutine

subroutine wrm_r8(x, nRows, nCols, fileUnit)
	real*8, dimension(:,:), intent(in) :: x
	integer*4, intent(in) :: nRows, nCols, fileUnit
	integer*4 :: ii,jj
	do ii=1, nRows
		do jj=1,nCols-1
			write(fileUnit,"(G0,A)", advance="no")x(ii,jj),','
		end do
		write(fileUnit, "(G0)")x(ii,nCols)
	end do
end subroutine

subroutine wrm_c16(x, nRows, nCols, fileUnit)
	complex*16, dimension(:,:), intent(in) :: x
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
 	real*8, dimension(:) :: x
	integer*4 :: sz, shift

	real*8, dimension(:), allocatable :: s
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
