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
public :: create_grid1d,propagate_hov_shifted,propogate_wf_hov_shifted,ho_groundstate
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
	complex*16, dimension(:), allocatable :: temp
	integer*4 :: idx,jdx
	complex*16 :: factor
	real*8 :: t
	type(C_PTR) :: planf,planb
! since fortran is column major better to fill the matrix column by column
	allocate(psi(grd%sz, tGrd%sz))
	allocate(temp(grd%sz))
	factor = complex(0,-2*hbar*tGrd%a*(acos(-1.0)**2)/(m*(grd%a*grd%sz)**2))
	psi(:,1)=psi0

	call dfftw_plan_dft_1d(planf,grd%sz,temp,temp,FFTW_FORWARD,FFTW_ESTIMATE);
	call dfftw_plan_dft_1d(planb,grd%sz,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
	
	temp = psi0
	t = tGrd%a
	do idx=2, tGrd%sz
! potential part, first half
		call propagate_hov_shifted(temp,grd,tGrd%a/2,hbar,m,omega,t,Tmax)
! momentum space
		call dfftw_execute_dft(planf,temp,temp)
! multiply by kinetic term
		do jdx=1,grd%sz
			temp(jdx)=temp(jdx)*exp(factor*jdx*jdx)
		end do
! position space	
		call dfftw_execute_dft(planb, temp, temp)
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

function ho_groundstate(grd, hbar, m, omega)result(psi)
	type(grid1), intent(in) :: grd
	real*8, intent(in) :: hbar,m,omega
	complex*16, dimension(:), allocatable :: psi
	real*8 :: factor, factor2, x
	integer*4 :: idx
	allocate(psi(grd%sz))
	factor = (m*omega)/hbar
	factor2 = ((factor/acos(-1.0))**0.25)/sqrt(2.0)
	x = grd%xmin
	do idx=1,grd%sz
		psi(idx)=complex(factor2*exp(-0.5*factor*x*x),0.0)
		x = x + grd%a
	end do
end function

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
