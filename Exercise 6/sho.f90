module quantumutil
implicit none
contains

! second derivative and exact up to eps^2
function dx2_2(psi,eps,N)result(phi)
! lattice size
	real*4 :: eps
! number of discretized points and loop index
	integer*4 :: N,ii
! wave function
	complex*8, dimension(:) :: psi
	complex*8, dimension(:), allocatable :: phi

	allocate(phi(N))
! assume psi(-1)=0
! assume psi(N+1)=0
	phi(1)=(psi(2)-2*psi(1))/(eps*eps)
	phi(N)=(-2*psi(N)+psi(N-1))/(eps*eps)
	do ii=2,N-1
		phi(ii)=(psi(ii+1)-2*psi(ii)+psi(ii-1))/(eps*eps)
	end do
	
end function

function hoscx_2(psi, V,eps, N, omega=sqrt(2),hbar=1.0,m=0.5)result(phi)
! lattice size
	real*4 :: eps
! number of discretized points and loop index
	integer*4 :: N,ii
! mass
	real*4 :: m
! hbar
	real*4 :: hbar
! hbar
	real*4 :: hbar
! wave function
	complex*8, dimension(:) :: psi
	complex*8, dimension(:), allocatable :: phi

	phi=dx2_2(psi, eps, N)
	phi = ((-hbar*hbar)/(2*m))*phi
	do ii=1,N
		phi(ii)=V(ii*eps)*psi(ii)
	end do
end function




end module

program shotest
implicit none

end program
