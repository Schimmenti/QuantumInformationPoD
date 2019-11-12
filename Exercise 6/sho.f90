module quantumutil
implicit none
contains

! second derivative and exact up to eps^2
function dx2_2(psi,eps,N, pbc=.TRUE.)result(phi)
! lattice size
	real*4 :: eps
! number of discretized points and loop index
	integer*4 :: N,ii
! wave function
	complex*8, dimension(:) :: psi
	complex*8, dimension(:), allocatable :: phi

	allocate(phi(N))

	if(pbc)then
		phi(1)=(psi(2)-2*psi(1)+psi(N))/(eps*eps)
		phi(N)=(psi(1)-2*psi(N)+psi(N-1))/(eps*eps)
	else
		phi(1)=(psi(2)-2*psi(1))/(eps*eps)
		phi(N)=(-2*psi(N)+psi(N-1))/(eps*eps)
	end if
	do ii=2,N-1
		phi(ii)=(psi(ii+1)-2*psi(ii)+psi(ii-1))/(eps*eps)
	end do
	
end function

function laplmatr(eps, N, pbc=.TRUE.)result(A)
! lattice size
	real*4 :: eps
! number of discretized points and loop index
	integer*4 :: N,ii

	real*4, dimension(:,:), allocatable :: A
	allocate(A(N,N))
	A=0
	if(pbc)
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
	
end function

function hoscx_2(psi,eps, L, omega=sqrt(2),hbar=1.0,m=0.5)result(phi)
! lattice size
	real*4 :: eps
! number of discretized points and loop index
	integer*4 :: N,L,ii
! mass
	real*4 :: m
! hbar
	real*4 :: hbar
! hbar
	real*4 :: hbar
! wave function
	complex*8, dimension(:) :: psi
	complex*8, dimension(:), allocatable :: phi
	N=2*L+1
	phi=dx2_2(psi, eps, N)
	phi = ((-hbar*hbar)/(2*m))*phi
	do ii=-L,L
		phi(ii)=phi(ii)+(0.5*omega*omega)*((ii*eps)**2)*psi(ii+L+1)
	end do
end function






end module

program shotest
implicit none

end program
