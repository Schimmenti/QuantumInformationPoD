module matrixutil
implicit none
integer :: ii,jj,kk

contains

! fills a NxN matrix with VAL
function FILLMATF(N, VAL)result(x)
	integer :: N
	real :: VAL
	real*4, dimension(N,N) :: x
	do ii=0,N*N-1
		x(mod(ii,N)+1,(ii/N)+1)=VAL
	end do
end function

! fills a NxN matrix with random [0,1] values
function RANDMATF(N)result(x)
	integer :: N
	real*4, dimension(N,N) :: x
	do ii=0,N*N-1
		x(mod(ii,N)+1,(ii/N)+1)=RAND(0)
	end do
end function

! matrix multiplication - c(i,j) = \sum_k a(i,k)*b(k,j) - the index i runs slowest
function IJKMULF(a, b, N)result(c)
	INTEGER  :: N
	REAL*4, DIMENSION(:,:) :: a, b
	REAL*4, DIMENSION(N,N) :: c
	do ii=1,N
		do jj=1,N
			do kk=1,N
				c(ii,jj)=c(ii,jj) + a(ii,kk)*b(kk,jj)
			end do
		end do
	end do
end function IJKMULF

! matrix multiplication - c(i,j) = \sum_k a(i,k)*b(k,j) - the index i runs fastest
function KJIMULF(a, b, N)result(c)
	INTEGER  :: N
	REAL*4, DIMENSION(:,:) :: a, b
	REAL*4, DIMENSION(N,N) :: c
	do kk=1,N
		do jj=1,N
			do ii=1,N
				c(ii,jj)=c(ii,jj) + a(ii,kk)*b(kk,jj)
			end do
		end do
	end do
end function KJIMULF

end module


program testmmul
use matrixutil
implicit none
  integer nn, stat
  character(10) :: nnchar
  real :: temp
  real :: start, finish
  real, dimension (:,:), allocatable :: matr1,matr2, res1, res2, res3
  IF(COMMAND_ARGUMENT_COUNT() < 1)THEN
!    Using default matrix size
     nn = 100
  ELSE
      CALL GET_COMMAND_ARGUMENT(1, nnchar)
      READ(nnchar,*)nn
  ENDIF
! we allocate the matrices with the given size; this is done in order
! to allow dynamic sizing
  allocate(matr1(nn,nn))
  allocate(matr2(nn,nn))
  allocate(res1(nn,nn))
  allocate(res2(nn,nn))
  allocate(res3(nn,nn))
! populate the the two matrices using random [0,1] values
! and the resulting ones using zeros
  matr1 = RANDMATF(nn)
  matr2 = RANDMATF(nn)
  res1 = FILLMATF(nn,0e0)
  res2 = FILLMATF(nn,0e0)
  res3 = FILLMATF(nn,0e0)
  start = 0
  finish = 0
! first method of multiplication
  call cpu_time(start)
  res1 = IJKMULF(matr1, matr2, nn)
  call cpu_time(finish)
  print*,finish-start
! second method of multiplication
  call cpu_time(start)
  res2 = KJIMULF(matr1, matr2, nn)
  call cpu_time(finish)
  print*,finish-start
! intrinsic method of multiplication
  call cpu_time(start)
  res3 = matmul(matr1,matr2)
  call cpu_time(finish)
  print*,finish-start
! deallocate all used matrices
  deallocate(matr1)
  deallocate(matr2)
  deallocate(res1)
  deallocate(res2)
  deallocate(res3)
  end program
