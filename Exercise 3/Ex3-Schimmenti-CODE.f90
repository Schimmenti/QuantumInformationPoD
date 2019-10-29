MODULE DEBUGUTIL
IMPLICIT NONE
CONTAINS
! ########### DOCUMENTATION ###########
! DEBUGPROC prints a debug message using different options
! ########### ARGUMENTS ###########
! DEBUG: logical. If .TRUE. means that we are in debug mode so we might want an error or debug message.
! CONDITION: logical, optional. If .TRUE. error has occurred; it not set it is assumed to be .TRUE..
! MSG: character(:). The error\debug message
! CONTENT: any primitive type (non vector\matrix). The variable to be printed additionally.
! FRAME: logical, optional. If .TRUE. a small frame is added to the message
! STOPONERROR: logical, optional. If .TRUE. and if DEBUG AND CONDITION are .TRUE. stops the program.
! ########### CONDITIONS ###########
! PRE: N is an integer*4, VAL is a real*4
! POST: x is a dimension(N,N) matrix with all values equal to VAL if N > 0 OR a zero-dimensional matrix if N <= 0
	SUBROUTINE DEBUGPROC(DEBUG,CONDITION,MSG,CONTENT,FRAME,STOPONERROR)
		LOGICAL, INTENT(IN) :: DEBUG
		LOGICAL :: PROCEED, USEFRAME
		LOGICAL, INTENT(IN), OPTIONAL :: CONDITION
		LOGICAL, INTENT(IN), OPTIONAL :: FRAME
		LOGICAL, INTENT(IN), OPTIONAL :: STOPONERROR
		CHARACTER(*), INTENT(IN) :: MSG
		CLASS(*), INTENT(IN), OPTIONAL :: CONTENT

		
		PROCEED = (PRESENT(CONDITION).AND.(DEBUG.AND.CONDITION)).OR.(DEBUG)
		USEFRAME = (PRESENT(FRAME).AND.FRAME).OR.(.NOT.(PRESENT(FRAME)))

		IF (PROCEED) THEN
			if(USEFRAME)then
				PRINT*,"------------------------------------------------"
			end if
			PRINT*,MSG
			if(PRESENT(CONTENT))then
			select type(CONTENT)
				type is (integer(1))
				PRINT*,CONTENT
				type is (integer(2))
				PRINT*,CONTENT
				type is (integer(4))
				PRINT*,CONTENT
				type is (integer(8))
				PRINT*,CONTENT
				type is (real(4))
				PRINT*,CONTENT
				type is (real(8))
				PRINT*,CONTENT
				type is (logical)
				PRINT*,CONTENT
			end select
			end if
			if(USEFRAME)then
				PRINT*,"------------------------------------------------"
			end if
			IF(PRESENT(CONDITION).AND.PRESENT(STOPONERROR).AND.STOPONERROR)THEN
				STOP
			END IF
		END IF
	END SUBROUTINE	
END MODULE


module matrixutil
use DEBUGUTIL
implicit none
! indices for looping through array elements
integer :: ii,jj,kk

! we set debugMode == .TRUE. if we want to get messages on errors
logical :: debugMode
contains

! ########### DOCUMENTATION ###########
! CREATEMATF creates a real*4 N*N matrix
! with all values equal to VAL
! ########### ARGUMENTS ###########
! N: integer*4
! VAL: real*4
! ########### CONDITIONS ###########
! PRE: N is an integer*4, VAL is a real*4
! POST: x is a dimension(N,N) matrix with all values equal to VAL if N > 0 OR a zero-dimensional matrix if N <= 0
function CREATEMATF(N, VAL)result(x)
	integer :: N
	real*4 :: VAL
	real*4, dimension(:,:),allocatable :: x
! the matrix dimension must be positive in order to allocate
! an actual matrix
	if(N <= 0)then
		allocate(x(0,0))
		call DEBUGPROC(DEBUG=debugMode,MSG="Dimension is not valid")
	else
		allocate(x(N,N))
		x = VAL
	end if
end function

! ########### DOCUMENTATION ###########
! RANDMATF creates a real*4 N*N matrix
! with all values extracted pseudo-randomly from [0,1] using RAND(0)
! ########### ARGUMENTS ###########
! N: integer*4
! ########### CONDITIONS ###########
! PRE: N is an integer*4, VAL is a real*4
! POST: x is a dimension(N,N) matrix with all values between 0 and 1 if N > 0 OR a zero-dimensional matrix if N <= 0
function RANDMATF(N)result(x)
	integer :: N
	real*4, dimension(:,:), allocatable :: x
! the matrix dimension must be positive in order to allocate
! an actual matrix
	if(N <= 0)then
		allocate(x(0,0))
		call DEBUGPROC(DEBUG=debugMode,MSG="Dimension is not valid")
	else
		allocate(x(N,N))
! since each element is different from the other we need to
! loop through all the elements; one cannot set "x = RAND(0)" since
! RAND is not vectorized and would set all the array elements equal to the same
! first random value
		do jj=1,N
			do ii=1,N
				x(ii,jj)=RAND(0)
			end do
		end do
	end if
end function

! ########### DOCUMENTATION ###########
! IJKMULF multiplies two N*N matrices
! using the order given by the argument 'mode':
! (leftmost index is the slowest running)
! 0 => IJK (default)
! 1 => IKJ
! 2 => JIK
! 3 => JKI
! 4 => KIJ
! 5 => KJI
! c(i,j)=\sum_{k=1}^N a(i,k)*b(k,j)
! ########### ARGUMENTS ###########
! a: real*4 N*N matrix
! b: real*4 N*N matrix
! N: integer*4
! mode: integer*1
! ########### CONDITIONS ###########
! PRE: 'mode' is an integer*1 AND (('N' is a positive integer*4
! 'a' and 'b' are two real*4 (at least) N*N matrices) OR ('N' is a non-positive integer
! and 'a' and 'b' are two real*4 unknown size matrices))
! POST: if N > 0 'c' is the result of matrix multiplication 'a'*'b'
! if N <= 0 'c' is an zero-dimensional matrix
function MMULF(a, b, N, mode)result(c)
	integer*4  :: N
	integer*1 :: mode
	real*4, dimension(:,:) :: a, b
	real*4, dimension(:,:), allocatable :: c
	real*4 :: temp
	if(N > 0)then
		allocate(c(N,N))
		c = 0E0
		select case(mode)
! IKJ
		case (1)
		do ii=1,N
		do kk=1,N
		temp = a(ii,kk)
		do jj=1,N
		c(ii,jj)=c(ii,jj) + temp*b(kk,jj)
		end do
		end do
		end do
! JIK
		case (2)
		do jj=1,N
		do ii=1,N
		temp = 0
		do kk=1,N
		temp = temp + a(ii,kk)*b(kk,jj)
		end do
		c(ii,jj)=temp
		end do
		end do
! JKI
		case (3)
		do jj=1,N
		do kk=1,N
		temp = b(kk,jj)
		do ii=1,N
		c(ii,jj)=c(ii,jj) + a(ii,kk)*temp
		end do
		end do
		end do
! KIJ
		case (4)
		do kk=1,N
		do ii=1,N
		temp = a(ii,kk)
		do jj=1,N
		c(ii,jj)=c(ii,jj) + temp*b(kk,jj)
		end do
		end do
		end do
! KJI
		case (5)
		do kk=1,N
		do jj=1,N
		temp = b(kk,jj)
		do ii=1,N
		c(ii,jj)=c(ii,jj) + a(ii,kk)*temp
		end do
		end do
		end do
! IJK
		case default
		do ii=1,N
		do jj=1,N
		temp = 0
		do kk=1,N
		temp = temp + a(ii,kk)*b(kk,jj)
		end do
		c(ii,jj)=temp
		end do
		end do
		end select
	else
		allocate(c(0,0))
		call DEBUGPROC(DEBUG=debugMode,MSG="Dimension is not valid")
	end if
end function MMULF





end module



program matrixmul
use matrixutil
implicit none
integer nn, stat
integer*1 :: mode
character(10) :: nnchar
real :: temp
real :: start, finish
real, dimension (:,:), allocatable :: mat1,mat2,res
character(len=20) :: modes(6)

debugMode = .TRUE.


! getting the matrix size from command line
! too few arguments i.e. we assume that we are going with the default matrix size which is 100*100
if(COMMAND_ARGUMENT_COUNT() < 1)then
	call DEBUGPROC(DEBUG=debugMode,MSG="No size specified: using default matrix size (100)")
	nn = 100
else
! get the only command line argument given
	call GET_COMMAND_ARGUMENT(1, nnchar)
	read(nnchar,FMT=*, IOSTAT=stat)nn
! if stat is different from zero so we fall back into invalid state and we use default matrix size
	if((stat.NE.0).OR.(nn <= 0))then
		call DEBUGPROC(DEBUG=debugMode,MSG="No size specified: using default matrix size (100)")
		nn = 100
	end if
	
endif

! matrix generation
mat1 = RANDMATF(nn)
mat2 = RANDMATF(nn)

! testing the different implementations of matrix multiplication
modes = ["IJK", "IKJ", "JIK", "JKI", "KIJ","KJI"]
do mode=0,5
call cpu_time(start)
res = MMULF(mat1,mat2,nn,mode)
call cpu_time(finish)
print*, "Mode: ", modes(mode+1),finish-start
end do

! measuring intrinsic fortran multiplication speed
call cpu_time(start)
res = matmul(mat1,mat2)
call cpu_time(finish)
print*, "Mode: intrinsic ",finish-start
! deallocate resources
deallocate(mat1)
deallocate(mat2)
deallocate(res)
  end program
