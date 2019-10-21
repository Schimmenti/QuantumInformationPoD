module dcmat

implicit none

type dcmatrix
integer, dimension(2) :: n
double complex, dimension(:,:), allocatable :: els
double complex :: mtrace
double complex :: mdet
end type dcmatrix

interface operator (.adj.)
	module procedure dcmatadjoint
end interface

interface operator (.trace.)
	module procedure dcmattrace
end interface


contains

function createmat(rows,cols)
	integer, intent(in) :: rows, cols
	type(dcmatrix) :: createmat
	createmat%n(1)=rows
	createmat%n(2)=cols
	allocate(createmat%els(rows, cols))
	return
end function

function createrndmat(rows,cols)
	integer, intent(in) :: rows, cols
	integer :: iindex
	type(dcmatrix) :: createrndmat
	createrndmat = createmat(rows, cols)
	do iindex=0,rows*cols-1
		createrndmat%els(mod(iindex,rows)+1, (iindex/rows)+1)=complex(RAND(0),RAND(0))
	end do
	return
end function


function dcmatadjoint(mat)result(matres)
	type(dcmatrix), intent(in) :: mat
	type(dcmatrix) :: matres
	matres%n = mat%n
	matres%els = conjg(transpose(mat%els))
	return
end function

function dcmattrace(mat)result(res)
	type(dcmatrix), intent(in) :: mat
	double complex :: res
	integer :: iindex
! we use a safe implementation to avoid errors in non square matrices
	integer :: maxindex
	maxindex = min(mat%n(1), mat%n(2))
	res = complex(0d0,0d0)
	do iindex=1,maxindex
		res = res + mat%els(iindex,iindex)
	end do
	return
end function

end module dcmat

program dcmatest
use dcmat
type(dcmatrix) :: A

A = createrndmat(10,10)

print*,.trace.(A)
end program
