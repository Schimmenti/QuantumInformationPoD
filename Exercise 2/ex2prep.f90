module matrices
type dmatrix
integer, dimension(2) :: N
double complex, dimension(:,:), allocatable :: elem

end type dmatrix



interface operator (.Adj.)
	module procedure matadjoint,vectoradjoint
end interface


contains

function matadjoint(A)
	
end function

function vectoradjoint(v)

end function

end module

program main
implicit none
use matrices
type(dmatrix) :: A,b

A%Elem = 0d0
B =.Adj.A

end program
