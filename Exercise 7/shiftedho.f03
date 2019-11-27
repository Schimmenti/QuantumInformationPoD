program shiftedho
use qiutil
implicit none
complex*16, dimension(:), allocatable :: psi0
complex*16, dimension(:,:), allocatable :: psi
type(grid1) :: grd
real*8 :: L
integer*4 :: sz

grd = create_grid1d(-L,L,sz)
psi0 = ho_groundstate(grd,1.0D1,1.0D1,1D4)

	

end program
