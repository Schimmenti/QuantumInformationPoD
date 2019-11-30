program shiftedho
use qiutil
implicit none
complex*16, dimension(:), allocatable :: psi0
complex*16, dimension(:,:), allocatable :: psi
real*8, dimension(:,:), allocatable :: psi_r, psi_i
type(grid1) :: grd, tGrd
real*8 :: L, Tmax, TT
integer*4 :: sz, tsz
L=1.0
Tmax=1.0
TT=30*Tmax
sz=1000
tsz=1000
grd = create_grid1d(-L,L,sz)
tGrd = create_grid1d(0D1,TT,tsz)
psi0 = ho_groundstate(grd,1.0D1,1.0D1,1D2)
! in rows we find position index, in column the time index
psi = propogate_wf_hov_shifted(psi0,grd,tGrd, 1.0D1, 1.0D1, 1D4, Tmax)
psi_r = dble(psi)
psi_i = dimag(psi)
open(unit=42, file="wf_real.dat")
open(unit=420, file="wf_imag.dat")
call wrm_r8(psi_r, sz, tsz, 42)
call wrm_r8(psi_i, sz, tsz, 420)
close(42)
close(420)
end program
