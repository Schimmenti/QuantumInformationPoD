program shiftedho
use qiutil
implicit none
complex*16, dimension(:), allocatable :: psi0
complex*16, dimension(:,:), allocatable :: psi, states
real*8, dimension(:,:), allocatable :: psi_r, psi_i
type(grid1) :: grd, tGrd
real*8 :: L, Tmax, TT, hbar, m, omega
integer*4 :: sz, tsz
L=5.0
Tmax=1.0
TT=3*Tmax
sz=1000
tsz=300
hbar=1.0
m=1.0
omega=1D1
grd = create_grid1d(-L,L,sz)
tGrd = create_grid1d(0D1,TT,tsz)
states = ho_diag(grd,hbar,m,omega)
psi0 = states(:,1)
deallocate(states)
! in rows we find position index, in column the time index
psi = propogate_wf_hov_shifted(psi0,grd,tGrd, hbar, m, omega, Tmax)
psi_r = dble(psi)
psi_i = dimag(psi)
open(unit=42, file="wf_real.dat")
open(unit=420, file="wf_imag.dat")
open(unit=4200, file="grids.dat")
call wrm_r8(psi_r, sz, tsz, 42)
call wrm_r8(psi_i, sz, tsz, 420)
call print_grid1(4200, grd)
call print_grid1(4200, tGrd)
write(4200, '(G0,A,G0,A,G0,A,G0)')hbar,',', m,',',omega,',',Tmax
close(42)
close(420)
close(4200)
end program
