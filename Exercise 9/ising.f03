program test
use nbody
use qiutil
implicit none
integer*4 :: Nmax, nn,ll,info
type(grid1) :: lGrd
complex*16, dimension(:,:), allocatable :: H
real*8, dimension(:), allocatable :: eigs
real*8 :: lmbd
character(:), allocatable :: fname
lGrd = create_grid1d(0.0D0,3.0D0,20)
Nmax=9
do nn=2,Nmax
    print*,'Doing',nn,'spins...'
    lmbd = lGrd%xmin
    fname = 'isingsp_'//trim(str(nn))//'.txt'
    open(unit=42, file=fname)
    do ll=1,lGrd%sz+1
        write(42,'(G0,A)', advance='no')lmbd,' '
        H = isingHamiltonian(nn, lmbd)
        call herm_diag(H, 2**nn, eigs, 'V', info)
        call write_vec_row(eigs, 2**nn, 42)
        deallocate(H)
        deallocate(eigs)
        lmbd = lmbd + lGrd%a
    end do
    close(42)
end do
end program