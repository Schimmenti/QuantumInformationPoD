program test
use nbody
use qiutil
implicit none
integer*4 :: Nmax, nn,ll,info
type(grid1) :: wGrd
complex*16, dimension(:,:), allocatable :: H
real*8, dimension(:), allocatable :: eigs
character(:), allocatable :: fname
real*8 :: ww
wGrd = create_grid1d(0.0D0,3.0D0,20)
Nmax=12
do nn=2,Nmax
    print*,'Doing',nn,'positions...'
    ww = wGrd%xmin
    fname = 'anderson_'//trim(str(nn))//'.txt'
    open(unit=42, file=fname)
    do ll=1,wGrd%sz+1
        write(42,'(G0,A)', advance='no')ww,' '
        H = andersonHamiltonian(nn, ww)
        call herm_diag(H, nn, eigs, 'V', info)
        call write_vec_row(eigs, nn, 42)
        deallocate(H)
        deallocate(eigs)
        ww = ww + wGrd%a
    end do
    close(42)
end do
end program