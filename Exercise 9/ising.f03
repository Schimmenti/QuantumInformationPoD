program test
use nbody
use qiutil
implicit none
integer*4 :: Nmax, nn,ll,info
type(grid1) :: lGrd
complex*16, dimension(:,:), allocatable :: H
real*8, dimension(:), allocatable :: eigs
real*8 :: lmbd
character(:), allocatable :: fname,fname2
lGrd = create_grid1d(0.0D0,3.0D0,50)
Nmax=10
do nn=2,Nmax
    print*,'Doing',nn,'spins...'
    lmbd = lGrd%xmin
    fname = 'isingsp_'//trim(str(nn))//'.txt'
    fname2 = 'isinggs_'//trim(str(nn))//'.txt'
    open(unit=42, file=fname)
    open(unit=24, file=fname2)
    do ll=1,lGrd%sz+1
        write(42,'(G0,A)', advance='no')lmbd,' '
        H = isingHamiltonian(nn, lmbd)
        call herm_diag(H, 2**nn, eigs, 'V', info)
        call write_vec_row(eigs, 2**nn, 42)

        write(24,'(G0,A)', advance='no')lmbd,' '
        call write_vec_row(real(H(:,1)), 2**nn, 24)
        write(24,'(G0,A)', advance='no')lmbd,' '
        call write_vec_row(aimag(H(:,1)), 2**nn, 24)
        
        deallocate(H)
        deallocate(eigs)
        lmbd = lmbd + lGrd%a
    end do
    close(42)
    close(24)
end do
end program