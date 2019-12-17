program test
use nbody
use qiutil
implicit none
integer*4 :: nAverages, nSites,nn,ll,info
type(grid1) :: wGrd
complex*16, dimension(:,:), allocatable :: H
real*8, dimension(:), allocatable :: eigs, avgEig
character(:), allocatable :: fname
real*8 :: ww
wGrd = create_grid1d(0.0D0,3.0D0,20)
nSites =250
nAverages = 100
fname = 'anderson_'//trim(str(nSites))//'.txt'
!loop though possible W
allocate(avgEig(nSites))
open(unit=42, file=fname)
do ll=1,wGrd%sz+1
    avgEig = 0
    do nn=1, nAverages
        H = andersonHamiltonian(nSites, ww)
        call herm_diag(H, nSites, eigs, 'N', info)
        ! call write_vec_row(eigs, nSites, 42)
        avgEig = avgEig + eigs
        deallocate(H)
        deallocate(eigs)
    end do
    write(42, '(G0,A)', advance='no')ww,' '
    call write_vec_row(avgEig/nAverages, nSites, 42)
    print*, ll,'/',wGrd%sz+1
    ww = ww + wGrd%a
end do
close(42)

end program