module heisenberg
use, intrinsic :: iso_c_binding
implicit none

contains

subroutine dsymm_diag(matr,nn,eigvs,job,info)
    real*8, dimension(:,:), intent(inout) :: matr
    integer*4,intent(out) :: info
    integer*4,intent(in) :: nn
    character*1, intent(in) :: job
    integer*4 :: lwork
    real*8, dimension(:),allocatable,intent(inout) :: eigvs
    real*8, dimension(:), allocatable :: work
    real*8, dimension(:), allocatable :: rwork
    ! optimal lwork
    lwork=-1
    allocate(work(1))
    allocate(rwork(max(1, 3*nn-2)))
    if(.not.allocated(eigvs))then
        allocate(eigvs(nn))
    end if
    call dsyev(job,'U',nn,matr,nn,eigvs,work,lwork,rwork,info)
    lwork = int(real(work(1)))
    deallocate(work)
    deallocate(rwork)   
    ! actual diag
    allocate(work(max(1,lwork)))
    allocate(rwork(max(1, 3*nn-2)))
    call dsyev(job,'U',nn,matr,nn,eigvs,work,lwork,rwork,info)  
    deallocate(work)
    deallocate(rwork)
end subroutine

function factorial(n)result(fn)
    integer*4, intent(in) :: n
    integer*4 :: fn,idx
    fn = 1
    if(n >1)then
        do idx=2,n
            fn = fn*idx
        end do
    end if
end function

function sigma_z(state, L)result(m_z)
    integer*4, intent(in) :: state, L
    integer*4, dimension(:), allocatable :: m_z
    integer*4 :: idx
    allocate(m_z(L))
    do idx=0,L-1
        if(btest(state,idx))then
            m_z(idx+1)=1
        else
            m_z(idx+1)=-1
        end if
    end do
end function

function create_magnetizations(L)result(m_state)
    integer*4, intent(in) :: L
    integer*4, dimension(:), allocatable :: m_state
    integer*4 :: N,state
    N=2**L
    allocate(m_state(N))
    do state=0,N-1
        m_state(state+1)=sum(sigma_z(state,L))
    end do
end function

function get_spin_map(L)result(map)
    integer*4, intent(in) :: L
    integer*4, dimension(:,:), allocatable :: map
    integer*4 :: state, N
    N = lshift(1,L)
    allocate(map(N,L))
    do state=0,N-1
        map(state+1,1:L) = sigma_z(state,L)
    end do
end function

function get_sector(sector_magn, magnetizations,L)result(sec)
    integer*4, intent(in) :: L, sector_magn
    integer*4, dimension(:), intent(in) :: magnetizations
    integer*4, dimension(:), allocatable :: sec
    integer*4 :: idx,N,state, imbalance, sector_size
    N = 2**L
    imbalance = (L-abs(sector_magn))/2
    sector_size = factorial(L)/(factorial(imbalance)*factorial(L-imbalance))
    allocate(sec(sector_size))
    idx = 1 ! index to fill 'sec'
    do state=1,N
        if(magnetizations(state)==sector_magn)then
            sec(idx) = state-1
            idx = idx + 1
        end if
    end do
end function

! Builds the Heisenberg XXZ chain Hamiltonian for the given spin sector.
function hamiltonian_heisenberg(sector,L,Jxy,Jzz)result(H)
    integer*4, intent(in) :: L
    integer*4, dimension(:), intent(in) :: sector
    real*8, intent(in) :: Jxy,Jzz
    real*8, dimension(:,:), allocatable :: H
    integer*4 :: idx, jdx, dim,state,state_p,site, next_site,mask
    dim = size(sector)
    allocate(H(dim,dim))
    H=0
    do idx=1,dim
        state = sector(idx)
        do site=0,L-1
            next_site = mod(site+1,L)
            ! xy part
            if(btest(state, site).neqv.btest(state, next_site))then
                ! target state for xy interaction
                mask = lshift(1,site)+lshift(1,next_site)
                state_p =  xor(state,mask)
                do jdx=1,dim
                    if(sector(jdx).eq.state_p)then
                        H(idx,jdx) = 2*Jxy/4
                    end if
                end do
                ! target state for z interaction
                state_p =  xor(state,mask)
            end if
            ! ising part (diagonal)
            if(btest(state, site).eqv.btest(state, next_site))then
                H(idx,idx) = H(idx,idx) + Jzz/4D0
            else
                H(idx,idx) = H(idx,idx) - Jzz/4D0
            end if

        end do
    end do
    
end function

! Adds a z-local field to the Hamiltonian of a given spin sector.
subroutine add_zfield(H,sector,L,field)
    integer*4, intent(in) :: L
    integer*4, dimension(:), intent(in) :: sector
    real*8, dimension(:,:), intent(inout) :: H
    real*8, dimension(:), intent(in) :: field
    integer*4 :: dim,idx
    dim = size(sector)
    do idx=1,dim
        H(idx,idx) = H(idx,idx) + sum(sigma_z(sector(idx), L)*field)
    end do
end subroutine


function m_wave_compute(prob,sector_map, L,dim,wave_cfs)result(f)
    !  f_njk = np.einsum('sn,sj,sk',p, spin_sector, spin_sector) #(L,L,sz)
    ! f_nj = np.einsum('sn,sj', p, spin_sector) #(L,sz)
    ! M_wave[k,s] = np.mean(1-np.abs(np.matmul(wave_cfs, f_nj))**2/np.einsum('jkn,j,k', f_njk, wave_cfs, np.conj(wave_cfs))).real

    real*8, dimension(:,:), intent(in) :: prob
    integer*4, dimension(:,:), intent(in) :: sector_map
    integer*4, intent(in) :: L, dim
    complex*16, dimension(:), intent(in) :: wave_cfs
    complex*16, dimension(:,:,:), allocatable :: fnjk
    complex*16, dimension(:,:), allocatable :: fnj
    complex*16, dimension(:), allocatable :: f
    integer*4 :: idx, pos1, pos2
    allocate(f(dim))
    allocate(fnjk(dim,L,L))
    allocate(fnj(dim,L))
    do idx=1,dim
        do pos1=1,L
            fnj(idx,pos1)=sum(prob(1:dim,idx)*sector_map(1:dim, pos1))
            do pos2=1,L
                fnjk(idx,pos1,pos2) = sum(prob(1:dim,idx)*sector_map(1:dim, pos1)*sector_map(1:dim, pos2))
            end do
        end do

        f(idx) = 1 -(abs(dot_product(fnj(idx,1:L), wave_cfs))**2)/dot_product(wave_cfs, matmul(fnjk(idx,1:L,1:L), wave_cfs))
    end do
end function

function differentiate(x)result(dx)
    real*8, dimension(:), intent(in) :: x
    real*8, dimension(:), allocatable :: dx
    integer*4 :: idx, sz
    sz = size(x)
    allocate(dx(sz-1))
    do idx=1,sz-1
        dx(idx) = x(idx+1)-x(idx)
    end do
end function

end module