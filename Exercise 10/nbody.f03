module nbody
implicit none
type pstate
    integer*4 :: dim
    integer*4 :: np
    logical :: sep
    complex*16, dimension(:), allocatable :: psi
end type

contains
! Allocates a pure separable state. The memory occupied is sizeof(complex*16)*dimH*nParticles
function allocate_pure_separable(dimH, nParticles)result(state)
    integer*4, intent(in) :: dimH, nParticles
    type(pstate) :: state
    state%dim = dimH
    state%np = nParticles
    allocate(state%psi(dimH*nParticles))
    state%sep =.TRUE.
end function
! Allocates a pure state. The memory occupied is sizeof(complex*16)*dimH**nParticles
function allocate_pure(dimH, nParticles)result(state)
    integer*4, intent(in) :: dimH, nParticles
    type(pstate) :: state
    state%dim = dimH
    state%np = nParticles
    allocate(state%psi(dimH**nParticles))
    state%sep =.FALSE.
end function
! Creates a random pure state.
function random_pure_state(dimH, nParticles, separable)result(state)
    integer*4, intent(in) :: dimH, nParticles
    logical, intent(in) :: separable
    real*8, dimension(:), allocatable :: temp, temp2
    real*8 :: norm
    type(pstate) :: state
    state%dim = dimH
    state%np = nParticles
    state%sep = separable
    if(separable)then
        allocate(temp(dimH*nParticles))
        allocate(temp2(dimH*nParticles))
    else
        allocate(temp(dimH**nParticles))
        allocate(temp2(dimH**nParticles))
    end if
    call random_number(temp)
    call random_number(temp2)
    norm = sum(temp**2+temp2**2)
    state%psi = temp/norm
    deallocate(temp)
    state%psi = state%psi+temp2*complex(0,1/norm)
    deallocate(temp2)
end function

function trace_c16(mat, sz)result(tr)
    complex*16, dimension(:,:), intent(in) :: mat
    integer*4, intent(in) :: sz
    complex*16 :: tr
    integer*4 :: ii
    tr = complex(0.0,0.0)
    do ii=1,sz
        tr = tr + mat(ii,ii)
    end do
end function

function get_element_pstate(sIndex,state)result(coeff)
! this stands for basis element
! it is a d-dimensional vector
! with d = pstate%dim
    integer*4, dimension(:), intent(in) :: sIndex
    type(pstate), intent(in) :: state
    complex*16 :: coeff
    integer*4 :: d,N,ii,jj,accum
    d = state%dim
    N = state%np
    ! if separable the state is made in the following way:
    ! psi1,1...psi1,d|psi2,1,...,psi2,d|...|psiN1,...,psiNd
    if(state%sep)then
        coeff = 1.0
        do ii=0,N-1
            coeff = coeff*state%psi(ii*d+sIndex(ii+1)+1)
        end do
    else
    ! from the general state, we can view the index as
    ! a_N ... a_1 -> (a_N * d^(N-1)+...+a_2*d+a_1)+1
    ! the +1 is due to fortran indexing
        accum=1
        do ii=1, N
            jj = jj + accum*sIndex(ii)
            accum = accum*d
        end do
        coeff = state%psi(jj+1)
    end if
end function

! set the element where the basis is |sIndex[1]>...|sIndex[N]>
! not working for separable states
subroutine set_element_pstate_gen(sIndex,state, coeff)
    ! this stands for basis element
    ! it is a d-dimensional vector
    ! with d = pstate%dim
        integer*4, dimension(:), intent(in) :: sIndex
        type(pstate), intent(inout) :: state
        complex*16, intent(in) :: coeff
        integer*4 :: d,N,ii,jj,accum
        d = state%dim
        N = state%np
        if(state%sep.EQV..FALSE.)then
            accum=1
            do ii=1, N
                jj = accum*sIndex(ii)
                accum = accum*d
            end do
            state%psi(jj+1)=coeff
        end if
end subroutine
! set the coefficents at the particle 'pIndex'
! only working for separable states
subroutine set_element_pstate_sep(pIndex,state, coeff)
    ! this stands for basis element
    ! it is a d-dimensional vector
    ! with d = pstate%dim
    integer*4, intent(in) :: pIndex
    type(pstate), intent(inout) :: state
    complex*16, dimension(:),intent(in) :: coeff
    integer*4 :: d,N
    d = state%dim
    N = state%np
    if(state%sep)then
        state%psi((pIndex-1)*d+1:pIndex*d)=coeff
    end if
end subroutine
! Performs the outer product between the vector 'x' and 'y'. The 'y' is conjugated.
function vector_oproduct(x,y)result(rho)
    complex*16, dimension(:), intent(in) :: x,y
    complex*16, dimension(:,:), allocatable :: rho
    integer*4 :: nn, mm
    nn = size(x)
    mm = size(y)
    rho = matmul(reshape(x, (/nn,1/)), reshape(conjg(y), (/1,mm/)))
end function
! Performs the tensor product between two matrices 'a' and 'b'
function tensor_product(a, b, aSz, bSz)result(tt)
    complex*16, dimension(:,:), intent(in) :: a,b
    integer*4, intent(in) :: aSz, bSz
    complex*16, dimension(:, :), allocatable :: tt
    integer*4 :: ii, jj
    allocate(tt(aSz*bSz, aSz*bSz))
    do ii=1, aSz
        do jj=1, aSz
            tt((ii-1)*bSz+1:ii*bSz, (jj-1)*bSz+1:jj*bSz)=a(ii,jj)*b
        end do
    end do
end function
! general tracing
! if we employ qubits one can simplify the code
! to make use of logical operations
function trace_dm(rho, nn, dim, sSystem)result(rhop)
    complex*16, dimension(:,:), intent(in) :: rho
    integer*4, intent(in) :: nn, dim, sSystem
    complex*16, dimension(:,:), allocatable :: rhop
    complex*16 :: temp
    integer*4 :: sz,szp, factor, rightMax, leftMax, lli, rri,llj, rrj, val, ii,jj,iip, jjp
    sz = dim**nn
    szp = dim**(nn-1)
    
    ! the index is created from 'nn' elements, each going from 0 to 'dim-1'
    ! to trace an index we loop to the left and right part of sSystem
    ! (leftPart)(sSystem)(rightPart)
    ! the left is made of nn-sSystem elements
    ! the right by sSystem - 1
    rightMax = dim**(sSystem-1)
    leftMax = dim**(nn-sSystem)
    allocate(rhop(szp, szp))
    ! the factor is used to joint the two indices
    factor = dim**(sSystem-1)
    ! the period, in the indexing, of the state sSystem is d**(sSystem-1) (assuming sSystem>=1)
    ! we have to sum the d state corresponding to the same remaining indices
    ! indeed when we have N indices and we trace one out d state collapse into one of N-1 particles
    
    do rri=0, rightMax-1
        do lli=0, leftMax-1
            do rrj=0, rightMax-1
                do llj =0, leftMax-1
                    ! index for the traced matrix
                    iip = lli*factor+rri+1
                    jjp = llj*factor+rrj+1
                    temp = complex(0.0,0.0)
                    do val=0,dim-1
                        ! indices for the original matrix
                        ii =(lli*dim+val)*factor+rri+1
                        jj =(llj*dim+val)*factor+rrj+1
                        temp = temp+rho(ii,jj)
                    end do
                    rhop(iip,jjp)=temp
                end do
            end do
        end do
    end do
end function

subroutine print_matrix_c16(mat, nRows, nCols, realPart)
    complex*16, dimension(:,:), intent(in) :: mat
    integer*4, intent(in) :: nRows, nCols
    logical, intent(in) :: realPart
    integer*4 :: ii, jj
    do ii=1, nRows
        do jj=1,nCols-1
            if(realPart)then
                write(*, '(A,F7.4)', advance="no")'|',real(mat(ii,jj))
            else
                write(*, '(A,F7.4)', advance="no")'|',aimag(mat(ii,jj))
            end if
        end do
        if(realPart)then
            write(*, '(A,F7.4,A)')'|',real(mat(ii,jj)),'|'
        else
            write(*, '(A,F7.4,A)')'|',aimag(mat(ii,jj)),'|'
        end if
    end do
end subroutine
function magnetization_integer(x, N)result(mgn)
    integer*4, intent(in) :: x,N
    integer*4 :: ii, mgn
    mgn = 0
    do ii=0,N-1
        if(BTEST(x, ii))then
            mgn = mgn - 1
        else
            mgn = mgn + 1
        end if
    end do
end function
function isingHamiltonian(N,lmbd)result(H)
    integer*4, intent(in) :: N
    real*8, intent(in) :: lmbd
    complex*16, dimension(:,:), allocatable :: H
    integer*4 :: sz,mm,nn,mask,kk
    sz = 2**N
    allocate(H(sz,sz))
    H = 0
    do mm=0,sz-1
        H(mm+1,mm+1)=lmbd*magnetization_integer(mm, N)
        mask = 3
        do kk=1,N-1
            nn = xor(mm,mask)
            H(mm+1,nn+1)=1
            H(nn+1,mm+1)=1
            mask = mask*2
        end do
    end do
end function

function andersonHamiltonian(N, W)result(H)
    integer*4, intent(in) :: N
    real*8, intent(in) :: W
    real*8, dimension(:), allocatable :: W_i
    complex*16, dimension(:,:), allocatable :: H
    integer*4 :: mm
    allocate(H(N,N))
    allocate(W_i(N))
    call random_number(W_i)
    W_i = (2D0*W_i-1.0D0)*W
    H = 0
    do mm=2, N-1
        H(mm,mm)=W_i(mm)
        H(mm,mm+1)=1
        H(mm,mm-1)=1
    end do
    ! boundary
    H(1,2)=1
    H(1,N)=1
    H(N,N-1)=1
    H(N,1)=1
end function

end module