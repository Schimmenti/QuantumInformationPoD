program main
use heisenberg
implicit none
! variables
real*8, dimension(:,:), allocatable :: H0,H
integer*4, dimension(:), allocatable :: magnetizations, sector
integer*4 :: L, info, idx, jdx, kdx, nSamples,dummy
real*8 :: field_strength
real*8, dimension(:), allocatable :: rnd_field
integer*4, dimension(:,:), allocatable :: spin_map, sector_map
complex*16, dimension(:), allocatable :: wave_cfs
!parameters
integer*4, dimension(3) :: Ls, sampleSize
integer*4 :: lwork, dim, field_verbosity
real*8, dimension(7) :: fields
real*8, dimension(:), allocatable :: work
real*8, dimension(:), allocatable :: rwork
!results
real*8, dimension(:), allocatable :: spectrum, temp_seq, temp_seq2, min_d, max_d
real*8, dimension(:,:), allocatable :: prob
real*8, dimension(:,:), allocatable :: f_coef
real*8, dimension(:,:), allocatable :: r_coef
real*8, dimension(:,:), allocatable :: dmagn,local_magn
real*8 :: math_pi
math_pi = 4.D0*datan(1.D0)
open(unit=42, file='heisenberg_results.txt')
Ls = [6,8,10]
! Big but reasonable sizes should be around 5000 for L <= 8, for L=10 around 1000
sampleSize = [10000,1000,100]
fields = [0.6D0, 1.0D0, 2.0D0, 2.7D0, 3.6D0, 5.0D0, 8.0D0]
! verbosity for printing percentage of completion at given L
field_verbosity = 1

! allocate results
allocate(f_coef(size(Ls), size(fields)))
allocate(r_coef(size(Ls), size(fields)))
allocate(dmagn(size(Ls), size(fields)))
f_coef = 0D0
r_coef = 0D0
dmagn = 0D0

! loop over system sizes
do idx=1,size(Ls)
    ! #####################################################
    ! ################  PRELIMINARIES  ####################
    ! #####################################################


    ! current system size
    L = Ls(idx)
    ! current number of samples to generate
    nSamples = sampleSize(idx)
    
    print*,"Probing size",L,"with",nSamples,"samples"

    ! initialize variables and Hamiltonian (non disordered)
    ! done just once per system size
    spin_map = get_spin_map(L)
    magnetizations = create_magnetizations(L)
    sector = get_sector(0,magnetizations,L)
    sector_map = spin_map(sector+1, 1:L)
    allocate(wave_cfs(L))
    do dummy=1,L
        wave_cfs(dummy) = exp(complex(0,1)*2D0*math_pi*(dummy-1)/L)
    end do
    ! sector size
    dim = size(sector)
    ! compute non disordered Hamiltonian (PBC)
    H0 = hamiltonian_heisenberg(sector,L,1.0D0,1.0D0)
    ! allocate the disordered one and the random field
    allocate(H(dim, dim))
    allocate(rnd_field(L))

    ! prepare diagoalization parameters for current system size    
    ! optimal lwork
    allocate(spectrum(dim))
    lwork=-1
    allocate(work(1))
    allocate(rwork(max(1, 3*dim-2)))
    call dsyev('V','U',dim,H0,dim,spectrum,work,lwork,rwork,info)
    lwork = int(real(work(1)))
    deallocate(work)
    deallocate(rwork)


    ! #####################################################
    ! ################ DIAGONALIZATION ####################
    ! #####################################################

    
    
    ! loop over field strength
    do jdx=1,size(fields)
        ! print status
        if((field_verbosity > 0).and.(mod(jdx-1,field_verbosity).eq.0))then
            print*,"Size",L,"Percentage",real(100*(jdx-1))/size(fields)
        end if


        ! current field strength
        field_strength = fields(jdx)

        ! loop over disorder realizations
        do kdx=1,nSamples
            H = H0 ! copy the base Hamiltonian, we overwrite the old realization
            call random_number(rnd_field)
            ! scale the field
            rnd_field = rnd_field*field_strength
            ! takes H and add in the diagonal the local field
            call add_zfield(H, sector, L, rnd_field)
            
            ! diagonalize -> hamiltonian is real and symmetric, we can use dsyev
            allocate(work(max(1,lwork)))
            allocate(rwork(max(1, 3*dim-2)))
            call dsyev('V','U',dim,H,dim,spectrum,work,lwork,rwork,info)
            if(info.ne.0)then
                print*,'Error in diagonalization. Current size', L, 'Current field strength',field_strength
            end if
            deallocate(work)
            deallocate(rwork)

            ! now inside H we have the eigenvectors (real!) while spectrum contains eigenvalues

            ! #####################################################
            ! ################    RESULTS      ####################
            ! #####################################################

            ! prob(s,n) = |v(s,n)|^2
            prob = abs(H)**2
            ! LOCAL MAGNETIZATION DELTA
            local_magn = matmul(transpose(sector_map), prob)
            allocate(temp_seq2(dim))
            do dummy=1,dim
                temp_seq2(dummy) = sum(abs(differentiate(local_magn(1:L,dummy))))/(L-1)
            end do
            dmagn(idx,jdx) = dmagn(idx,jdx) + sum(temp_seq2)/dim
            deallocate(temp_seq2)
            ! WAVE RELAXATION
            f_coef(idx,jdx) = f_coef(idx,jdx) + sum(real(m_wave_compute(prob, sector_map, L, dim, wave_cfs)))/dim
            ! ENERGY SPACING DISTRIBUTION
            temp_seq = abs(differentiate(spectrum))
            min_d =  min(temp_seq(1:dim-2),temp_seq(2:dim-1))
            max_d =  max(temp_seq(1:dim-2),temp_seq(2:dim-1))
            r_coef(idx,jdx) = r_coef(idx,jdx) + sum(min_d/max_d)/dim
        end do
        ! divide by number of samples to get the average
        f_coef(idx,jdx) = f_coef(idx,jdx)/nSamples
        r_coef(idx,jdx) = r_coef(idx,jdx)/nSamples
        dmagn(idx,jdx) = dmagn(idx,jdx)/nSamples
        write(42,*)L,field_strength,f_coef(idx,jdx),r_coef(idx,jdx), dmagn(idx,jdx)
    end do

    ! deallocate everything since we are changing system size
    deallocate(H)
    deallocate(H0)
    deallocate(rnd_field)
    deallocate(magnetizations)
    deallocate(sector)
    deallocate(spectrum)
    deallocate(wave_cfs)
    deallocate(prob)
    deallocate(temp_seq, min_d, max_d, local_magn)
end do
close(42)
end program