program test
use nbody
implicit none
type(pstate), dimension(2) :: states
type(pstate) :: testState
complex*16, dimension(:,:), allocatable :: rho1, rho2, rho, testRho
complex*16 :: trace
real*8 :: cBfrX, cBfrY, start, finish
integer*4 :: ii, mode,Nmax
character(len=32) :: arg
mode=0
Nmax=15
open(unit=42, file='timing.txt')
do ii=2,Nmax
    testState = random_pure_state(2, ii, .FALSE.)
    write(42, '(I0,A)', advance='no')ii,' '
    call cpu_time(start)
    testRho = vector_oproduct(testState%psi, testState%psi)
    call cpu_time(finish)
    deallocate(testState%psi)

    write(42, '(G0)')(finish-start)
    
    deallocate(testRho)
end do
close(42)
if(COMMAND_ARGUMENT_COUNT()>=2)then
    call GET_COMMAND_ARGUMENT(1, arg)
    open(unit=42,file=arg)
    call get_command_argument(2, arg)
    read(arg,*)mode
    if(mode.EQ.0)then
        print*,'Reading two separable states...'
        do ii=1,2
            states(ii)=allocate_pure(2, 1)
        end do
        do ii=0,3
            read(42,*)cBfrX
            read(42,*)cBfrY
            states(ii/2+1)%psi(mod(ii,2)+1)=complex(cBfrX, cBfrY)
        end do
        
        print*,"State for system A",states(1)%psi
        print*,"State for system B",states(2)%psi
        rho1= vector_oproduct(states(1)%psi,states(1)%psi )
        rho2= vector_oproduct(states(2)%psi,states(2)%psi )
        rho = tensor_product(rho2, rho1, 2, 2)
    else
        print*,'Reading qubits density matrix for two particles...'
        allocate(rho(4,4))
        ! read the real part
        do ii=1,4
            read(42, *)rho(ii,1:4)
        end do
    end if
    close(42)
else
    
    
end if
trace = trace_c16(rho, 4)
if((abs(real(trace)-1)>=1D-10).OR.(abs(dimag(trace))>=1D-10))then
    print*,"Fatal error, stopping..."
    stop
end if
rho1 = trace_dm(rho, 2, 2, 2)
rho2 = trace_dm(rho, 2, 2, 1)
print*,'-----------------------------'
print*,'Real part of the density matrix A+B'
call print_matrix_c16(rho, 4, 4, .TRUE.)
print*,'Imaginary part of the density matrix A+B'
call print_matrix_c16(rho, 4, 4, .FALSE.)
print*,'-----------------------------'
print*,'Real part of the reduced density matrix A'
call print_matrix_c16(rho1, 2, 2, .TRUE.)
print*,'Imaginary part of the reduced density matrix A'
call print_matrix_c16(rho1, 2, 2, .FALSE.)
print*,'-----------------------------'
print*,'Real part of the reduced density matrix B'
call print_matrix_c16(rho2, 2, 2, .TRUE.)
print*,'Imaginary part of the reduced density matrix B'
call print_matrix_c16(rho2, 2, 2, .FALSE.)
print*,'-----------------------------'


end program