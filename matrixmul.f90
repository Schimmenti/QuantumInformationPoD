program matrixmul
  integer N, stat
  integer i,j,k
  real :: temp
  real :: start, finish
  real, dimension (:,:), allocatable ::  a, b, c,d, e
  N=100
  print*,"Matrix size is:",N
! we allocate the matrices with the given size; this is done in order
! to allow dynamic sizing for N scaling testing
  allocate(a(N,N))
  allocate(b(N,N))
  allocate(c(N,N))
  allocate(d(N,N))
  allocate(e(N,N))
! populate the array using random [0,1] values
  do  j=1,N
     do i = 1,N
        a(i,j) = RAND(0)
        b(i,j) = RAND(0)
        c(i,j)=0
        d(i,j)=0
     end do
  end do
  print*,"Starting..."
  
  call cpu_time(start)
  do i=1,N
     do j=1,N
        do k=1,N
           d(i,k) = d(i,k) + a(i,j)*b(j,k)
        end do
     end do
  end do
  call cpu_time(finish)
  print*, "Time elapsed for row-slow multiplication is ", finish-start
  
  call cpu_time(start)
  do i=1,N
     do j=1,N
        do k=1,N
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
  call cpu_time(finish)
  print*, "Time elapsed for column-fast multiplication is ", finish-start
! checking whether the two resulting matrices are the same
  do j=1,N
     do i = 1,N
        if (c(i,j) /= d(i,j)) then
           print*,"ERROR"
           stop
        end if
     end do
  end do
  call cpu_time(start)
  e = matmul(a,b)
  call cpu_time(finish)
  print*, "Time elapsed for intrinsic multiplication is ", finish-start
! goodbye
  deallocate(a)
  deallocate(b)
  deallocate(c)
  deallocate(d)
  deallocate(e)
  end program
