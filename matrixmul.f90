program matrixmul
  integer N, stat
  integer i,j,k
  real :: temp
  real :: start, finish
  real, dimension (100,100) ::  a, b, c,d
  N=100
  print*,"Matrix size is:",N
  
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
! D^T_{i,j}=\sum_k B^T_{i,k} A^T_{k,j}=\sum_k B_{k,i} A_{j,k}
  call cpu_time(start)
  do i=1,N
     do j=1,N
        do k=1,N
           d(i,k) = d(i,k) + a(i,j)*b(j,k)
        end do
     end do
  end do
  call cpu_time(finish)
  print*, "Time elapsed for row-fast multiplication is ", finish-start
  
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
  do  j=1,N
     do i = 1,N
        if (c(i,j) /= d(i,j)) then
           print*,"ERROR"
           stop
        end if
     end do
  end do
  end program
