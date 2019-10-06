program matrixvec
  real, dimension(1000,1000) :: a
  real, dimension(1000) :: x,y
  real :: start, finish
  integer :: i,j
  N = 1000
  do j=1,N
     do i=1,N
        a(i,j)=RAND(0)
     end do
     x(j)=RAND(0)
     y(j)=0
  end do
  print*,"Starting..."
  call cpu_time(start)
  do j=1,N
     do i=1,N
        y(i)=y(i)+a(i,j)*x(j)
     end do
  end do
  call cpu_time(finish)
  print*,finish-start
   call cpu_time(start)
  do i=1,N
     do j=1,N
        y(i)=y(i)+a(i,j)*x(j)
     end do
  end do
  call cpu_time(finish)
  print*,finish-start
  end program
