program numprecision
  integer*2 :: num1s,num2s
  integer*4 :: num1,num2
  real*4:: real1f,real2f
  real*8 :: real1d,real2d
  num1s = 2000000
  num2s = 1
  num1 = 2000000
  num2 = 1
  real1f = acos(-1e0)*(1e32)
  real2f  = sqrt(2e0)*(1e21)
  real1d = acos(-1d0)*(1d32)
  real2d  = sqrt(2d0)*(1d21)
  print*,"Summing 2.000.000 and 1 using integer2",num1s+num2s
  print*,"Summing 2.000.000 and 1 using integer4",num1+num2
  print*,"Summing pi*10^32 and sqrt(2)*10^21 using real4:",real1f,"+",real2f,"=",real1f+real2f
  print*,"Summing pi*10^32 and sqrt(2)*10^21 using real8",real1d,"+",real2d,"=",real1d+real2d
  stop
  end program
