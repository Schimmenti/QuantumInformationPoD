program numprecision
  integer*2 :: a,b
  integer*4 :: c,d
  real*4:: e,f
  real*8 :: g,h
  a = 2000000
  b = 1
  c = 2000000
  d = 1
  e = 3.14159265359E32
  f  = 1.41421356237E21
  g = 3.14159265359E32
  h  = 1.41421356237E21
  print*,"Summing 2.000.000 and 1 using integer2",a+b
  print*,"Summing 2.000.000 and 1 using integer4",c+d
  print*,"Summing pi*10^32 and sqrt(2)*10^21 using real4",e+f
  print*,"Summing pi*10^32 and sqrt(2)*10^21 using real8",g+h
  stop
  end program
