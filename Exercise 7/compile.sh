gfortran -c  -I/usr/local/include -L/usr/lib/ -llapack -lfftw3 -lm qiutil.f03 
gfortran -c  shiftedho.f03
gfortran qiutil.o shiftedho.o -o shiftedho.out
