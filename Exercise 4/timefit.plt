fname = system("cat file_to_fit.txt")
f(x)=a+b*x+c*x**2+d*x**3
fit f(x) fname via a,b,c,d
set print fname.'_fit_params.txt'
print "#f(x)=a+b*x+c*x**2+d*x**3"
print a,a_err
print b,b_err
print c,c_err
print d,d_err
