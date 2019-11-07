fname = 'rmt_pdf.txt'
f(x)=a*(x**b)*exp(-c*(x**d))
fit f(x) fname via a,b,c,d
set print fname.'_fit_params.txt'
print "#f(x)=a*x**b*exp(-c*x)"
print a,a_err
print b,b_err
print c,c_err
print d,d_err
