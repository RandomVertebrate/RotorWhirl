clear all
syms x k h theta C1 C2 C3

C1 = 1000
C2 = 1e8
C3 = 1.7

r = 0.0375
n = 6

phi0 = pi/6

phi = 2*pi/n

f(x) = C1/(1+C2*x^(C3))

F(h,theta) = symsum(f(h+theta*r*sin(k*phi+phi0)),k,0,n-1)

M(h, theta) = symsum(r*sin(k*phi+phi0)*f(h+theta*r*sin(k*phi+phi0)),k,0,n-1)

latex(F)

latex(M)