function w_exact = advdiff_exact_solution(m,t,a,d,u0)
%exact solution for fix coefficient
ik = [0:m/2-1 0 -m/2+1:-1]';
u0h = fft(u0);
A_ediff = (2*pi*1i*ik.*(2*pi*1i*ik*d));
A_exact = (2*pi*1i*ik.*(-a));
wh = u0h.*exp(t*(A_exact+A_ediff));
w_exact = ifft(wh);