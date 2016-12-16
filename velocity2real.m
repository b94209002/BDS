function [u v] = velocity2real(vavg,uh,vh)

u = vavg + real(ifft2(uh));
v = vavg + real(ifft2(vh));
