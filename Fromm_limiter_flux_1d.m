function cx = Fromm_limiter_flux_1d(dt,dx,u,c)

nu = dt/dx;
% first order term
up = zeros(size(u)); um = up; us = um;
us(u >= 0) = 1; us(u < 0 ) = -1;
up(u >= 0) = 1; up = up.*u;
um(u < 0) = 1; um=um.*u;
ux = up.*c + um.*circshift(c,-1);
% second order term
s = compute_slope(up,um,c);
% Fromm's slope  
s1 = .5*(1+s);
s2 = min(2,2*s);
%MC limiter for Fromm's scheme
s =max(0,min(s1,s2));
nu2 = .5*us.*(1-us.*u*dt/dx);
cx = ux + s.*nu2.*upwind_flux_1d(1,1,u,circshift(c,-1));
cx = nu*cx;