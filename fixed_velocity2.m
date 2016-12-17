function [u v] = fixed_velocity2(xf,yf,xx,yy,u0,v0)

u = u0*cos(2*pi*(xf-yy))';
v = v0*cos(2*pi*(xx-yf))';

%u = u0*cos(2*pi*(yy+.5))';
%v = v0*cos(2*pi*(xx+.5))';