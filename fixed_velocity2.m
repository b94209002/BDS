function [u v] = fixed_velocity2(xx,yy,u0,v0)

u = -2*u0*cos(2*pi*xx).*sin(2*pi*yy);
v = 2*v0*sin(2*pi*xx).*cos(2*pi*yy);