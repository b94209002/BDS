function [u v] = fixed_velocity2(xx,yy,u0,v0)

%u = u0*cos(2*pi*(xx-yy+.5));
%v = v0*cos(2*pi*(xx-yy+.5));

u = u0*cos(2*pi*(yy+.5));
v = v0*cos(2*pi*(yy+.5));