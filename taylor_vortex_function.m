function [ufun vfun pxfun pyfun] = taylor_vortex_function(mu,vavg,L)
ufun =@(x,y,t) vavg - 2*exp(-8*pi^2*mu*t/L^2).*(cos(2*pi*(x-vavg*t)/L).*sin(2*pi*(y-vavg*t)/L))';
vfun =@(x,y,t) vavg + 2*exp(-8*pi^2*mu*t/L^2).*(sin(2*pi*(x-vavg*t)/L).*cos(2*pi*(y-vavg*t)/L))'; 
pxfun = @(x,y,t) 4*pi/L*exp(-16*pi^2*mu*t/L^2).*sin(4*pi*(x-vavg*t)/L)';
pyfun = @(x,y,t) 4*pi/L*exp(-16*pi^2*mu*t/L^2).*sin(4*pi*(y-vavg*t)/L)';
%pxfun = @(x,y,t) -exp(-16*pi^2*mu*t/L^2).*cos(4*pi*(x-vavg*t)/L)';
%pyfun = @(x,y,t) -exp(-16*pi^2*mu*t/L^2).*cos(4*pi*(y-vavg*t)/L)';