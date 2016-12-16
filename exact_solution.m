function [u v w p] = taylor_vortex(mu,vavg,L,xx,yy,t)
u_exact =@(x,y,t) vavg -2*exp(-8*pi^2*mu*t/L^2).*cos(2*pi*(x-vavg*t)/L).*sin(2*pi*(y-vavg*t)/L);
v_exact =@(x,y,t) vavg +2*exp(-8*pi^2*mu*t/L^2).*sin(2*pi*(x-vavg*t)/L).*cos(2*pi*(y-vavg*t)/L); 
w_exact =@(x,y,t) 8*pi/L*exp(-8*pi^2*mu*t/L^2).*cos(2*pi*(x-vavg*t)/L).*cos(2*pi*(y-vavg*t)/L);
p_exact= @(x,y,t) -exp(-16*pi^2*mu*t/L^2).*(cos(4*pi*(x-vavg*t)/L)+cos(4*pi*(y-vavg*t)/L));

u = u_exact(xx,yy,t);
v = v_exact(xx,yy,t);
w = w_exact(xx,yy,t);
p = p_exact(xx,yy,t);