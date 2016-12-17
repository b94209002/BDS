function [u v] = taylor_vortex(ufun,vfun,xf,yf,xx,yy,t)

u = ufun(xf,yy,t);
v = vfun(xx,yf,t);
