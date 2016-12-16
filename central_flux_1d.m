function ux = central_flux_1d(dt,dx,u,c)
nu = dt/dx;
ux = nu*u.*(circshift(c,-1)-circshift(c,1));