function ux = upwind_flux_1d(dt,dx,u,c)
nu =dt/dx.*u;
ux = nu.*(c-circshift(c,1));