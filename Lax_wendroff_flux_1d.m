function ux = Lax_wendroff_flux_1d(dt,dx,u,c)

nu = dt/dx;
% first order term
ux = u.*c;
% second order term
ux = ux + upwind_flux_1d(.5*(dx-dt),dx,u,circshift(c,-1));
ux = nu*ux;