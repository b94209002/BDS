function ux = Fromm_flux_1d(dt,dx,u,c)

nu = dt/dx;
% first order term
ux = u.*c;
% second order term
ux = ux + .25*central_flux_1d(dx-dt,dx,u,c);
ux = nu*ux;