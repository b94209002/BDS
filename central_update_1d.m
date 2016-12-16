function c= central_update_1d(dt,dx,u,c)

ucx = central_flux_1d(dt,dx,u,c);

c = c - .5*ucx;


