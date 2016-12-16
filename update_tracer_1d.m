function c= update_tracer_1d(dt,dx,u,c)

ucx = upwind_flux_1d(dt,dx,u,c);

c= c - (ucx);


