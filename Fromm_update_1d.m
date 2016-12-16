function c = Fromm_update_1d(dt,dx,u,c)

%cx = Fromm_flux_1d(dt,dx,u,c);
cx = Fromm_limiter_flux_1d(dt,dx,u,c);
c = c - (circshift(cx,0) -circshift(cx,1));

