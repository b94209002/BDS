function c = Lax_wendroff_update_1d(dt,dx,u,c)

cx = Lax_wendroff_limiter_flux_1d(dt,dx,u,c);
%cx = Lax_wendroff_flux_1d(dt,dx,u,c);
c = c - (circshift(cx,0) - circshift(cx,1));