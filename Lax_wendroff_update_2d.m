function c = Lax_wendroff_update_2d(dt,dx,dy,u,v,c)

c = Lax_wendroff_update_1d(dt,dx,u',c')';
c = Lax_wendroff_update_1d(dt,dy,v,c);

