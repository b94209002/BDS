function c = Fromm_update_2d(dt,dx,dy,u,v,c)

c = Fromm_update_1d(dt,dx,u',c')';
c = Fromm_update_1d(dt,dy,v,c);

