function c= update_tracer_2d(dt,dx,dy,u,v,c)

c= update_tracer_1d(dt,dx,u,c);
c= update_tracer_1d(dt,dy,v',c')';


