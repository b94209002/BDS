function F = BDS_update_2d(dt,dx,dy,u,v,f,c)
%function c = BDS_update_2d(dt,dx,dy,u,v,f,c)


[sxy sx sy sh] = BDS_bilinear_poly(dx,dy,c);

%f0 = zeros(size(f));
[usl usr vsl vsh] = BDS_compute_flux(dt,dx,dy,u,v,sxy,sx,sy,sh,f);

ua = circshift(u,[-1 0]);va = circshift(v,[0 -1]);
%c = c - dt/dx*(ua.*usr-u.*usl) - dt/dy*(va.*vsh - v.*vsl) + f*dt;
F =  - dt/dx*(ua.*usr-u.*usl) - dt/dy*(va.*vsh - v.*vsl) ;
return