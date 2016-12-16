function [ch vch] = update_tracer(dt,KX,KY,filter,u,v,ch,vch)
vch0 = -.5*vch;
[chx chy] = gradient_omega(KX,KY,ch);
[cx cy] =gradient2real(chx,chy);
vc = advection(u,v,cx,cy);
vch = fft2(vc); vch = vch.*filter;
ch = ch - dt*(1.5*vch +vch0);