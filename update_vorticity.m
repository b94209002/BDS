function [vwh wh u v] = update_vorticity(m,dt,mu,vavg,KX,KY,delsq,poisson,dealias,wh,vwh)

vwh0 = vwh;
psi = compute_stream_function(wh,poisson);
[uh vh] = gradient_psi(KX,KY,psi);
[whx why] = gradient_omega(KX,KY,wh);
[u v] =velocity2real(vavg,uh,vh);
[wx wy] = gradient2real(whx,why);
vw =advection(u,v,wx,wy);
vwh = fft2(vw); vwh = vwh.*dealias;
ddw = forward_diffusion(m,dt,mu,delsq,wh);
adv = Adams_Bashforth(vwh,vwh0);
wh = backward_diffusion(m,dt,mu,delsq,ddw,adv);