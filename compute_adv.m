function vwh = compute_adv(KX,KY,dealias,u,v,wh)
[whx why] = gradient_omega(KX,KY,wh);
[wx wy] = gradient2real(whx,why);
vw =advection(u,v,wx,wy);
vwh = fft2(vw); vwh = vwh.*dealias;
