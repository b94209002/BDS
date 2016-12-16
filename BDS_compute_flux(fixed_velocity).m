function [usl usr vsl vsh] = BDS_compute_flux(dt,dx,dy,u,v,sxy,sx,sy,sh)
% Positive velocity up, vp at interface and minous velocity um,vm
up = .5*(u+abs(u)); um = .5*(u-abs(u));
vp = .5*(v+abs(v)); vm = .5*(v-abs(v));

% corner index dx/2 or -dx/2 us, vs
us = ones(size(u)); vs = ones(size(v));
us(u<0) = -1; vs(v<0) = -1;

[fxy,fx,fy,fh] = BDS_compute_flux_source(u,v,sxy,sx,sy,sh);

%note here suppose u> 0
sxp = .5*sx.*(dx-u*dt) + sh;
sxm = .5*circshift(sx,[-1,0]).*(-dx-u*dt) + circshift(sh,[-1,0]);
usr = up.*sxp + um.*sxm;
usl = circshift(usr,[1 0]);

% DEF area in BDS paper
Gammay = fh + 1/6*(fx.*(3*us*dx-4*dt*u) + fy.*(3*vs*dy-2*dt*v));
Gammay = Gammay + 1/12*fxy.*(3*dx*dy*us.*vs - 4*dt*dy*vs.*u - 2*dt*dx*us.*v + 3*dt*dt*u.*v);

usly = .5*dt*dt*u.*v.*Gammay;

usr = usr - 1/dt/dy*(circshift(usly,[0 0]) - circshift(usly,[0 1]));
usl = usl - 1/dt/dy*(circshift(usly,[1 0]) - circshift(usly,[1 1]));

%note here suppose v > 0 
syp = .5*sy.*(dy-v*dt) + sh;
sym = .5*circshift(sy,[0,-1]).*(-dy-v*dt) + circshift(sh,[0,-1]);
vsh = vp.*syp+vm.*sym;
vsl = circshift(vsh,[0 1]);

Gammax = fh + 1/6*(fy.*(3*vs*dy-4*dt*v) + fx.*(3*us*dx-2*dt*u));
Gammax = Gammax + 1/12*fxy.*(3*us.*vs*dx*dy - 4*dt*dx*us.*v - 2*dt*dy*vs.*u + 3*dt*dt*u.*v);

vshx = .5*dt*dt*u.*v.*Gammax;

vsh = vsh - 1/dt/dx*(circshift(vshx,[0 0]) - circshift(vshx,[1 0]));
vsl = vsl - 1/dt/dx*(circshift(vshx,[0 1]) - circshift(vshx,[1 1]));

% Gamma index for u 
% u > 0, v > 0 
%GammaP = sh+ 1/12*sxy.*(3*dx*dy - 4*dt*dy*u -2*dt*dx*v+3*dt*dt*u.*v);
%GammaP = GammaP + 1/6*(sx.*(3*dx-4*dt*u) +sy.*(3*dy-2*dt*v)); 
% u < 0, v > 0 
%GammaM = circshift(sh,[-1,0])+ 1/12*circshift(sxy,[-1,0]).*(-3*dx*dy - 4*dt*dy*u + 2*dt*dx*v + 3*dt*dt*u.*v);
%GammaM = GammaM + 1/6*(circshift(sx,[-1 0]).*(-3*dx-4*dt*u) +circshift(sy,[-1,0]).*(3*dy-2*dt*v)); 
% u > 0, v < 0 
%GammaP = circshift(sh,[0,-1])+ 1/12*circshift(sxy,[0,-1]).*(-3*dx*dy + 4*dt*dy*u - 2*dt*dx*v + 3*dt*dt*u.*v);
%GammaP = GammaP + 1/6*(circshift(sx,[0,-1]).*(3*dx-4*dt*u)+ circshift(sy,[0,-1]).*(-3*dy-2*dt*v)); 
% u < 0, v < 0 
%GammaM = circshift(sh,[0,-1])+ 1/12*circshift(sxy,[0,-1]).*(3*dx*dy + 4*dt*dy*u + 2*dt*dx*v + 3*dt*dt*u.*v);
%GammaM = GammaM + 1/6*(circshift(sx,[0,-1]).*(-3*dx-4*dt*u) +circshift(sy,[0,-1]).*(-3*dy-2*dt*v)); 

