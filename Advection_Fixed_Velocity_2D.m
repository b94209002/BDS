clear

exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;
for m = [64 128 256 512]
m
%m = 100; % number of grid points
h = 1/m; % Grid size
x = linspace(h/2,1-h/2,m)';y = linspace(h/2,1-h/2,m)';
[xx yy] = meshgrid(x,y);
x = linspace(0,1-h,m)';y = linspace(0,1-h,m)';
[xf yf] = meshgrid(x,y);

% x-coordinate
%xi = linspace(h/2,1-h/2,m)';% x-coordinate on the interface
a = 1; % advection coiefficient
T = 1;
mu = 0.01; % diffusion coefficient
vavg =0;
L = 1;
nu = 0.5; 
dt = nu*h;

% Use FFT approach diffusion
k = (2*pi)*[0:(m/2-1) (-m/2):(-1)]; % Vector of wavenumbers
[KX KY]  = meshgrid(k,k); % Matrix of (x,y) wavenumbers corresponding
                          % to Fourier mode (m,n)
delsq = -(KX.^2 + KY.^2);

% u(1,1) is the velocity at (0,0)
%[u v] = fixed_velocity(m,m,2*a,-a);
[u v] = fixed_velocity2(xf,yf,xx,yy,a,a);
[ufun vfun pxfun pyfun] = taylor_vortex_function(mu,vavg,1);


D = mu*dt*delsq;
CNf = ones(m) +.5*mu*dt*delsq;
CNb = 1./(ones(m) -.5*mu*dt*delsq);

[c c2]=taylor_vortex(ufun,vfun,xx,yy,xx,yy,0);
ch = fft2(c);
%c = (sin(pi*xx/L).*sin(pi*yy/L)).^100; %tracer initial condition.
%ch = fft2(c);
%c = ones(m);
%c((xx' < .6 & xx' >.4)) = 1;
%[uc vc] = fixed_velocity2(xx,yy,xx,yy,a,a);
%c0 =c;
%c0 =((sin(pi*(xx-uc')/L).*sin(pi*(yy-vc')/L)).^100)';
for t =dt:dt:T

 [u v] = taylor_vortex(ufun,vfun,xf,yf,xx,yy,t);
% 
% % advection tracer 
% %c1 = Lax_wendroff_update_2d(dt,h,h,circshift(u,[-1 0]),circshift(v,[0 -1]),c1);
% %c1 = Fromm_update_2d(dt,h,h,u,v,c1);
px = -pxfun(xx,yy,t+.25*dt);
c = c + .5*dt*px;
ch = fft2(c);
S = ifft2(D.*ch);
%px = -(pxfun(xx,yy,t+.5*dt));
F = BDS_update_2d(dt,h,h,u,v,S,c);
Fh = fft2(F);

ch = Crank_Nicolson(CNb,CNf,ch,Fh);
%c2 = BDS_update_2d(dt,h,h,u,v,c2);
%c0 = (sin(pi*(xx-t)/L).*sin(pi*(yy-t)/L)).^100;
c = real(ifft2(ch));
px = -pxfun(xx,yy,t+.75*dt);
c = c + .5*dt*px;
pcolor(xx,yy,c');shading flat;colorbar;%hold on;quiver(xx,yy,u',v','w');hold off;
title([' t = ' num2str(t) ' max = ' num2str(max(max(c))) ', min = ' num2str(min(min(c))) ] );drawnow
end
c = real(ifft2(ch));
[c0 vc] = taylor_vortex(ufun,vfun,xx,yy,xx,yy,t);
%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h^2,c,c0);']);

end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
