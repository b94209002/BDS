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
mu = 0.00; % diffusion coefficient
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
%[u v] = fixed_velocity(m,m,a,a);
[u v] = fixed_velocity2(xf,yf,xx,yy,a,a);
%[ufun vfun pxfun pyfun] = taylor_vortex_function(mu,vavg,1);

D = mu*dt*delsq;
CNf = ones(m) +.5*mu*dt*delsq;
CNb = 1./(ones(m) -.5*mu*dt*delsq);

%[c c2]=taylor_vortex(ufun,vfun,xx,yy,xx,yy,0);
%ch = fft2(c);
c = (sin(pi*xx/L).*sin(pi*yy/L)).^100; %tracer initial condition.
%c = sin(2*pi*yy/L) ; 
%ch = fft2(c);
%c = ones(m);
%c((xx' < .6 & xx' >.4)) = 1;
[uc vc] = fixed_velocity2(xx,yy,xx,yy,a,a);
%[uc vc] = fixed_velocity(m,m,a,a);
%c0 =c;
c0 =((sin(pi*(xx-uc')/L).*sin(pi*(yy-vc')/L)).^100)';
for t =dt:dt:T

S = 0;
%BDS update
%c = c+BDS_update_2d(dt,h,h,u,v,S,c);
%Fromm's method
ua = circshift(u,[-1,0]);
va = circshift(v,[0,-1]);
%c = Fromm_update_2d(dt,h,h,u,v,c);

%Lax-Wrendroff method
c = Lax_wendroff_update_2d(dt,h,h,ua,va,c);

%MPDATA update
%c = MPDATA_update_2D(dt,h,h,u,v,S,c);
%c0 = sin(2*pi*(yy-a*t)/L) ; 
pcolor(xx,yy,(c)');shading flat;colorbar;%hold on;quiver(xx,yy,u',v','w');hold off;
title([' t = ' num2str(t) ' max = ' num2str(max(max(c))) ', min = ' num2str(min(min(c))) ] );drawnow
end
%c = real(ifft2(ch));
%[c0 vc] = taylor_vortex(ufun,vfun,xx,yy,xx,yy,t);
%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h^2,c,c0);']);

eval(['w_max(' num2str(log2(m/64)+1) ')= max(max(c));']);
eval(['w_min(' num2str(log2(m/64)+1) ')= min(min(c));']);

eval(['w' num2str(log2(m/64)+1) '= c;']);
eval(['x' num2str(log2(m/64)+1) '= xx;']);
end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
