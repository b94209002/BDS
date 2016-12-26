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
k = (2*pi*1i)*[0:(m/2-1) (-m/2):(-1)]; % Vector of wavenumbers
[KX KY]  = meshgrid(k,k); % Matrix of (x,y) wavenumbers corresponding
                          % to Fourier mode (m,n)
delsq = (KX.^2 + KY.^2);

[u v] = fixed_velocity(m,m,a,a);

%m2 = m*m;e = ones(m,1);I = speye(m);
%d0 = spdiags([e e -2*e e e],[-m+1 -1:1 m-1],m,m);
%D = mu*dt/h/h*(kron(I,d0)+kron(d0,I));
%CNf = speye(m2) +.5*D;
%CNb = speye(m2) -.5*D;

D = mu*delsq;
CNf = ones(m) +.5*mu*dt*delsq;
CNb = ones(m) -.5*mu*dt*delsq;

%c = sin(2*pi*xx'/L); %tracer initial condition.
c = sin(2*pi*(xx+yy)'/L); %tracer initial condition.
ch = fft2(c);


for t =dt:dt:T

%S = reshape(D*reshape(c,m2,1),m,m);
%F = BDS_update_2d(dt,h,h,u,v,S,c);
%c = Crank_Nicolson(m,m2,CNb,CNf,c,F);

S = real(ifft2(D.*ch));
F = BDS_update_2d(dt,h,h,u,v,S,c);
Fh = fft2(F);
ch = Crank_Nicolson(CNb,CNf,ch,Fh);
c = real(ifft2(ch));

pcolor(xx,yy,c');shading flat;colorbar;%hold on;quiver(xx,yy,u',v','w');hold off;
title([' t = ' num2str(t) ' max = ' num2str(max(max(c))) ', min = ' num2str(min(min(c))) ] );drawnow
end
%c = real(ifft2(ch));
%c0 = sin(2*pi*(xx)'/L)*exp(-4*pi^2*mu);
c0 = sin(2*pi*(xx+yy)'/L)*exp(-8*pi^2*mu);
%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h^2,c,c0);']);

eval(['w_max(' num2str(log2(m/64)+1) ')= max(max(c));']);
eval(['w_min(' num2str(log2(m/64)+1) ')= min(min(c));']);


eval(['w' num2str(log2(m/64)+1) '= c;']);
eval(['x' num2str(log2(m/64)+1) '= xx;']);
end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
