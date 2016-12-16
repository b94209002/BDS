clear

exact_sol =1;
%exact_sol = 0;mx = fliplr([100 200 400 800 1600]/2); hx= 1./mx;
for m = [64 128 256 512 1024]
%m = 100; % number of grid points
h = 1/m; % Grid size
x = linspace(h/2,1-h/2,m)';
y = linspace(h/2,1-h/2,m)';
[xx yy] = meshgrid(x,y);
% x-coordinate
%xi = linspace(h/2,1-h/2,m)';% x-coordinate on the interface
a = 1; % advection coiefficient
T = 1;
mu = 0.05; % diffusion coefficient
L = 1;
nu = 0.5; 
vavg = 1;
dt = nu*h;

%KX(1,1)=1;KY(1,1)=1;
%[u v w p] = exact_solution(mu,vavg,L,xx,yy,0);
[u v] = fixed_velocity(m,m,1,1);

c = (sin(pi*xx/L).*sin(pi*yy/L)).^100; %tracer initial condition.
c0 =c ;

for t =dt:dt:T

% advection tracer 
c = Lax_wendroff_update_2d(dt,h,h,u,v,c);
%c = Fromm_update_2d(dt,h,h,u,v,c);

%c0 = (sin(pi*(xx-t)/L).*sin(pi*(yy-t)/L)).^100;

%pcolor(xx,yy,c-c0);shading flat;colorbar;title(num2str(t))
%drawnow
end

%eval(['w' num2str(log2(m/64)+1) '= c;'])
eval(['w_err(' num2str(log2(m/64)+1) ')= compute_error(h^2,c,c0);']);

end


%compute_norm(wnorm1,wnorm2,wnorm3,wnorm4,wnorm5)

return
